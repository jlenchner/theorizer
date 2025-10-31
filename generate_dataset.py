"""
SynPAT dataset generator: builds random, dimensionally consistent axiom systems,
optionally generates replacement axioms, computes algebraic consequences (via M2),
and emits noiseless/noisy datasets.

Notes
-----
- Timeout uses SIGALRM (POSIX only). On Windows, set --timeout None or run on macOS/Linux.
"""

import argparse
import ast
import io
import itertools
import os
import random
import re
import shutil
import signal
import time
from contextlib import contextmanager, redirect_stdout
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import psutil
from sympy import Integer

# Project imports
from consequence_data_generation import (
    run_consequence_noisy_data_generation,
    run_consequence_noiseless_data_generation,
)
from consequence_generation import run_consequence_generation
from constant import *
from derivative import *
from equation import *
from equationSystem import *
from m2_functions import *
from system_data_generation import (
    run_noiseless_system_data_generation,
    run_noisy_system_data_generation,
)
from theorizer import *
from unitOfMeasure import *

# ==============================
# Exceptions & Timeout Context
# ==============================
class TimeoutException(Exception):
    """Raised when `time_limit` context times out."""
    pass

@contextmanager
def time_limit(seconds: int):
    """
    POSIX-only alarm-based timeout context.

    Args:
        seconds: Maximum wall time for the wrapped block. Use None to disable externally.

    Raises:
        TimeoutException: If the time limit is exceeded.
    """
    if seconds is None:
        # No timeout requested: run block normally
        yield
        return

    def _handler(signum, frame):
        raise TimeoutException("Timed out!")

    # Register alarm; ensure it’s cleared afterward
    signal.signal(signal.SIGALRM, _handler)
    signal.alarm(int(seconds))
    try:
        yield
    finally:
        signal.alarm(0)

# ===================================
# CLI Parsing & Simple Validations
# ===================================
def parse_args():
    """
    CLI for generating systems and data.

    Returns:
        argparse.Namespace with fields used throughout generation:
          vars, derivs, numVars, numDerivs, numEquations, numReplacements,
          numSystems, numConstConseq, genReplacements, genSysData,
          genConsequence, genConsequenceData, conseqDataRange, sysDataRange,
          maxConseqTerms, timeout, seed.
    """
    parser = argparse.ArgumentParser(description="Generate dimensionally consistent systems for dataset.")

    parser.add_argument('--numVars', type=ast.literal_eval, default="[6,7,8,9]",
                        help="List or single number of variables, e.g., '[6,7]' or '[6]'")
    parser.add_argument('--numDerivs', type=ast.literal_eval, default="[2,3,4]",
                        help="List or single number of derivatives, e.g., '[2,3]' or '[2]'")
    parser.add_argument('--numEquations', type=ast.literal_eval, default="[4,5,6]",
                        help="List of number of equations to try")
    parser.add_argument('--vars', type=ast.literal_eval,
                        default="['Fc', 'Fg', 'W', 'd1', 'd2', 'm1', 'm2', 'p', 'T', 'E']",
                        help="List of variables to use, e.g. ['Fc','Fg','d1']")
    parser.add_argument('--derivs', type=ast.literal_eval,
                        default="['dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2']",
                        help="List of derivatives to use, e.g. ['dx1dt','d2x1dt2']")

    parser.add_argument('--numReplacements', type=int, default=5)
    parser.add_argument('--numSystems', type=int, default=3)

    parser.add_argument('--numConstConseq', type=ast.literal_eval, default="1",
    help="Maximum number of constants allowed in the consequence. E.g 0 or 1"
)


    parser.add_argument('--conseqDataRange', type=ast.literal_eval, default="[1, 10]",
                        help='Range over which to sample variables when generating consequence data, e.g., "[1, 10]"'
    )
    parser.add_argument("--sysDataRange", type=ast.literal_eval, default="[1, 10]",
                        help='Variable range for system data generation, e.g., "[1, 10]"'
    )

    parser.add_argument('--genReplacements', type=ast.literal_eval, default=True)
    parser.add_argument('--genSysData', type=ast.literal_eval, default=True)
    parser.add_argument('--genConsequence', type=ast.literal_eval, default=True)
    parser.add_argument('--genConsequenceData', type=ast.literal_eval, default=True)

    parser.add_argument("--maxConseqTerms", type=ast.literal_eval, default="8")

    parser.add_argument('--timeout', type=int, default=3600,
                        help="Max number of seconds to allow for generating a single system. None for no timeout.")

    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    return parser.parse_args()

def parse_and_validate_range(range_str, flag_name):
    """
    Validate a 2-element [low, high] range.

    Args:
        range_list: Parsed value from CLI (already list via ast.literal_eval).
        flag_name: Name to print in error messages.

    Returns:
        The validated [low, high] list.

    Raises:
        ValueError: On invalid shape or order.
    """
    try:
        r = range_str
        assert (
            isinstance(r, list)
            and len(r) == 2
            and all(isinstance(x, (int, float)) for x in r)
            and r[0] < r[1]
        )
        return r
    except Exception:
        raise ValueError(f"Invalid {flag_name}. Use format like \"[1, 5]\" with start < end.")

# ==========================================
# Compact save helpers for EquationSystem
# ==========================================
def _split_top_level_commas(s: str):
    """Split by commas ignoring those inside parentheses."""
    parts, depth, start = [], 0, 0
    for i, ch in enumerate(s):
        if ch == '(':
            depth += 1
        elif ch == ')':
            if depth > 0:
                depth -= 1
        elif ch == ',' and depth == 0:
            parts.append(s[start:i].strip())
            start = i + 1
    tail = s[start:].strip()
    if tail:
        parts.append(tail)
    return parts

def _parse_name_unit_pairs(section_text: str):
    """
    Parse pairs like 'a (1), b (m), c (kg)' or spaced 'a (1) b (m) c (kg)'.

    Returns:
        names, units (balanced parentheses allowed inside units).
    """
    text = section_text.replace("\n", " ").strip()
    names, units = [], []
    i, n = 0, len(text)

    while i < n:
        while i < n and text[i] in " ,":
            i += 1
        if i >= n:
            break
        m = re.match(r'[A-Za-z_]\w*', text[i:])
        if not m:
            break
        name = m.group(0)
        i += m.end()

        while i < n and text[i].isspace():
            i += 1

        if i >= n or text[i] != '(':
            break
        depth = 0
        start_unit = i + 1
        while i < n:
            ch = text[i]
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0:
                    end_unit = i
                    i += 1 
                    break
            i += 1
        unit = text[start_unit:end_unit].strip()
        names.append(name)
        units.append(unit)

        while i < n and text[i] in " ,":
            i += 1

    return names, units

def _parse_constants(section_text: str):
    """
    Parse constants lines like:
      'G = 6.67e-11 (1/kg*m**3*s**(-2)), c = 299792458.0 (1/s*m)'
    Returns:
      (['G','c'], ['1/kg*m**3*s**(-2)', '1/s*m'])
    """
    text = section_text.replace("\n", " ").strip()
    items = _split_top_level_commas(text)
    names, units = [], []
    for it in items:
        m = re.match(r'\s*([A-Za-z_]\w*)\s*=\s*(.*)$', it)
        if not m:
            continue
        nm, rest = m.group(1), m.group(2).strip()
        i = rest.find('(')
        if i == -1:
            continue
        # find matching ')'
        depth, j = 0, i
        while j < len(rest):
            if rest[j] == '(':
                depth += 1
            elif rest[j] == ')':
                depth -= 1
                if depth == 0:
                    break
            j += 1
        if j >= len(rest):
            continue
        unit = rest[i+1:j].strip()
        names.append(nm)
        units.append(unit)
    return names, units

def _format_units(u: str) -> str:
    """Normalize unit strings: '**' -> '^' and collapse whitespace."""
    u = (u or "").replace("**", "^")
    return re.sub(r'\s+', '', u)

def save_equation_system_to_file(eqn_system_or_str, filename="temp.txt"):
    """
    Write a compact, machine-parsable system file.

    Accepts either an EquationSystem object or its pretty string.

    Args:
        eqn_system_or_str: EquationSystem or str(equation system).
        filename: Destination path.
    """
    text_full = str(eqn_system_or_str) if not isinstance(eqn_system_or_str, str) else eqn_system_or_str
    text = (text_full or "").strip().replace("\r\n", "\n")

    # Locate top-level sections
    header_re = re.compile(r'^\s*(Variables|Derivatives|Constants|Equations)\s*:\s*$', flags=re.MULTILINE)
    matches = list(header_re.finditer(text))
    sections = {}
    for i, m in enumerate(matches):
        name = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        sections[name] = text[start:end].strip("\n")

    if not sections:
        with open(filename, "w", encoding="utf-8") as f:
            f.write(text_full)
        return

    # Variables
    vars_names, vars_units = [], []
    if sections.get("Variables"):
        vn, vu = _parse_name_unit_pairs(sections["Variables"])
        vars_names = vn
        vars_units = [_format_units(x) for x in vu]

    # Derivatives
    deriv_names, deriv_units = [], []
    if sections.get("Derivatives"):
        dn, du = _parse_name_unit_pairs(sections["Derivatives"])
        deriv_names = dn
        deriv_units = [_format_units(x) for x in du]

    # Constants
    const_names, const_units = [], []
    if sections.get("Constants"):
        cn, cu = _parse_constants(sections["Constants"])
        const_names = cn
        const_units = [_format_units(x) for x in cu]

    # Equations (optional "(U of M: ...)" per line)
    eq_exprs, eq_units = [], []
    if sections.get("Equations"):
        lines = [ln for ln in sections["Equations"].splitlines() if ln.strip()]
        eq_pat = re.compile(r'^(.*?)(?:\s*\(U of M:\s*(.*)\)\s*)?$')
        for ln in lines:
            m = eq_pat.match(ln)
            if not m:
                continue
            expr = (m.group(1) or "").replace("**", "^").strip()
            if not expr or expr.lower() == "equations:":
                continue
            unit = _format_units(m.group(2)) if m.group(2) else ""
            eq_exprs.append(expr)
            eq_units.append(unit)

    # Fall back to pretty string if nothing parsed
    if not (vars_names or deriv_names or const_names or eq_exprs):
        with open(filename, "w", encoding="utf-8") as f:
            f.write(text_full)
        return

    # Writing to file
    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"Variables: {vars_names}\n")
        f.write(f"Constants: {const_names}\n")
        f.write(f"Derivatives: {deriv_names}\n")
        f.write("Equations:\n")
        for expr in eq_exprs:
            f.write(expr + "\n")
        f.write(f"Units of Measure of Variables: {vars_units}\n")
        f.write(f"Units of Measure of Constants: {const_units}\n")
        f.write(f"Units of Measure of Derivatives: {deriv_units}\n")
        f.write("Units of Measure of Equations:\n")
        for unit in eq_units:
            f.write((unit or "") + "\n")

# ======================================================================
# Temporary files utilities (Managing scripts generated by subprocesses)
# ======================================================================
def create_replacement_files(replacement_info, num_replacements, original_file="temp.txt"):
    """
    Create `temp_replacement_#.txt` files by applying provided replacements
    to a compact system file.

    `replacement_info` format (sequence of lines):
        <equation_index (1-based)>
        <equation string> (U of M: <unit string with possible nested parens>)
        <equation_index>
        <equation string> (U of M: <unit>)
        ...

    Args:
        replacement_info: Text block produced by `replaceRandomDimensionallyConsistentEqn()`.
        num_replacements: How many replacements/files to emit (prefix of the provided list).
        original_file: Source compact file (typically the just-saved temp system).
    """
    with open(original_file, 'r') as f:
        original_content = f.read() 

    # Parse replacement_info into structured entries   
    replacements = []
    current_replacement = {}
    lines = replacement_info.split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.isdigit(): 
            if current_replacement:  
                replacements.append(current_replacement)
            current_replacement = {'index': int(line) - 1} 
        elif '(U of M:' in line:  
            parts = line.split('(U of M:')
            current_replacement['equation'] = parts[0].strip().replace('**', '^')
            
            # Find matching ')' for nested unit parens
            unit_part = parts[1]
            open_parens = 1
            end_pos = 0
            for i, c in enumerate(unit_part):
                if c == '(':
                    open_parens += 1
                elif c == ')':
                    open_parens -= 1
                    if open_parens == 0:
                        end_pos = i
                        break
            
            unit = unit_part[:end_pos].strip()
            current_replacement['unit'] = unit.replace('**', '^')
    
    if current_replacement:  # Add the last replacement
        replacements.append(current_replacement)
    
    # Ensure we have exactly num_replacements
    replacements = replacements[:num_replacements]
    
    # Re-scan original into sections while preserving order
    sections = {}
    current_section = None
    sections_order = []  
    
    for line in original_content.split('\n'):
        if line.startswith('Variables:') or line.startswith('Constants:') or line.startswith('Derivatives:'):
            current_section = line.split(':')[0]
            sections[current_section] = line + '\n'
            sections_order.append(current_section)
        elif line.startswith('Equations:'):
            current_section = 'Equations'
            sections[current_section] = [line + '\n']  
            sections_order.append(current_section)
        elif line.startswith('Units of Measure of'):
            current_section = line
            sections[current_section] = line + '\n'
            sections_order.append(current_section)
        elif current_section == 'Equations' and line.strip():
            sections['Equations'].append(line + '\n')
        elif current_section and current_section.startswith('Units of Measure of'):
            if line.strip() and not line.startswith('Units of Measure of'):
                sections[current_section] += line + '\n'
    
    # Apply replacements and write out single-file variants
    for k, replacement in enumerate(replacements[:num_replacements]):

        new_content = {k: v[:] if isinstance(v, list) else v for k, v in sections.items()}
        # Swap equation line
        eqn_index = replacement['index']
        if 'Equations' in new_content and eqn_index < len(new_content['Equations']) - 1:
            new_content['Equations'][eqn_index + 1] = replacement['equation'] + '\n'
        
        # Swap corresponding unit line
        units_line = None
        for section in new_content:
            if section.startswith('Units of Measure of Equations'):
                units_line = section
                break
        
        if units_line:
            unit_lines = new_content[units_line].split('\n')
            header = unit_lines[0]
            unit_values = unit_lines[1:-1]  # Skip header and last empty line
            
            if eqn_index < len(unit_values):
                unit_values[eqn_index] = replacement['unit']
            
            # Rebuild the units section
            new_content[units_line] = header + '\n'
            for unit in unit_values:
                if unit.strip():  # Only add non-empty lines
                    new_content[units_line] += unit + '\n'
        
        # Reconstruct in original section order
        output_lines = []
        for section in sections_order:
            if section in new_content:
                if isinstance(new_content[section], list):
                    output_lines.extend(new_content[section])
                else:
                    output_lines.append(new_content[section])
        
        replacement_file = f"temp_replacement_{k+1}.txt"
        with open(replacement_file, 'w') as f:
            f.write(''.join(output_lines))
        
        print(f"Created replacement file: {replacement_file}")

def cleanup_temp_files(system_num):
    """
    Remove temporary files emitted during a generation attempt for a given system index.
    """
    patterns = [
        f"temp_system_{system_num}.txt",
        f"temp_system_{system_num}.dat",
        f"temp_system_{system_num}_*.dat",
        f"temp_consequence_{system_num}.txt",
        f"temp_consequence_{system_num}_reordered.txt",
        f"temp_data_{system_num}.dat",
        f"temp_data_{system_num}_*.dat",
        "temp_replacement_*.txt"
    ]
    
    for pattern in patterns:
        for f in Path(".").glob(pattern):
            try:
                f.unlink()
            except:
                pass

def parse_eqnSystem_from_file(file_path):
    """
    Read a compact system file and reconstruct an EquationSystem.

    Args:
        file_path: Path to a compact 'system.txt'.

    Returns:
        EquationSystem instance ready for downstream use.
    """
    data = {
        'variables': [],
        'constants': [],
        'derivatives': [],
        'equations': [],
    }
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    def extract_list(section, prefix):
        list_str = section.split(prefix, 1)[1].split('\n', 1)[0].strip()
        try:
            return eval(list_str)
        except:
            return []
    
    # Split into sections and collect lists/equations
    sections = re.split(r'\n(?=\w+:)', content)
    for section in sections:
        if section.startswith('Variables:'):
            data['variables'] = extract_list(section, 'Variables:')
        elif section.startswith('Constants:'):
            data['constants'] = extract_list(section, 'Constants:')
        elif section.startswith('Derivatives:'):
            data['derivatives'] = extract_list(section, 'Derivatives:')
        elif section.startswith('Equations:'):
            eqn_content = section.split('Equations:', 1)[1]
            if 'Units of Measure of Variables:' in eqn_content:
                eqn_content = eqn_content.split('Units of Measure of Variables:', 1)[0]            
            equations = []
            for line in eqn_content.split('\n'):
                line = line.strip()
                if line and not line.startswith('Units of Measure of'):
                    equations.append(line)
            data['equations'] = equations

    # Build symbol objects
    vars = variables(','.join(data['variables'])) if data['variables'] else []
    derivs = derivatives(','.join(data['derivatives'])) if data['derivatives'] else []
    consts = constants(','.join(data['constants'])) if data['constants'] else []
    
    # Evaluate equation strings in the symbol context
    eqns = []
    for eqn_str in data['equations']:
        try:
            if 'Units of Measure of' in eqn_str:
                continue                
            eqn_str = eqn_str.replace('^', '**')
            
            # Create a dictionary of all symbols for evaluation
            symbols_dict = {}
            for var in vars:
                symbols_dict[var.name] = var
            for deriv in derivs:
                symbols_dict[deriv.name] = deriv
            for const in consts:
                symbols_dict[const.name] = const
            
            # Evaluate the equation string in the context of our symbols
            eqn_expr = eval(eqn_str, {'__builtins__': None}, symbols_dict)
            eqn = Equation(eqn_expr)
            eqns.append(eqn)

        except Exception as e:
            print(f"Error parsing equation '{eqn_str}': {str(e)}")
            continue

    eqnSys = EquationSystem(
        vars=vars,
        derivatives=derivs,
        constants=consts,
        equations=eqns,
        max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
    )
    return eqnSys

def get_next_system_num():
    """
    Return next available global System_N index under dataset/ (legacy helper).
    Prefer the per-config incrementer used inside run_generation().
    """
    system_num = 1
    while Path(f"dataset/System_{system_num}").exists():
        system_num += 1
    return system_num


# =========================================
# Rerun pipeline on an existing system dir
# =========================================
def run_generation_from_prior_file(system_directory: str, args: argparse.Namespace) -> None:
    """
    Rerun the downstream pipeline (replacements, data, consequence) for an existing system directory.

    Args:
        system_directory: Path like dataset/<config>/System_<n>.
        args: Parsed CLI args (flags determine which stages are run).
    """
    genReplacement = args.genReplacements
    genSysData = args.genSysData
    genConsequence = args.genConsequence
    genConsequenceData = args.genConsequenceData
    conseqDataRange = args.conseqDataRange
    sysDataRange = args.sysDataRange
    timeout_limit = args.timeout

    start_time = time.process_time()
    process = psutil.Process(os.getpid())
    Equation.SetLogging()

    stage_times: Dict[str, float] = {}
    mem_usage: Dict[str, float] = {}

    if not Path(system_directory).exists():
        print(f"Directory {system_directory} does not exist")
        return

    config_key = Path(system_directory).parent.name
    system_num = int(Path(system_directory).name.split('_')[-1])

    print(f"\n=== Processing existing system: {system_directory} ===")

    temp_system_file = f"temp_system_{system_num}.txt"
    temp_system_dat = f"temp_system_{system_num}.dat"
    temp_consequence = f"temp_consequence_{system_num}.txt"
    temp_data_dat = f"temp_data_{system_num}.dat"

    try:
        if not Path(f"{system_directory}/system.txt").exists():
            print("System file not found in directory")
            return

        eqnSystem = parse_eqnSystem_from_file(f"{system_directory}/system.txt")
        save_equation_system_to_file(str(eqnSystem), temp_system_file)

        with time_limit(timeout_limit):
            # Replacement Axioms
            t = time.process_time()
            if genReplacement:
                print(" Now generating replacement axioms. \n")
                replacement_info = ""
                for j in range(3):
                    index, eqn = eqnSystem.copy().replaceRandomDimensionallyConsistentEqn()
                    replacement_info += f"{index+1}\n{eqn}\n"
                    print(f"Replacement {j+1}: Equation {index+1} with {eqn}")
                create_replacement_files(replacement_info, 3, temp_system_file)
                print("Replacement Axioms Created. \n")
            stage_times['replacement_generation'] = time.process_time() - t
            mem_usage['replacement_generation'] = process.memory_info().rss / (1024 * 1024)

            # System data
            t = time.process_time()
            if genSysData:
                print("Generating system data. Searching for roots. \n")
                found_roots = run_noiseless_system_data_generation(temp_system_file, temp_system_dat, sysDataRange)
                if not found_roots:
                    print("Could not find roots to the system. Skipping. \n")
                    cleanup_temp_files(system_num)
                    return
                run_noisy_system_data_generation(temp_system_file, temp_system_dat)
                print("Noisy System Data Generated. \n")
            stage_times['system_root_search'] = time.process_time() - t
            mem_usage['system_root_search'] = process.memory_info().rss / (1024 * 1024)

            # Consequence
            t = time.process_time()
            if genConsequence:
                print("Generating consequence. \n")
                # Pass the same knobs as in run_generation for consistency
                found_consequence = run_consequence_generation(
                    temp_system_file,
                    temp_consequence,
                    args.numConstConseq,
                    args.maxConseqTerms
                )
                if not found_consequence:
                    print("Could not find a consequence. Skipping. \n")
                    cleanup_temp_files(system_num)
                    return
                print("Consequence generated. \n")
            stage_times['consequence_generation'] = time.process_time() - t
            mem_usage['consequence_generation'] = process.memory_info().rss / (1024 * 1024)

            # Consequence data
            t = time.process_time()
            if genConsequenceData:
                print("Generating noiseless data for consequence. ")
                if not Path(temp_consequence).exists():
                    print("No consequence found in dataset. Skipping")
                    cleanup_temp_files(system_num)
                    return
                found_roots = run_consequence_noiseless_data_generation(temp_consequence, temp_data_dat, conseqDataRange)
                if not found_roots:
                    print("Could not solve equation for roots. Skipping")
                    cleanup_temp_files(system_num)
                    return
                run_consequence_noisy_data_generation(temp_data_dat)
                print("Noisy consequence data generated. System Complete. \n")
            stage_times['consequence_data_generation'] = time.process_time() - t
            mem_usage['consequence_data_generation'] = process.memory_info().rss / (1024 * 1024)

            total_time = time.process_time() - start_time
            peak_memory = max(mem_usage.values())

            print("\nPerformance Summary:")
            print(f"Total time: {total_time:.2f} seconds")
            print(f"Peak memory: {peak_memory:.2f} MB")
            print("Stage times:")
            for stage, tv in stage_times.items():
                print(f"  {stage}: {tv:.2f}s")

            stats_dir = f"dataset_gen_statistics/{config_key}/System_{system_num}"
            Path(stats_dir).mkdir(parents=True, exist_ok=True)
            with open(f"{stats_dir}/performance.txt", "a") as f:
                f.write(f"\n--- RERUN STATS ---\n")
                f.write(f"Total time: {total_time:.2f} seconds\n")
                f.write(f"Peak memory: {peak_memory:.2f} MB\n")
                for stage, tv in stage_times.items():
                    f.write(f"{stage}_time: {tv:.2f}\n")
                for stage, mv in mem_usage.items():
                    f.write(f"{stage}_memory: {mv:.2f} MB\n")

    except TimeoutException:
        print(f"Timeout reached after {timeout_limit} seconds, skipping this rerun.")
        cleanup_temp_files(system_num)

    except Exception as e:
        print(f"Error during rerun: {str(e)}")
        cleanup_temp_files(system_num)

# =========================
# Main generation pipeline
# =========================        
def run_generation(args):
    """
    Top-level pipeline to enumerate configurations, generate systems, and run selected stages.

    Args:
        args: Parsed CLI args.
    """

    Equation.SetLogging()

    variable_options = args.vars
    derivative_options = args.derivs
    num_vars_options = args.numVars
    num_derivs_options = args.numDerivs
    num_eqns_options = args.numEquations
    numReplacements = args.numReplacements
    numSystems = args.numSystems
    numConstConseq = args.numConstConseq
    conseqDataRange = args.conseqDataRange
    sysDataRange = args.sysDataRange
    genReplacement = args.genReplacements
    genSysData = args.genSysData
    genConsequence = args.genConsequence
    genConsequenceData = args.genConsequenceData
    max_conseq_terms = args.maxConseqTerms
    seed = args.seed
    timeout_limit = args.timeout

    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
        print(f"Using random seed: {seed}")

    # output dir
    Path("dataset").mkdir(exist_ok=True)

    # Count existing systems per configuration (so that code may resume between runs)
    config_counts = {}

    # All parameter combinations
    param_combinations = list(itertools.product(num_vars_options, num_derivs_options, num_eqns_options))

    for num_vars, num_derivs, num_eqns in param_combinations:
        config_key = f"vars_{num_vars}_derivs_{num_derivs}_eqns_{num_eqns}"
        config_dir = f"dataset/{config_key}"
        Path(config_dir).mkdir(exist_ok=True)

        existing_systems = 0
        if Path(config_dir).exists():
            existing_systems = len([f for f in os.listdir(config_dir) if f.startswith("System_")])
        
        config_counts[config_key] = existing_systems
        
        print(f"\n=== Configuration: {num_vars} vars, {num_derivs} derivs, {num_eqns} eqns ===")
        print(f"Found {existing_systems} existing systems in {config_dir}")
        print(f"Adding {numSystems} new systems to this configuration")
        
        systems_added = 0
        attempts = 1

        # Keep trying until we add the requested number or hit a soft attempt cap
        while config_counts[config_key] < (numSystems + existing_systems) and attempts < 10:
            start_time = time.process_time()
            process = psutil.Process(os.getpid())
            print(f"\nAttempt {attempts} (Adding system: {systems_added+1}/{numSystems})")
            attempts+=1
        
            # Randomly sample symbols for this system
            var_names   = random.sample(variable_options, k=num_vars)
            deriv_names = random.sample(derivative_options, k=num_derivs)

            # If sin/cos are present without theta, swap in theta for a non-trig var (ensures identity validity)
            deps = {'sinTheta', 'cosTheta', 'eTheta'}
            if any(d in var_names for d in deps) and 'theta' not in var_names and 'theta' in variable_options:
                i = next((i for i, v in enumerate(var_names) if v not in deps), None)
                var_names[i if i is not None else random.randrange(len(var_names))] = 'theta'
            
            # Construct EquationSystem
            vars = variables(','.join(var_names))
            derivs = derivatives(','.join(deriv_names))
            G, c, pi = constants('G,c,pi')

            need_trig_identity = any(str(v) == "sinTheta" for v in vars) and any(str(v) == "cosTheta" for v in vars)
            base_num_eqns = num_eqns - 1 if need_trig_identity else num_eqns
            if base_num_eqns < 0:
                base_num_eqns = 0  
            
            print(f"\n Generating axiom system {config_counts[config_key]+1}/{numSystems+existing_systems} for config {config_key}")

            # Track time and memory at key stages
            stage_times = {}
            mem_usage = {}

            # Find the next available system number within this configuration
            def get_next_system_num_with_config():
                system_num = 1
                while Path(f"{config_dir}/System_{system_num}").exists():
                    system_num += 1
                return system_num
            
            system_num = get_next_system_num_with_config()
            
            try:
                with time_limit(timeout_limit):
                    # ----- 1) Generate random system (and optionally tack on trig identity) -----
                    t = time.process_time()

                    if base_num_eqns > 0:
                        eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(
                            vars=vars,
                            derivatives=derivs,
                            constants=[G, c],
                            measuredVars=vars,
                            numEqns=base_num_eqns,
                            max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
                        )
                    else:
                        # If num_eqns was 1 and we need the trig axiom, start with an empty system.
                        eqnSystem = EquationSystem(
                            vars=vars,
                            derivatives=derivs,
                            constants=[G, c],
                            measuredVars=vars,
                            equations=[],
                            max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
                        )

                    if need_trig_identity:
                        sin_v = EquationSystem._get_var_by_name(eqnSystem._vars, "sinTheta")
                        cos_v = EquationSystem._get_var_by_name(eqnSystem._vars, "cosTheta")
                        if sin_v is not None and cos_v is not None:
                            trig_identity = Equation(sin_v**2 + cos_v**2 - Integer(1))
                            eqnSystem.addEquation(trig_identity, keep_count=False)

                    eqn_system_str = str(eqnSystem)
                    stage_times['generation'] = time.process_time() - t
                    mem_usage['generation'] = process.memory_info().rss / (1024 * 1024)

                    print(f"Axiom System {system_num} generated.")
                    print(eqn_system_str)
                    print(f"Variables: {var_names}")
                    print(f"Derivatives: {deriv_names}")
                    print(f"Number of Equations: {num_eqns}")

                    temp_system_file = f"temp_system_{system_num}.txt"
                    save_equation_system_to_file(eqn_system_str, temp_system_file)

                    # ----- 2) Replacement axioms -----
                    t = time.process_time()
                    if genReplacement:
                        print("\n Now generating replacement axioms. \n")
                        replacement_info = ""
                        
                        if numReplacements > 0:
                            for j in range(numReplacements):
                                index, eqn = eqnSystem.copy().replaceRandomDimensionallyConsistentEqn()
                                replacement_info += f"{index+1}\n{eqn}\n"  # +1 to make it 1-based index
                                print(f"Replacement {j+1}: Equation {index+1} with {eqn}")

                        if numReplacements > 0:
                            create_replacement_files(replacement_info, numReplacements, temp_system_file)

                        print(" \n Replacement Axioms Created. \n")
                    stage_times['replacement_generation'] = time.process_time() - t
                    mem_usage['replacement_generation'] = process.memory_info().rss / (1024 * 1024)
                    
                    # ----- 3) System data -----
                    t = time.process_time()
                    if genSysData:
                        print("\n Generating system data. Searching for roots. \n")
                        temp_system_dat = f"temp_system_{system_num}.dat"
                        found_system_roots = run_noiseless_system_data_generation(temp_system_file, temp_system_dat, sysDataRange)
                        if not found_system_roots:
                            print("\n Could not find roots to the system. Skipping. \n")
                            cleanup_temp_files(system_num)
                            continue
                        print("\n System data generated. Now adding noise. \n")
                        
                        added_system_noise = run_noisy_system_data_generation(temp_system_file, temp_system_dat)
                        if not added_system_noise:
                            print("Failed to add noise")
                            cleanup_temp_files(system_num)
                            continue
                        print("\n Noisy System Data Generated. \n")
                    stage_times['system_root_search'] = time.process_time() - t
                    mem_usage['system_root_search'] = process.memory_info().rss / (1024 * 1024)
                    
                    # ----- 4) Algebraic consequence (M2 elimination) -----
                    t = time.process_time()
                    if genConsequence:
                        print("Generating consequence. \n")
                        temp_consequence = f"temp_consequence_{system_num}.txt"
                        found_consequence = run_consequence_generation(temp_system_file, temp_consequence, numConstConseq, max_conseq_terms)
                        if not found_consequence:
                            print("Could not find a consequence. Skipping. \n")
                            cleanup_temp_files(system_num)
                            continue
                        print("\n Consequence generated. \n")
                    stage_times['consequence_generation'] = time.process_time() - t
                    mem_usage['consequence_generation'] = process.memory_info().rss / (1024 * 1024)

                    # ----- 5) Data for consequence -----
                    t = time.process_time()
                    if genConsequenceData:
                        print("\n Generating noiseless data for consequence. ")
                        temp_data_dat = f"temp_data_{system_num}.dat"
                        temp_consequence = f"temp_consequence_{system_num}.txt"
                        if not os.path.exists(temp_consequence):
                            print("No consequence found in dataset. Skipping")
                            cleanup_temp_files(system_num)
                            continue
                        found_roots = run_consequence_noiseless_data_generation(temp_consequence, temp_data_dat, conseqDataRange)
                        if not found_roots:
                            print("Could not solve equation for roots. Skipping")
                            cleanup_temp_files(system_num)
                            continue
                        print("Noiseless consequence data generated. Adding noise \n")

                        run_consequence_noisy_data_generation(temp_data_dat)
                        print("Noisy consequence data generated. System Complete. \n")
                    stage_times['consequence_data_generation'] = time.process_time() - t
                    mem_usage['consequence_data_generation'] = process.memory_info().rss / (1024 * 1024)

                    # ----- 6) Write outputs and stats -----
                    total_time = time.process_time() - start_time
                    peak_memory = max(mem_usage.values())
                    print("\nPerformance Summary:")
                    print(f"Total time: {total_time:.2f} seconds")
                    print(f"Peak memory: {peak_memory:.2f} MB")
                    print("Stage times:")
                    for stage, t in stage_times.items():
                        print(f"  {stage}: {t:.2f}s")
                    
                    # All steps succeeded - create system directory within config directory
                    system_dir = f"{config_dir}/System_{system_num}"
                    Path(system_dir).mkdir(exist_ok=True)
                    
                    def safe_move(src, dst):
                        if Path(src).exists():
                            shutil.move(src, dst)

                    # Move system files
                    safe_move(temp_system_file, f"{system_dir}/system.txt")
                    safe_move(f"temp_system_{system_num}.dat", f"{system_dir}/system.dat")

                    # Move noisy system data files
                    for f in Path(".").glob(f"temp_system_{system_num}_*.dat"):
                        noise_level = f.stem.split('_')[-1]
                        safe_move(str(f), f"{system_dir}/system_{noise_level}.dat")

                    # Move consequence files
                    safe_move(f"temp_consequence_{system_num}.txt", f"{system_dir}/consequence.txt")
                    safe_move(f"temp_data_{system_num}.dat", f"{system_dir}/consequence.dat")

                    # Move noisy consequence data files
                    for f in Path(".").glob(f"temp_data_{system_num}_*.dat"):
                        noise_level = f.stem.split('_')[-1]
                        safe_move(str(f), f"{system_dir}/consequence_{noise_level}.dat")

                    # Move replacement files
                    for f in Path(".").glob("temp_replacement_*.txt"):
                        replacement_num = f.stem.split('_')[-1]
                        safe_move(str(f), f"{system_dir}/replacement_{replacement_num}.txt")

                    
                    stats_dir = f"dataset_gen_statistics/{config_key}/System_{system_num}"
                    Path(stats_dir).mkdir(parents=True, exist_ok=True)

                    with open(f"{stats_dir}/performance.txt", "w") as f:
                        f.write(f"Total time: {total_time:.2f} seconds\n")
                        f.write(f"Peak memory: {peak_memory:.2f} MB\n")
                        for stage, t in stage_times.items():
                            f.write(f"{stage}_time: {t:.2f}\n")
                        for stage, m in mem_usage.items():
                            f.write(f"{stage}_memory: {m:.2f} MB\n")
                    
                    config_counts[config_key] += 1
                    systems_added += 1
                    print(f"\n Moved files to: {system_dir}")
                    print(f"Successfully completed system {config_counts[config_key]}/{numSystems+existing_systems} for this configuration")
                    attempts = 1

            except TimeoutException:
                print(f"Timeout reached after {timeout_limit} seconds, skipping this attempt.")
                cleanup_temp_files(system_num)
                continue

            except Exception as e:
                print(f"Error generating system: {str(e)}")
                cleanup_temp_files(system_num)
                continue

if __name__ == "__main__":

    args = parse_args()
    args.conseqDataRange = parse_and_validate_range(args.conseqDataRange, "--conseqDataRange")
    args.sysDataRange = parse_and_validate_range(args.sysDataRange, "--sysDataRange")
    run_generation(args)

    # Batch rerun utility (disabled by default)
    """
    base_dir = "dataset"
    all_system_dirs: List[Path] = []
    for config_dir in Path(base_dir).iterdir():
        if config_dir.is_dir():
            for system_dir in config_dir.iterdir():
                if system_dir.is_dir() and system_dir.name.startswith("System_"):
                    all_system_dirs.append(system_dir)

    for current_directory in all_system_dirs:
        print(f"\nProcessing directory: {current_directory} \n")
        run_generation_from_prior_file(str(current_directory), args)
    """
                