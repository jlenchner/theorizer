import random
from theorizer import *
from consequence_generation import *
from consequence_data_generation import *
from system_data_generation import *
from equationSystem import *
from equation import *
from m2_functions import *
from unitOfMeasure import *
from derivative import *
from constant import *
import re
import os
from pathlib import Path
import shutil
import time
import psutil 
import argparse
import ast
import signal
from contextlib import contextmanager

class TimeoutException(Exception):
    pass

@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

def parse_args():
    parser = argparse.ArgumentParser(description="Generate dimensionally consistent systems for benchmarking.")

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

    parser.add_argument('--genReplacements', type=ast.literal_eval, default=True)
    parser.add_argument('--genSysData', type=ast.literal_eval, default=True)
    parser.add_argument('--genConsequence', type=ast.literal_eval, default=True)
    parser.add_argument('--genConsequenceData', type=ast.literal_eval, default=True)

    parser.add_argument('--conseqDataRange', type=ast.literal_eval, default="[1, 10]",
        help='Range over which to sample variables when generating consequence data, e.g., "[1, 10]"'
    )
    parser.add_argument("--sysDataRange", type=ast.literal_eval, default="[1, 10]",
        help='Variable range for system data generation, e.g., "[1, 10]"'
    )

    parser.add_argument('--timeout', type=int, default=3600,
    help="Max number of seconds to allow for generating a single system. None for no timeout.")

    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')

    return parser.parse_args()

def parse_and_validate_range(range_str, flag_name):
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

def save_equation_system_to_file(eqn_system_str, filename="temp.txt"):
    # Initialize dictionaries to store parsed data
    data = {
        'Variables': {'names': [], 'units': []},
        'Constants': {'names': [], 'units': []},
        'Derivatives': {'names': [], 'units': []},
        'Equations': {'expressions': [], 'units': []}
    }
    
    # Split the input string into sections
    sections = re.split(r'\n\n+', eqn_system_str.strip())
    
    # Helper function to parse units with proper handling of parentheses
    def parse_unit(unit_str):
        # Replace ** with ^ first
        unit_str = unit_str.replace('**', '^')
        return unit_str
    
    # Parse Variables section
    vars_section = sections[0]
    var_matches = re.finditer(r'(\w+) \((.*?)\)(?=\s*(?:,|$))', vars_section)
    for match in var_matches:
        data['Variables']['names'].append(match.group(1))
        data['Variables']['units'].append(parse_unit(match.group(2)))
    
    # Parse Derivatives section - improved regex to handle parentheses
    derivs_section = sections[1]
    deriv_matches = re.finditer(r'(\w+) \((.*?)\)(?=\s*(?:\w|$))', derivs_section)
    for match in deriv_matches:
        data['Derivatives']['names'].append(match.group(1))
        data['Derivatives']['units'].append(parse_unit(match.group(2)))
    
    # Parse Constants section - improved regex to handle parentheses
    consts_section = sections[2]
    const_matches = re.finditer(r'(\w+) = [\d.e+-]+ \((.*?)\)(?=\s*(?:,|$))', consts_section)
    for match in const_matches:
        data['Constants']['names'].append(match.group(1))
        data['Constants']['units'].append(parse_unit(match.group(2)))
    
    # Parse Equations section
    eqns_section = sections[3]
    # Split equations while preserving the U of M information
    eqn_blocks = re.split(r'\n(?=\S)', eqns_section.split('\n', 1)[1])
    for block in eqn_blocks:
        if '(U of M:' in block:
            # Split equation expression and unit
            parts = re.split(r'\(U of M:\s*(.*?)\)\s*$', block)
            expr = parts[0].strip()
            unit = parse_unit(parts[1])
            data['Equations']['expressions'].append(expr.replace('**', '^'))
            data['Equations']['units'].append(unit)
    
    # Write to file
    with open(filename, 'w') as f:
        # Write Variables
        f.write(f"Variables: {data['Variables']['names']}\n")
        
        # Write Constants
        f.write(f"Constants: {data['Constants']['names']}\n")
        
        # Write Derivatives
        f.write(f"Derivatives: {data['Derivatives']['names']}\n")
        
        # Write Equations
        f.write("Equations:\n")
        for expr in data['Equations']['expressions']:
            f.write(expr + '\n')
        
        # Write Units of Measure
        f.write(f"Units of Measure of Variables: {data['Variables']['units']}\n")
        f.write(f"Units of Measure of Constants: {data['Constants']['units']}\n")
        f.write(f"Units of Measure of Derivatives: {data['Derivatives']['units']}\n")
        
        # Write Equations Units
        f.write("Units of Measure of Equations:\n")
        for unit in data['Equations']['units']:
            f.write(unit + '\n')

def create_replacement_files(replacement_info, num_replacements, original_file="temp.txt"):
    # Read the original file
    with open(original_file, 'r') as f:
        original_content = f.read()
    
    # Parse the replacement info
    replacements = []
    current_replacement = {}
    lines = replacement_info.split('\n')
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.isdigit():  # This is an equation index
            if current_replacement:  # Save previous replacement if exists
                replacements.append(current_replacement)
            current_replacement = {'index': int(line) - 1}  # Convert to 0-based index
        elif '(U of M:' in line:  # Unit of measure
            # Extract the equation and unit - handle nested parentheses in units
            parts = line.split('(U of M:')
            current_replacement['equation'] = parts[0].strip().replace('**', '^')
            
            # Find the matching closing parenthesis for the unit
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
    
    # Parse the original file structure
    sections = {}
    current_section = None
    sections_order = []  # To maintain original section order
    
    for line in original_content.split('\n'):
        if line.startswith('Variables:') or line.startswith('Constants:') or line.startswith('Derivatives:'):
            current_section = line.split(':')[0]
            sections[current_section] = line + '\n'
            sections_order.append(current_section)
        elif line.startswith('Equations:'):
            current_section = 'Equations'
            sections[current_section] = [line + '\n']  # Start equations list
            sections_order.append(current_section)
        elif line.startswith('Units of Measure of'):
            current_section = line
            sections[current_section] = line + '\n'
            sections_order.append(current_section)
        elif current_section == 'Equations' and line.strip():
            sections['Equations'].append(line + '\n')
        elif current_section and current_section.startswith('Units of Measure of'):
            # Only add to units section if it's not the header line
            if line.strip() and not line.startswith('Units of Measure of'):
                sections[current_section] += line + '\n'
    
    # For each replacement, create a new file
    for k, replacement in enumerate(replacements[:num_replacements]):
        # Create a copy of the original content
        new_content = {k: v[:] if isinstance(v, list) else v for k, v in sections.items()}
        
        # Replace the specified equation
        eqn_index = replacement['index']
        if 'Equations' in new_content and eqn_index < len(new_content['Equations']) - 1:
            new_content['Equations'][eqn_index + 1] = replacement['equation'] + '\n'
        
        # Replace the corresponding unit of measure
        units_line = None
        for section in new_content:
            if section.startswith('Units of Measure of Equations'):
                units_line = section
                break
        
        if units_line:
            # Split units section into lines
            unit_lines = new_content[units_line].split('\n')
            # The first line is the header "Units of Measure of Equations:"
            header = unit_lines[0]
            unit_values = unit_lines[1:-1]  # Skip header and last empty line
            
            if eqn_index < len(unit_values):
                unit_values[eqn_index] = replacement['unit']
            
            # Rebuild the units section
            new_content[units_line] = header + '\n'
            for unit in unit_values:
                if unit.strip():  # Only add non-empty lines
                    new_content[units_line] += unit + '\n'
        
        # Reconstruct the file content in original order
        output_lines = []
        for section in sections_order:
            if section in new_content:
                if isinstance(new_content[section], list):
                    output_lines.extend(new_content[section])
                else:
                    output_lines.append(new_content[section])
        
        # Write to replacement file
        replacement_file = f"temp_replacement_{k+1}.txt"
        with open(replacement_file, 'w') as f:
            f.write(''.join(output_lines))
        
        print(f"Created replacement file: {replacement_file}")

def cleanup_temp_files(system_num):
    """Clean up all temporary files for a given system number"""
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
    data = {
        'variables': [],
        'constants': [],
        'derivatives': [],
        'equations': [],
    }
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    # Helper function to extract list data
    def extract_list(section, prefix):
        list_str = section.split(prefix, 1)[1].split('\n', 1)[0].strip()
        try:
            return eval(list_str)
        except:
            return []
    
    # Split into main sections
    sections = re.split(r'\n(?=\w+:)', content)
    
    # Process each section
    for section in sections:
        if section.startswith('Variables:'):
            data['variables'] = extract_list(section, 'Variables:')
        elif section.startswith('Constants:'):
            data['constants'] = extract_list(section, 'Constants:')
        elif section.startswith('Derivatives:'):
            data['derivatives'] = extract_list(section, 'Derivatives:')
        elif section.startswith('Equations:'):
            # Get everything after "Equations:" and before "Units of Measure"
            eqn_content = section.split('Equations:', 1)[1]
            if 'Units of Measure of Variables:' in eqn_content:
                eqn_content = eqn_content.split('Units of Measure of Variables:', 1)[0]
            
            # Take all non-empty lines that don't look like units sections
            equations = []
            for line in eqn_content.split('\n'):
                line = line.strip()
                if line and not line.startswith('Units of Measure of'):
                    equations.append(line)
            data['equations'] = equations

    # Create variables
    vars = variables(','.join(data['variables'])) if data['variables'] else []
    
    # Create derivatives
    derivs = derivatives(','.join(data['derivatives'])) if data['derivatives'] else []
    
    # Create constants
    consts = constants(','.join(data['constants'])) if data['constants'] else []
    
    # Create equations by parsing the equation strings
    eqns = []
    for eqn_str in data['equations']:
        try:
            # Skip any lines that look like units
            if 'Units of Measure of' in eqn_str:
                continue
                
            # Convert ^ to ** for Python exponentiation
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
    
    # Create and return the EquationSystem
    eqnSys = EquationSystem(
        vars=vars,
        derivatives=derivs,
        constants=consts,
        equations=eqns,
        max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
    )
    
    return eqnSys

def get_next_system_num():
    system_num = 1
    while Path(f"Benchmarking/System_{system_num}").exists():
        system_num += 1
    return system_num

def run_generation_from_prior_file(system_directory, args):

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

    stage_times = {}
    mem_usage = {}

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
                found_consequence = run_consequence_generation(temp_system_file, temp_consequence)
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
            for stage, t in stage_times.items():
                print(f"  {stage}: {t:.2f}s")

            with open(f"{system_directory}/performance.txt", "a") as f:
                f.write(f"\n--- RERUN STATS ---\n")
                f.write(f"Total time: {total_time:.2f} seconds\n")
                f.write(f"Peak memory: {peak_memory:.2f} MB\n")
                for stage, t in stage_times.items():
                    f.write(f"{stage}_time: {t:.2f}\n")
                for stage, m in mem_usage.items():
                    f.write(f"{stage}_memory: {m:.2f} MB\n")

    except TimeoutException:
        print(f"Timeout reached after {timeout_limit} seconds, skipping this rerun.")
        cleanup_temp_files(system_num)

    except Exception as e:
        print(f"Error during rerun: {str(e)}")
        cleanup_temp_files(system_num)

def run_generation(args):
    Equation.SetLogging()

    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Using random seed: {args.seed}")

    variable_options = args.vars
    derivative_options = args.derivs
    num_vars_options = args.numVars
    num_derivs_options = args.numDerivs
    num_eqns_options = args.numEquations
    numReplacements = args.numReplacements
    numSystems = args.numSystems
    genReplacement = args.genReplacements
    genSysData = args.genSysData
    genConsequence = args.genConsequence
    genConsequenceData = args.genConsequenceData
    conseqDataRange = args.conseqDataRange
    sysDataRange = args.sysDataRange
    timeout_limit = args.timeout

    # Create Benchmarking directory if it doesn't exist
    Path("Benchmarking").mkdir(exist_ok=True)

    # Track completed systems per configuration
    config_counts = {}

    # Generate all parameter combinations
    param_combinations = list(itertools.product(num_vars_options, num_derivs_options, num_eqns_options))

    for num_vars, num_derivs, num_eqns in param_combinations:
        config_key = f"vars_{num_vars}_derivs_{num_derivs}_eqns_{num_eqns}"
        config_dir = f"Benchmarking/{config_key}"
        Path(config_dir).mkdir(exist_ok=True)

        existing_systems = 0
        if Path(config_dir).exists():
            existing_systems = len([f for f in os.listdir(config_dir) if f.startswith("System_")])
        
        config_counts[config_key] = existing_systems
        
        print(f"\n=== Starting configuration: {num_vars} vars, {num_derivs} derivs, {num_eqns} eqns ===")
        attempts = 1

        while config_counts[config_key] < numSystems and attempts < 10:
            # Start timing
            start_time = time.process_time()
            process = psutil.Process(os.getpid())
            print("Attempt ", attempts, "\n")
            attempts+=1
            # Generate variables and derivatives dynamically based on current config
            var_names = variable_options[:num_vars]
            deriv_names = derivative_options[:num_derivs]
            
            # Create the symbols
            vars = variables(','.join(var_names))
            derivs = derivatives(','.join(deriv_names))
            G, c, pi = constants('G,c,pi')
            
            print(f"\nGenerating system {config_counts[config_key]+1}/{numSystems} for config {config_key}")

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
                    ##### Generating system #####
                    t = time.process_time()
                    eqnSystem = EquationSystem.GenerateRandomDimensionallyConsistent(
                        vars=vars, 
                        derivatives=derivs,
                        constants=[G, c],
                        measuredVars=vars, 
                        numEqns=num_eqns,
                        max_vars_derivatives_and_constants_per_eqn=Equation.NO_MAX
                    )
                    eqn_system_str = str(eqnSystem)
                    stage_times['generation'] = time.process_time() - t
                    mem_usage['generation'] = process.memory_info().rss / (1024 * 1024)

                    print(f"System {system_num} generated.")
                    print(f"Variables: {var_names}")
                    print(f"Derivatives: {deriv_names}")
                    print(f"Equations: {num_eqns}")

                    temp_system_file = f"temp_system_{system_num}.txt"
                    save_equation_system_to_file(eqn_system_str, temp_system_file)

                    ##### Generating replacement axioms #####
                    t = time.process_time()
                    if genReplacement:
                        print(" Now generating replacement axioms. \n")
                        replacement_info = ""
                        
                        if numReplacements > 0:
                            for j in range(numReplacements):
                                index, eqn = eqnSystem.copy().replaceRandomDimensionallyConsistentEqn()
                                replacement_info += f"{index+1}\n{eqn}\n"  # +1 to make it 1-based index
                                print(f"Replacement {j+1}: Equation {index+1} with {eqn}")

                        # Generate replacement files if needed
                        if numReplacements > 0:
                            create_replacement_files(replacement_info, numReplacements, temp_system_file)

                        print("Replacement Axioms Created. \n")
                    stage_times['replacement_generation'] = time.process_time() - t
                    mem_usage['replacement_generation'] = process.memory_info().rss / (1024 * 1024)
                    
                    ##### Generating system data #####
                    t = time.process_time()
                    if genSysData:
                        print("Generating system data. Searching for roots. \n")
                        temp_system_dat = f"temp_system_{system_num}.dat"
                        found_system_roots = run_noiseless_system_data_generation(temp_system_file, temp_system_dat, sysDataRange)
                        if not found_system_roots:
                            print("Could not find roots to the system. Skipping. \n")
                            cleanup_temp_files(system_num)
                            continue
                        print("System data generated. Now adding noise. \n")
                        
                        added_system_noise = run_noisy_system_data_generation(temp_system_file, temp_system_dat)
                        if not added_system_noise:
                            print("Failed to add noise")
                            cleanup_temp_files(system_num)
                            continue
                        print("Noisy System Data Generated. \n")
                    stage_times['system_root_search'] = time.process_time() - t
                    mem_usage['system_root_search'] = process.memory_info().rss / (1024 * 1024)
                    
                    ##### Generating consequence #####
                    t = time.process_time()
                    if genConsequence:
                        print("Generating consequence. \n")
                        temp_consequence = f"temp_consequence_{system_num}.txt"
                        found_consequence = run_consequence_generation(temp_system_file, temp_consequence)
                        if not found_consequence:
                            print("Could not find a consequence. Skipping. \n")
                            cleanup_temp_files(system_num)
                            continue
                        print("Consequence generated. \n")
                    stage_times['consequence_generation'] = time.process_time() - t
                    mem_usage['consequence_generation'] = process.memory_info().rss / (1024 * 1024)

                    ##### Generating data for consequence #####
                    t = time.process_time()
                    if genConsequenceData:
                        print("Generating noiseless data for consequence. ")
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

                    
                    stats_dir = f"Data_Gen_Statistics/{config_key}/System_{system_num}"
                    Path(stats_dir).mkdir(parents=True, exist_ok=True)

                    with open(f"{stats_dir}/performance.txt", "w") as f:
                        f.write(f"Total time: {total_time:.2f} seconds\n")
                        f.write(f"Peak memory: {peak_memory:.2f} MB\n")
                        for stage, t in stage_times.items():
                            f.write(f"{stage}_time: {t:.2f}\n")
                        for stage, m in mem_usage.items():
                            f.write(f"{stage}_memory: {m:.2f} MB\n")
                    
                    config_counts[config_key] += 1
                    print(f"Successfully completed system {config_counts[config_key]}/{numSystems} for this configuration")
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

    """
    base_dir = "Benchmarking"
    all_system_dirs = []

    # First collect all directories
    for config_dir in Path(base_dir).iterdir():
        if config_dir.is_dir():
            for system_dir in config_dir.iterdir():
                if system_dir.is_dir() and system_dir.name.startswith("System_"):
                    all_system_dirs.append(system_dir)

    # Then process them
    for current_directory in all_system_dirs:
        print(f"\nProcessing directory: {current_directory} \n")
        run_generation_from_prior_file(current_directory)
    """
                    
