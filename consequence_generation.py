from m2_functions import *
import numpy as np
import os
import re
from collections import defaultdict
import random
from contextlib import redirect_stdout
import io

# ---------------------------------------------------------------------
# Utility: functions for loading data and polynomial processing
# ---------------------------------------------------------------------
def parse_data(file_path):
    """
    Parse a compact system file (your standardized `system.txt` format).

    Expected sections (examples):
      Variables: ['x','y',...]
      Constants: ['G','c',...]
      Derivatives: ['dx1dt','d2x2dt2',...]
      Equations:
      <eqn_1>
      <eqn_2>
      ...
      Units of Measure of Variables: ['m','s^(-1)',...]
      Units of Measure of Constants: [...]
      Units of Measure of Derivatives: [...]
      Units of Measure of Equations:
      <unit_eqn_1>
      <unit_eqn_2>
      ...

    Args:
        file_path: Path to the system file.

    Returns:
        dict with keys:
          - 'variables', 'constants', 'derivatives', 'equations'
          - 'var_units', 'const_units', 'deriv_units', 'eqn_units'
    """

    data = {
        'variables': [],
        'constants': [],
        'derivatives': [],
        'equations': [],
        'var_units': [],
        'const_units': [],
        'deriv_units': [],
        'eqn_units': []
    }
    
    with open(file_path, 'r') as file:
        content = file.read()
    
    def extract_list(section, prefix):
        """Safely eval a Python-list literal after a section prefix."""
        list_str = section.split(prefix, 1)[1].split('\n', 1)[0].strip()
        try:
            return eval(list_str)
        except:
            return []
    
    
    def extract_units(content, section_name):
        """
        Find a 'Units of Measure of <section_name>:' block and parse its list.
        Stops at the next 'Units of Measure of ...' or next <Word>:' section.
        """
        pattern = r'Units of Measure of ' + section_name + r':(.*?)(?=\n\s*Units of Measure of |\n\s*\w+:|$)'
        match = re.search(pattern, content, re.DOTALL)
        if match:
            units_str = match.group(1).strip()
            try:
                return eval(units_str)
            except:
                return []
        return []
    
    # Split at top-level section headers (e.g., "Variables:", "Equations:", etc.)
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
            # Take everything after "Equations:"; extract equations lines and (optionally) eqn-unit lines.
            eqn_content = section.split('Equations:', 1)[1]            
            if 'Units of Measure of Equations:' in eqn_content:
                eqn_part, units_part = eqn_content.split('Units of Measure of Equations:', 1)                
                equations = []
                for line in eqn_part.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('Units of Measure of'):
                        equations.append(line)
                data['equations'] = equations                
                eqn_units = [line.strip() for line in units_part.split('\n') if line.strip()]
                data['eqn_units'] = eqn_units
            else:
                # If no units section, just take all non-empty lines as equations
                equations = [line.strip() for line in eqn_content.split('\n') if line.strip()]
                data['equations'] = equations
    
    # Units for variables/constants/derivatives via separate precise matcher
    data['var_units'] = extract_units(content, 'Variables')
    data['const_units'] = extract_units(content, 'Constants')
    data['deriv_units'] = extract_units(content, 'Derivatives')
    
    return data

def extract_variables(equation):
    """
    Extract symbol-like tokens from a single equation, excluding numeric literals.

    Args:
        equation: A single line polynomial string.

    Returns:
        Unique variable-like names (order not guaranteed).
    """
    # Matches identifiers like d2x1dt2, G, m1, sinTheta; avoids numbers.
    variables = re.findall(r'\b(?!\d+\.?\d*\b)[a-zA-Z_][a-zA-Z0-9_]*\b', equation)
    return list(set(variables))

def sort_measured_components_by_degree(polynomial, parsed_data):
    """
    Order variables/constants/derivatives by their highest exponent in a polynomial.
    Only exact token matches are credited (no partial-name collisions).

    Args:
        polynomial: Single polynomial string (caret exponents allowed).
        parsed_data: Output from `parse_data()` with keys 'variables', 'constants', 'derivatives'.

    Returns:
        (vars_sorted, consts_sorted, derivs_sorted) — each list sorted by desc. degree.
    """
    degree_dict = defaultdict(int)
    
    # Get all possible components
    all_vars = parsed_data['variables']
    all_consts = parsed_data['constants']
    all_derivs = parsed_data['derivatives']
    all_components = all_vars + all_consts + all_derivs
    
    # Create a pattern that matches whole words only (exact matches)
    word_pattern = re.compile(r'\b(' + '|'.join(map(re.escape, all_components)) + r')\b')
    
    # Split polynomial into additive terms; normalize spaces out.
    terms = re.split(r'(?=[+-])', polynomial.replace(' ', ''))
    
    # For each term, find factors with optional ^exponent and record max exponent per component.
    for term in terms:
        if not term or term in ['+', '-']:
            continue
    
        factors = re.findall(r'([a-zA-Z_][a-zA-Z0-9_]*(?:\^\d+)?)', term)
        
        for factor in factors:
            if '^' in factor:
                base, exp = factor.split('^')
                exp = int(exp)
            else:
                base = factor
                exp = 1
            
            if base in all_components:
                if exp > degree_dict.get(base, 0):
                    degree_dict[base] = exp
    
    # Sort each category by degree (descending)
    def sort_by_degree(components):
        appearing = [c for c in components if c in degree_dict]
        return sorted(appearing, key=lambda x: degree_dict[x], reverse=True)
    
    return (sort_by_degree(all_vars),
            sort_by_degree(all_consts),
            sort_by_degree(all_derivs))

def count_poly_terms(poly_str: str) -> int:
    """
    Count additive terms in a single polynomial string (M2 GB line style).

    Args:
        poly_str: String like "x^2 - 3xy + y".

    Returns:
        Number of terms (3 for the example above).
    """
    s = poly_str.strip()
    if not s:
        return 0
    # Normalize so each additive term starts with optional sign
    # Example: "x^2 - 3xy + y" -> ["x^2", "- 3xy", "+ y"]
    parts = [p.strip() for p in re.split(r'(?=[+-])', s) if p.strip() not in ('+', '-')]
    return len(parts)

def reorder_measured_vars_for_save(sorted_vars):
    """
    Enforce the measured-variable ordering convention before saving:
      1) d1, d2 (if present)
      2) theta
      3) sinTheta, cosTheta, eTheta
      4) remaining vars in their existing order

    Args:
        sorted_vars: Variables already degree-sorted.

    Returns:
        New list with required prefix ordering applied.
    """
    prefix = []
    for v in ('d1', 'd2'):
        if v in sorted_vars:
            prefix.append(v)
    if 'theta' in sorted_vars:
        prefix.append('theta')
    for v in ('sinTheta', 'cosTheta', 'eTheta'):
        if v in sorted_vars:
            prefix.append(v)
    rest = [v for v in sorted_vars if v not in set(prefix)]
    return prefix + rest

# ---------------------------------------------------------------------
# Consequence generation: functions for generating consequences of the axiom system
# ---------------------------------------------------------------------
def run_consequence_generation(input_filepath, output_filepath, numConstConseq=1, max_terms=8):
    """
    Compute an algebraic consequence via Macaulay2 elimination and emit a compact
    consequence spec file.

    Strategy (repeated up to `max_attempts`):
      1) Randomly sample a budgeted subset of constants (<= numConstConseq).
      2) Shuffle all symbols; sweep prefixes of that order to propose measured sets.
         - Enforce: only selected constants can appear in measured set.
         - Avoid trivial projections: reject if measured set is subset/superset of any axiom's vars.
      3) Call `projection(...)` (M2) to eliminate unmeasured symbols.
      4) From the output, pick the first GB polynomial line (if any).
      5) Accept if the number of additive terms <= max_terms.
      6) Sort measured components by inferred degree; apply display-order convention; write file.

    Args:
        input_filepath: Path to the source system file (your `system.txt` format).
        output_filepath: Path to write the consequence file (compact, machine-parsable).
        numConstConseq: Max number of constants allowed to appear in the measured set.
        max_terms: Reject consequences with more than this many additive terms.

    Returns:
        True if a consequence was written; False otherwise.
    """
    parsed_data = parse_data(input_filepath)

    all_vars = parsed_data['variables'] + parsed_data['derivatives'] + parsed_data['constants']
    equations = parsed_data['equations']
    axiom_variables = [extract_variables(eq) for eq in equations]

    max_attempts = 10
    attempt = 0

    while attempt < max_attempts:
        print(f"Consequence generation attempt {attempt + 1} out of {max_attempts}.")
        constants = parsed_data['constants']

        # (1) Sample which constants are allowed to participate in the measured set
        if constants:
            selected_constants = random.sample(constants, min(numConstConseq, len(constants)))
            other_constants = [c for c in constants if c not in selected_constants]
        else:
            selected_constants = []
            other_constants = []

        # Shuffle symbol order and push unselected constants to the end
        np.random.shuffle(all_vars)
        all_vars_reordered = [v for v in all_vars if v not in other_constants] + other_constants

        # (2) Sweep prefixes as candidate measured sets
        for j in range(1, len(all_vars_reordered)):
            temp_measured_vars = all_vars_reordered[:j]
            num_constants_used = sum(1 for c in selected_constants if c in temp_measured_vars)

            # Enforce constant budget
            if num_constants_used > numConstConseq:
                continue
            # Ensure no unselected constants leak in
            if any(c in temp_measured_vars for c in other_constants):
                continue
            # Avoid trivial projections
            if any(set(temp_measured_vars).issubset(set(axiom_vars)) for axiom_vars in axiom_variables) \
               or any(set(axiom_vars).issubset(set(temp_measured_vars)) for axiom_vars in axiom_variables):
                continue
            
            # (3) Run M2 projection silently; eliminate unmeasured symbols
            with redirect_stdout(io.StringIO()):
                projection(all_vars_reordered, equations, temp_measured_vars, 
                           list(set(all_vars_reordered) - set(temp_measured_vars)), 
                           filename='temp_proj.txt')

            # (4) Parse the temporary M2 output
            with open('temp_proj.txt', 'r') as temp_file:
                content = temp_file.read()
                os.remove('temp_proj.txt')

                if "matrix {}" in content or "map(R^1, R^0, 0)" in content:
                    # empty elimination; try another slice
                    continue

                # Extract the first polynomial line after the label
                polynomial = ""
                write_poly = False
                for line in content.split("\n"):
                    if write_poly:
                        if line.strip():
                            polynomial = line.strip()
                            break
                    if "Polynomials of the Gröbner basis of the eliminated ideal:" in line:
                        write_poly = True

                if not polynomial:
                    # no polynomial recovered; try again
                    continue
                
                # (5) Sparsity check
                n_terms = count_poly_terms(polynomial)
                if n_terms > max_terms:
                    # skip this consequence; seek a sparser one
                    continue

                # (6) Degree-based reordering for display before saving the result
                sorted_vars, sorted_consts, sorted_derivs = sort_measured_components_by_degree(polynomial, parsed_data)
                sorted_vars = reorder_measured_vars_for_save(sorted_vars)

                with open(output_filepath, 'w') as out_file:
                    out_file.write(f"Variables: {parsed_data['variables']}\n")
                    out_file.write(f"Constants: {parsed_data['constants']}\n")
                    out_file.write(f"Derivatives: {parsed_data['derivatives']}\n")
                    out_file.write("Equations:\n")
                    for eq in equations:
                        out_file.write(eq + "\n")

                    # Units of measure
                    out_file.write(f"Units of Measure of Variables: {parsed_data['var_units']}\n")
                    out_file.write(f"Units of Measure of Constants: {parsed_data['const_units']}\n")
                    out_file.write(f"Units of Measure of Derivatives: {parsed_data['deriv_units']}\n")
                    out_file.write("Units of Measure of Equations:\n")
                    for unit in parsed_data['eqn_units']:
                        out_file.write(unit + "\n")

                    # Measured components
                    out_file.write(f"\nMeasured Variables: {sorted_vars}\n")
                    out_file.write(f"Observed Constants: {sorted_consts}\n")
                    out_file.write(f"Measured Derivatives: {sorted_derivs}\n")

                    # Target polynomial
                    out_file.write("\nTarget Polynomial:\n")
                    out_file.write(polynomial + "\n")

                return True # Success

        # If we get here, try another shuffle/constant selection
        if parsed_data['constants'] and len(parsed_data['constants']) > 1:
            attempt += 1
            continue
        # If there are no (or only one) constants, reshuffling may not help; still advance attempt.
        else:
            break

    return False #Generation failed
