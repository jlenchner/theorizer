from m2_functions import *
import numpy as np
import os
import re
from collections import defaultdict
import random
from contextlib import redirect_stdout
import io
import signal

# ADD THESE CLASSES
class TimeoutException(Exception):
    pass

from contextlib import contextmanager

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

# Mapping from derivative names to their base variables
# This is used to ensure that when a derivative is in measured vars,
# its base variable is also included (required for ODE-based data generation)
DERIV_TO_BASE = {
    'dx1dt': 'd1',
    'd2x1dt2': 'd1',
    'dx2dt': 'd2',
    'd2x2dt2': 'd2',
}

# All known derivative names
ALL_DERIVATIVE_NAMES = {'dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2'}


def validate_derivative_in_polynomial(polynomial, derivs_in_poly):
    """
    Validate that derivatives in the polynomial satisfy:
    - At least one term does NOT contain the derivative (so derivative=0 is not trivial)
    
    Args:
        polynomial: String representation of the polynomial
        derivs_in_poly: List of derivative names that appear in this polynomial
    
    Returns:
        (is_valid, reason) tuple
    """
    if not derivs_in_poly:
        return True, "No derivatives"
    
    # Split polynomial into terms
    poly_normalized = polynomial.replace(' ', '').replace('^', '**')
    terms = re.split(r'(?=[+-])', poly_normalized)
    terms = [t for t in terms if t and t not in ['+', '-']]
    
    for deriv in derivs_in_poly:
        # Check: At least one term does NOT contain this derivative
        # (so that derivative=0 is not a trivial solution)
        has_term_without_deriv = False
        deriv_pattern = rf'\b{re.escape(deriv)}\b'
        for term in terms:
            if not re.search(deriv_pattern, term):
                has_term_without_deriv = True
                break
        
        if not has_term_without_deriv:
            return False, f"All terms contain {deriv} (trivial solution: {deriv}=0)"
    
    return True, "Valid"

def validate_dependent_variable_in_polynomial(polynomial, measured_vars):
    """
    Validate that dependent variables (d1, d2) don't appear in ALL terms.
    If a dependent variable appears in every term, it can be factored out,
    leading to trivial solutions or unsolvable constraints.
    
    Args:
        polynomial: String representation of the polynomial
        measured_vars: List of measured variable names
    
    Returns:
        (is_valid, reason) tuple
    """
    # Check for d1 and d2
    dependent_vars = [v for v in measured_vars if v in ['d1', 'd2']]
    
    if not dependent_vars:
        return True, "No dependent variables"
    
    # Split polynomial into terms
    poly_normalized = polynomial.replace(' ', '').replace('^', '**')
    terms = re.split(r'(?=[+-])', poly_normalized)
    terms = [t for t in terms if t and t not in ['+', '-']]
    
    for dep_var in dependent_vars:
        # Check: At least one term does NOT contain this dependent variable
        has_term_without_dep_var = False
        dep_var_pattern = rf'\b{re.escape(dep_var)}\b'
        for term in terms:
            if not re.search(dep_var_pattern, term):
                has_term_without_dep_var = True
                break
        
        if not has_term_without_dep_var:
            return False, f"All terms contain {dep_var} (can be factored out, leads to trivial/degenerate solutions)"
    
    return True, "Valid"

def validate_derivative_power(polynomial):
    """
    Validate that derivatives don't appear with power > 2.
    Higher powers lead to complex solutions and numerical instability.
    
    Args:
        polynomial: String representation of the polynomial
    
    Returns:
        (is_valid, reason) tuple
    """
    # Check for derivatives with exponent > 2
    for deriv in ALL_DERIVATIVE_NAMES:
        # Look for patterns like dx1dt^3, dx1dt^4, etc.
        pattern = rf'\b{re.escape(deriv)}\^([3-9]|\d{{2,}})\b'
        match = re.search(pattern, polynomial)
        if match:
            power = match.group(1)
            return False, f"Derivative {deriv} has power {power} > 2 (leads to complex/unstable solutions)"
    
    return True, "Valid derivative powers"

def parse_data(file_path):
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
    
    # Helper function to extract list data
    def extract_list(section, prefix):
        list_str = section.split(prefix, 1)[1].split('\n', 1)[0].strip()
        try:
            return eval(list_str)
        except:
            return []
    
    # Helper function to extract units
    def extract_units(content, section_name):
        pattern = r'Units of Measure of ' + section_name + r':(.*?)(?=\n\s*Units of Measure of |\n\s*\w+:|$)'
        match = re.search(pattern, content, re.DOTALL)
        if match:
            units_str = match.group(1).strip()
            try:
                return eval(units_str)
            except:
                return []
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
            # Get everything after "Equations:" and before next section
            eqn_content = section.split('Equations:', 1)[1]
            
            # Split into equations and units sections
            if 'Units of Measure of Equations:' in eqn_content:
                eqn_part, units_part = eqn_content.split('Units of Measure of Equations:', 1)
                
                # Process equations - take only lines that don't start with "Units of Measure"
                equations = []
                for line in eqn_part.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('Units of Measure of'):
                        equations.append(line)
                data['equations'] = equations
                
                # Process equation units
                eqn_units = [line.strip() for line in units_part.split('\n') if line.strip()]
                data['eqn_units'] = eqn_units
            else:
                # If no units section, just take all non-empty lines as equations
                equations = [line.strip() for line in eqn_content.split('\n') if line.strip()]
                data['equations'] = equations
    
    # Extract units separately with more precise matching
    data['var_units'] = extract_units(content, 'Variables')
    data['const_units'] = extract_units(content, 'Constants')
    data['deriv_units'] = extract_units(content, 'Derivatives')
    
    return data

def extract_variables(equation):
    # Improved variable extraction that handles exponents and derivatives
    variables = re.findall(r'\b(?!\d+\.?\d*\b)[a-zA-Z_][a-zA-Z0-9_]*\b', equation)
    return list(set(variables))

def sort_measured_components_by_degree(polynomial, parsed_data):
    """Sort all components by their highest degree in the polynomial with exact matching"""
    degree_dict = defaultdict(int)
    
    # Get all possible components
    all_vars = parsed_data['variables']
    all_consts = parsed_data['constants']
    all_derivs = parsed_data['derivatives']
    all_components = all_vars + all_consts + all_derivs
    
    # Create a pattern that matches whole words only (exact matches)
    word_pattern = re.compile(r'\b(' + '|'.join(map(re.escape, all_components)) + r')\b')
    
    # Split polynomial into terms
    terms = re.split(r'(?=[+-])', polynomial.replace(' ', ''))
    
    for term in terms:
        if not term or term in ['+', '-']:
            continue
            
        # Find all components in this term with their exponents
        components_in_term = {}
        # Split into factors while preserving exponents
        factors = re.findall(r'([a-zA-Z_][a-zA-Z0-9_]*(?:\^\d+)?)', term)
        
        for factor in factors:
            # Split into base and exponent
            if '^' in factor:
                base, exp = factor.split('^')
                exp = int(exp)
            else:
                base = factor
                exp = 1
            
            # Only count if it's an exact match to a known component
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
    Robustly count additive terms of a Macaulay2 polynomial line.
    We assume GB polynomials are sums of monomials with +/-, no negative exponents.
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
    Apply required ordering for measured variables before saving:
      1) d1, d2 (if present)
      2) theta (if present)
      3) sinTheta, cosTheta, eTheta (if present, in that order)
      4) remaining vars in the order they already appear (sorted by degree desc)
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

def run_consequence_generation(input_filepath, output_filepath, numConstConseq=1, max_terms=8, numDerivConseq=1):
    parsed_data = parse_data(input_filepath)

    all_vars = parsed_data['variables'] + parsed_data['derivatives'] + parsed_data['constants']
    equations = parsed_data['equations']
    axiom_variables = [extract_variables(eq) for eq in equations]
    
    # Define which symbols are derivatives for limiting
    all_derivative_names = {'dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2'}
    system_derivatives = parsed_data['derivatives']  # derivatives in this system

    max_attempts = 20
    attempt = 0

    while attempt < max_attempts:
        print(f"Consequence generation attempt {attempt + 1} out of {max_attempts}.")
        
        # Handle constants - same as before
        constants = parsed_data['constants']
        if constants:
            selected_constants = random.sample(constants, min(numConstConseq, len(constants)))
            other_constants = [c for c in constants if c not in selected_constants]
        else:
            selected_constants = []
            other_constants = []
        
        # Handle derivatives - similar logic to constants
        if system_derivatives:
            selected_derivatives = random.sample(system_derivatives, min(numDerivConseq, len(system_derivatives)))
            other_derivatives = [d for d in system_derivatives if d not in selected_derivatives]
        else:
            selected_derivatives = []
            other_derivatives = []

        np.random.shuffle(all_vars)
        # Push both unselected constants AND unselected derivatives to the end
        all_vars_reordered = [v for v in all_vars if v not in other_constants and v not in other_derivatives] + other_derivatives + other_constants
        
        # Ensure base variables (d1, d2) come before their derivatives in the ordering
        # This increases the likelihood of valid measured variable combinations
        for deriv, base in DERIV_TO_BASE.items():
            if deriv in all_vars_reordered and base in all_vars_reordered:
                deriv_idx = all_vars_reordered.index(deriv)
                base_idx = all_vars_reordered.index(base)
                if base_idx > deriv_idx:
                    # Move base variable to just before the derivative
                    all_vars_reordered.remove(base)
                    all_vars_reordered.insert(deriv_idx, base)

        for j in range(1, len(all_vars_reordered)):
            temp_measured_vars = all_vars_reordered[:j]
            num_constants_used = sum(1 for c in selected_constants if c in temp_measured_vars)
            num_derivatives_used = sum(1 for d in selected_derivatives if d in temp_measured_vars)

            # Enforce constant budget
            if num_constants_used > numConstConseq:
                continue
            # Ensure no unselected constants leak in
            if any(c in temp_measured_vars for c in other_constants):
                continue
            # Enforce derivative budget
            if num_derivatives_used > numDerivConseq:
                continue
            # Ensure no unselected derivatives leak in
            if any(d in temp_measured_vars for d in other_derivatives):
                continue
            # Avoid trivial projections
            if any(set(temp_measured_vars).issubset(set(axiom_vars)) for axiom_vars in axiom_variables) \
               or any(set(axiom_vars).issubset(set(temp_measured_vars)) for axiom_vars in axiom_variables):
                continue
            
            # Ensure that if a derivative is in measured vars, its base variable is also included
            # This is required for ODE-based data generation
            derivs_in_measured = [v for v in temp_measured_vars if v in DERIV_TO_BASE]
            missing_base = False
            for deriv in derivs_in_measured:
                base_var = DERIV_TO_BASE[deriv]
                if base_var not in temp_measured_vars:
                    missing_base = True
                    break
            if missing_base:
                continue

            try:
                with time_limit(90):  # 90 second timeout per projection attempt
                    with redirect_stdout(io.StringIO()):
                        projection(all_vars_reordered, equations, temp_measured_vars, 
                                list(set(all_vars_reordered) - set(temp_measured_vars)), 
                                filename='temp_proj.txt')
                print("Projection Computed. Analyzing.")
            except TimeoutException:
                print(f"Projection timed out after 60 seconds, trying next variable combination.")
                # Clean up temp file if it exists
                if os.path.exists('temp_proj.txt'):
                    os.remove('temp_proj.txt')
                continue  # Skip to next iteration

            # Inspect projection result
            with open('temp_proj.txt', 'r') as temp_file:
                content = temp_file.read()
                os.remove('temp_proj.txt')

                if "matrix {}" in content or "map(R^1, R^0, 0)" in content:
                    # empty elimination; try another slice
                    print("No polynomial in projection. Skipping.")
                    continue

                # Extract the FIRST polynomial line after the label
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
                    print("Polynomial was not recovered.")
                    continue

                # --- NEW: only accept if #terms <= max_terms ---
                n_terms = count_poly_terms(polynomial)
                if n_terms > max_terms:
                    # skip this consequence; seek a sparser one
                    print("Number of terms too large. Skipping.")
                    continue
                
                # --- Validate derivative constraints ---
                # Find which derivatives appear in the polynomial
                derivs_in_poly = [d for d in ALL_DERIVATIVE_NAMES if re.search(rf'\b{re.escape(d)}\b', polynomial)]
                is_valid, reason = validate_derivative_in_polynomial(polynomial, derivs_in_poly)
                if not is_valid:
                    print(f"Derivative validation failed: {reason}. Skipping.")
                    continue

                # --- NEW: Validate derivative powers ---
                is_valid, reason = validate_derivative_power(polynomial)
                if not is_valid:
                    print(f"Derivative power validation failed: {reason}. Skipping.")
                    continue

                # Degree-based sorts
                sorted_vars, sorted_consts, sorted_derivs = sort_measured_components_by_degree(polynomial, parsed_data)

                temp_measured_vars_for_check = sorted_vars  # These are the measured vars for this polynomial
                is_valid, reason = validate_dependent_variable_in_polynomial(polynomial, temp_measured_vars_for_check)
                if not is_valid:
                    print(f"Dependent variable validation failed: {reason}. Skipping.")
                    continue


                # --- NEW: reorder measured variables per your rule ---
                sorted_vars = reorder_measured_vars_for_save(sorted_vars)

                # Save in the desired format
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

                print("Consequence found: ", polynomial)
                return True

        # Try another shuffle/constant selection/derivative selection
        if (parsed_data['constants'] and len(parsed_data['constants']) > 1) or \
           (parsed_data['derivatives'] and len(parsed_data['derivatives']) > 1):
            attempt += 1
            continue
        else:
            break

    return False

"""
def run_consequence_generation(input_filepath, output_filepath, numConstConseq=1):
    parsed_data = parse_data(input_filepath)

    all_vars = parsed_data['variables'] + parsed_data['derivatives'] + parsed_data['constants']
    equations = parsed_data['equations']
    axiom_variables = [extract_variables(eq) for eq in equations]

    max_attempts = 10
    attempt = 0

    while attempt < max_attempts:
        print(f"Consequence generation attempt {attempt + 1} out of {max_attempts}.")
        constants = parsed_data['constants']
        if constants:
            # Select a random subset of constants up to the allowed maximum
            selected_constants = random.sample(constants, min(numConstConseq, len(constants)))
            other_constants = [c for c in constants if c not in selected_constants]
        else:
            selected_constants = []
            other_constants = []

        np.random.shuffle(all_vars)
        all_vars_reordered = [v for v in all_vars if v not in other_constants] + other_constants

        for j in range(1, len(all_vars_reordered)):
            temp_measured_vars = all_vars_reordered[:j]
            num_constants_used = sum(1 for c in selected_constants if c in temp_measured_vars)

            # Enforce that the number of constants used is within the allowed maximum
            if num_constants_used > numConstConseq:
                continue

            # Also ensure that none of the unselected constants sneak in
            if any(c in temp_measured_vars for c in other_constants):
                continue
            
            # Check validity
            if any(set(temp_measured_vars).issubset(set(axiom_vars)) for axiom_vars in axiom_variables) or any(set(axiom_vars).issubset(set(temp_measured_vars)) for axiom_vars in axiom_variables):
                continue
            
            with redirect_stdout(io.StringIO()):
                projection(all_vars_reordered, equations, temp_measured_vars, 
                        list(set(all_vars_reordered) - set(temp_measured_vars)), 
                        filename='temp_proj.txt')
            
            # Check if projection was successful
            with open('temp_proj.txt', 'r') as temp_file:
                content = temp_file.read()
                os.remove('temp_proj.txt')
                
                if "matrix {}" not in content and "map(R^1, R^0, 0)" not in content:
                    # Extract just the first polynomial from the GrÃƒÂ¶bner basis
                    poly_lines = []
                    write_poly = False
                    for line in content.split("\n"):
                        if write_poly:
                            if line.strip():  # Only process non-empty lines
                                # Take the first polynomial line and stop
                                polynomial = line.strip()
                                break
                        if "Polynomials of the GrÃƒÂ¶bner basis of the eliminated ideal:" in line:
                            #print(line)  # Keep the print statement
                            write_poly = True
                    else:
                        polynomial = ""  # If no polynomial found
                    
                    # Sort measured components by degree in polynomial
                    sorted_vars, sorted_consts, sorted_derivs = sort_measured_components_by_degree(polynomial, parsed_data)
                    
                    # Skip if no observed constants in polynomial
                    #if not sorted_consts and exists_constant:
                    #    continue
                    
                    # Save in the desired format with target polynomial
                    with open(output_filepath, 'w') as out_file:
                        # Write original system information
                        out_file.write(f"Variables: {parsed_data['variables']}\n")
                        out_file.write(f"Constants: {parsed_data['constants']}\n")
                        out_file.write(f"Derivatives: {parsed_data['derivatives']}\n")
                        out_file.write("Equations:\n")
                        for eq in equations:
                            out_file.write(eq + "\n")
                        
                        # Write units of measure
                        out_file.write(f"Units of Measure of Variables: {parsed_data['var_units']}\n")
                        out_file.write(f"Units of Measure of Constants: {parsed_data['const_units']}\n")
                        out_file.write(f"Units of Measure of Derivatives: {parsed_data['deriv_units']}\n")
                        out_file.write("Units of Measure of Equations:\n")
                        for unit in parsed_data['eqn_units']:
                            out_file.write(unit + "\n")
                        
                        # Write measured components sorted by degree
                        out_file.write(f"\nMeasured Variables: {sorted_vars}\n")
                        out_file.write(f"Observed Constants: {sorted_consts}\n")
                        out_file.write(f"Measured Derivatives: {sorted_derivs}\n")
                        
                        # Write target polynomial
                        out_file.write("\nTarget Polynomial:\n")
                        out_file.write(polynomial + "\n")
                    
                    return True
        
        # Try with a different constant if available
        if parsed_data['constants'] and len(parsed_data['constants']) > 1:
            attempt += 1
            continue
        else:
            break
    
    return False  # No valid consequence found
"""
