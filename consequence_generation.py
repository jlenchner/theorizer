from m2_functions import *
import numpy as np
import os
import re
from collections import defaultdict
import random
from contextlib import redirect_stdout
import io

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

def run_consequence_generation(input_filepath, output_filepath, numConstConseq=[1]):
    parsed_data = parse_data(input_filepath)

    all_vars = parsed_data['variables'] + parsed_data['derivatives'] + parsed_data['constants']
    equations = parsed_data['equations']
    axiom_variables = [extract_variables(eq) for eq in equations]

    max_attempts = 10
    attempt = 0

    while attempt < max_attempts:
        print(f"Consequence generation attempt {attempt+1} out of {max_attempts}.")
        constants = parsed_data['constants']
        if constants:
            # Select a random subset of constants (could be more than 1)
            max_allowed = max(numConstConseq)
            selected_constants = random.sample(constants, min(max_allowed, len(constants)))
            other_constants = [c for c in constants if c not in selected_constants]
        else:
            selected_constants = []
            other_constants = []

        np.random.shuffle(all_vars)
        all_vars_reordered = [v for v in all_vars if v not in other_constants] + other_constants

        for j in range(1, len(all_vars_reordered)):
            temp_measured_vars = all_vars_reordered[:j]
            num_constants_used = sum(1 for c in selected_constants if c in temp_measured_vars)

            # If the number of selected constants used is not in the allowed list, skip
            if num_constants_used not in numConstConseq:
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
                
                if "matrix {}" not in content:
                    # Extract just the first polynomial from the Gröbner basis
                    poly_lines = []
                    write_poly = False
                    for line in content.split("\n"):
                        if write_poly:
                            if line.strip():  # Only process non-empty lines
                                # Take the first polynomial line and stop
                                polynomial = line.strip()
                                break
                        if "Polynomials of the Gröbner basis of the eliminated ideal:" in line:
                            print(line)  # Keep the print statement
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
