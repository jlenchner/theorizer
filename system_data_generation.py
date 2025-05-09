import numpy as np
import sympy as sp
from sympy import symbols, sympify, Poly
import os
from collections import defaultdict
from m2_functions import *
import re
import ast

def parse_system_data(file_path):
    """
    Parses a system description file with variables, constants, derivatives, equations, and units of measure.
    Returns a dictionary with all extracted information.
    """
    with open(file_path, 'r') as file:
        content = file.read()

    system_data = {}
    
    # System number (if present)
    system_number = re.search(r'System number (\d+)', content)
    system_data['system_number'] = int(system_number.group(1)) if system_number else None
    
    # Variables
    variables = re.search(r'Variables:\s*(\[.*?\])', content)
    system_data['variables'] = ast.literal_eval(variables.group(1)) if variables else []
    
    # Constants
    constants = re.search(r'Constants:\s*(\[.*?\])', content)
    system_data['constants'] = ast.literal_eval(constants.group(1)) if constants else []
    
    # Derivatives
    derivatives = re.search(r'Derivatives:\s*(\[.*?\])', content)
    system_data['derivatives'] = ast.literal_eval(derivatives.group(1)) if derivatives else []
    
    # Equations
    equations = []
    equations_match = re.search(r'Equations:\n(.*?)(?=\nUnits of Measure of|$)', content, re.DOTALL)
    if equations_match:
        raw_eqs = equations_match.group(1).strip().split('\n')
        equations = [
            eq.strip().replace('^', '**').replace(' ', '')
            for eq in raw_eqs
            if eq.strip() and not eq.startswith('Units')
        ]
    system_data['equations'] = equations
    
    # Units of Measure
    def extract_units(section_name):
        match = re.search(rf'Units of Measure of {section_name}:\s*(\[.*?\]|.*?(?=\n\w+:|$))', content, re.DOTALL)
        if match:
            units_str = match.group(1).strip()
            if units_str.startswith('['):
                return ast.literal_eval(units_str)
            else:
                # Handle multi-line units for equations
                return [line.strip() for line in units_str.split('\n') if line.strip()]
        return []
    
    system_data['var_units'] = extract_units('Variables')
    system_data['const_units'] = extract_units('Constants')
    system_data['deriv_units'] = extract_units('Derivatives')
    system_data['eqn_units'] = extract_units('Equations')
    
    """
    # Measured components (if present)
    measured_vars = re.search(r'Measured Variables:\s*(\[.*?\])', content)
    system_data['measured_variables'] = ast.literal_eval(measured_vars.group(1)) if measured_vars else []
    
    observed_consts = re.search(r'Observed Constants:\s*(\[.*?\])', content)
    system_data['observed_constants'] = ast.literal_eval(observed_consts.group(1)) if observed_consts else []
    
    measured_derivs = re.search(r'Measured Derivatives:\s*(\[.*?\])', content)
    system_data['measured_derivatives'] = ast.literal_eval(measured_derivs.group(1)) if measured_derivs else []
    
    # Target Polynomial (if present)
    target_poly_match = re.search(r'Target Polynomial:\n(.*?)(?=\n\w+:|$)', content, re.DOTALL)
    system_data['target_polynomial'] = target_poly_match.group(1).strip() if target_poly_match else None
    """

    return system_data

def save_dataset(dataset, filename, delimiter='\t'):
    """
    Save dataset to .dat file with variable headers.
    
    Parameters:
        dataset (dict): Dictionary of {variable: numpy_array}
        filename (str): Output file path
        delimiter (str): Column separator (default: tab)
    """
    # Validate dataset
    if not dataset:
        raise ValueError("Dataset is empty")
    
    # Check consistent array lengths
    lengths = [len(arr) for arr in dataset.values()]
    if len(set(lengths)) > 1:
        raise ValueError("Inconsistent array lengths in dataset")
    
    # Prepare header and data
    headers = list(dataset.keys())
    data = np.column_stack([dataset[var] for var in headers])
    
    # Write to file
    with open(filename, 'w') as f:
        # Write header
        f.write(delimiter.join(headers) + '\n')
        
        # Write data
        np.savetxt(f, data, delimiter=delimiter, fmt='%.8e')

def get_variables_from_equations(equations):
    # Apply extract_variables to each equation in the list
    return [extract_variables(eq) for eq in equations]

def extract_variables(equation):
    variables = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', equation)
    return list(set(variables))

def get_base_variable(derivative):
    if derivative.startswith('dx') and 'dt' in derivative:
        parts = derivative.split('x')
        num_part = parts[1].split('dt')[0]
        num = num_part[0] if num_part else '1'
        return f'd{num}'
    elif derivative.startswith('d2x') and 'dt2' in derivative:
        parts = derivative.split('x')
        num_part = parts[1].split('dt')[0]
        num = num_part[0] if num_part else '1'
        return f'd{num}'
    else:
        return None

def collect_variables_in_sequence_and_reverse(equations):
    """
    Collects variables from each equation in sequence and reverses the list.
    Variables from equation 1 will be at the rightmost side,
    variables from equation 2 will be to the left of those from equation 1, and so on.
    
    Parameters:
        equations (list): List of equations as strings.
    
    Returns:
        list: Reversed list of variables collected in sequence.
    """
    all_vars = []
    
    # Get variables from each equation in sequence
    equation_variables = get_variables_from_equations(equations)
    
    # Add variables in sequence, avoiding duplicates
    for vars_list in equation_variables:
        for var in vars_list:
            if var not in all_vars:
                all_vars.append(var)
    
    # Reverse the list so variables from equation 1 are at the rightmost side
    all_vars.reverse()
    
    return all_vars

def evaluate_polynomials(dataset, equations, variable_order, tolerance=1e-5):
    """
    Evaluate how well the dataset satisfies the given polynomial equations.
    
    Parameters:
        dataset (dict): Dictionary of {variable: numpy_array}
        equations (list or str): Equation(s) as a list of strings or a single equation string.
                                 Each equation should be written as an expression equal to 0.
        variable_order (list): List of variable names in the expected order.
        tolerance (float): Maximum allowed absolute error.
    
    Returns:
        dict: Evaluation results with statistics including total points, number of valid points,
              pass rate (%), maximum absolute error, and indices of valid points.
    """
    # Ensure equations is a list, even if a single string is provided.
    if isinstance(equations, str):
        equations = [equations]
    
    # Create sympy symbols in the provided order.
    sym_vars = sp.symbols(' '.join(variable_order))
    
    # Build lambdified functions for each equation.
    eval_fns = []
    for eq in equations:
        try:
            # Replace caret with ** for exponentiation
            expr = sp.sympify(eq.replace('^', '**'))
            fn = sp.lambdify(sym_vars, expr, modules="numpy")
            eval_fns.append(fn)
        except Exception as e:
            raise ValueError(f"Could not parse equation '{eq}': {str(e)}")
    
    # Construct the data matrix using the provided variable order.
    try:
        data_matrix = np.column_stack([dataset[var] for var in variable_order])
    except KeyError as e:
        raise KeyError(f"Variable {e} is missing from the dataset. Check your variable_order list.")
    
    total_points = data_matrix.shape[0]
    valid_mask = np.ones(total_points, dtype=bool)
    max_error = 0.0
    
    # Evaluate each equation over all data points.
    for fn in eval_fns:
        try:
            # Pass each column as separate arguments.
            results = fn(*data_matrix.T)
            results = np.array(results)
            # If the results are complex, compute the absolute error as the sum of absolute real and imaginary parts.
            if np.iscomplexobj(results):
                error = np.abs(np.real(results)) + np.abs(np.imag(results))
            else:
                error = np.abs(results)
            max_error = max(max_error, np.nanmax(error))
            valid_mask &= (error <= tolerance)
        except Exception as e:
            print(f"Error during evaluation of an equation: {e}")
            valid_mask &= False

    valid_points = np.sum(valid_mask)
    pass_rate = valid_points / total_points * 100

    return {
        'total_points': total_points,
        'valid_points': valid_points,
        'pass_rate': pass_rate,
        'max_absolute_error': max_error,
        'valid_indices': np.where(valid_mask)[0].tolist()
    }


def generate_dataset_from_grobner_basis(grobner_basis, variables, observed_constants, measured_derivatives,
                                        constant_data=True, derivative_data=True, region=[1, 5],
                                        num_points=10000, seed=42, tol=1e-6):
    np.random.seed(seed)
    data = {}
    sym_vars = symbols(' '.join(variables))
    sorted_basis = sorted(grobner_basis, key=lambda eq: len(extract_variables(eq)))
    num_valid_points = num_points

    # Generate data for constants
    for const in observed_constants:
        if constant_data:
            value = 1
            data[const] = np.full(num_valid_points, value)
        else:
            data[const] = np.random.uniform(region[0], region[1], num_valid_points)

    # Generate data for derivatives
    for deriv in measured_derivatives:
        base_var = get_base_variable(deriv)
        if base_var in variables:
            if base_var in data:
                base_data = data[base_var]
                deriv_data = np.diff(base_data, append=base_data[0])
                data[deriv] = deriv_data
            else:
                data[deriv] = np.random.uniform(region[0], region[1], num_valid_points)
        else:
            data[deriv] = np.random.uniform(region[0], region[1], num_valid_points)

    # Process Gröbner basis equations
    for eq in sorted_basis:
        eq_expr = sympify(eq)
        eq_variables = extract_variables(eq)
        unknown_vars = [var for var in eq_variables if var not in data]

        if not unknown_vars:
            continue

        # Generate data for all but one unknown variable
        for var in unknown_vars[:-1]:
            data[var] = np.random.uniform(region[0], region[1], num_valid_points)

        last_var = unknown_vars[-1]
        last_var_sym = sym_vars[variables.index(last_var)]

        # Prepare polynomial coefficients
        try:
            poly = Poly(eq_expr, last_var_sym)
            coefficients = poly.all_coeffs()
        except (sp.PolynomialError, ValueError):
            print(f"Equation {eq} is not a valid polynomial in {last_var}. Skipping.")
            continue

        last_var_data = np.full(num_valid_points, np.nan)
        valid_indices = []

        for i in range(num_valid_points):
            subs_dict = {}
            for var in eq_variables:
                if var != last_var and var in data:
                    subs_dict[sym_vars[variables.index(var)]] = data[var][i]
            try:
                coeffs = [c.subs(subs_dict) for c in coefficients]
                coeffs = [complex(c) for c in coeffs]
                roots = np.roots(coeffs)
                real_roots = roots[np.isreal(roots)].real
                if real_roots.size > 0:
                    last_var_data[i] = real_roots[0]
                    valid_indices.append(i)
            except Exception as e:
                pass

        if valid_indices:
            valid_indices = np.array(valid_indices)
            for var in data:
                data[var] = data[var][valid_indices]
            data[last_var] = last_var_data[valid_indices]
            num_valid_points = len(valid_indices)
        else:
            print(f"No valid solutions for {last_var} in equation {eq}")
            return None

    # Generate remaining variables
    for var in variables:
        if var not in data:
            data[var] = np.random.uniform(region[0], region[1], num_valid_points)

    return data

def generate_data_for_system(equations, observed_constants, measured_derivatives, num_points=1000, seed=42,
                             constant_data=True, derivative_data=True, region=[1, 5]):
    variables = collect_variables_in_sequence_and_reverse(equations)
    projection(variables, equations, variables, [], filename='temp_grobner.txt')

    with open('temp_grobner.txt', 'r') as file:
        content = file.read()
        # Match everything after "Polynomials..." line until next section or end
        grobner_match = re.search(
            r'Polynomials of the Gröbner basis of the eliminated ideal:\s*\n(.*?)(?=\n\S+:|$)',
            content, 
            re.DOTALL
        )
        
        if grobner_match:
            # Split captured group into lines and clean
            grobner_basis = [
                line.strip().replace('^', '**') 
                for line in grobner_match.group(1).split('\n') 
                if line.strip()
            ]
        else:
            grobner_basis = []

    os.remove('temp_grobner.txt')

    if not grobner_basis:
        raise ValueError("No Gröbner basis found in the projection output.")

    return generate_dataset_from_grobner_basis(grobner_basis, variables, observed_constants, measured_derivatives,
                                              constant_data, derivative_data, region, num_points, seed)

def run_noiseless_system_data_generation(input_file, output_file):
    
    # Parse system data from original file
    parsed_system = parse_system_data(input_file)
    
    try:
        dataset = generate_data_for_system(
            equations=parsed_system['equations'],
            observed_constants=parsed_system['constants'],
            measured_derivatives=parsed_system['derivatives'],
            num_points=1000,
            seed=42,
            constant_data=True,
            derivative_data=True,
            region=[1, 10]
        )
        
        # Check if dataset has sufficient points
        if len(dataset) == 0 or len(next(iter(dataset.values()))) < 100:  # Check first array's length
            print("Dataset too small, generation failed")
            return False
            
        print("Successfully generated dataset with shapes:")
        for key, value in dataset.items():
            print(f"{key}: {value.shape}")
        
        save_dataset(dataset, output_file, delimiter=' ')
        
        return True  # Success
        
    except ValueError as e:
        print(f"Error generating dataset: {e}")
        return False

def add_gaussian_noise(input_file, output_file, epsilon, num_constant_cols=0):
    """
    Adds Gaussian noise to each non-constant column of a .dat file.
    
    Parameters:
        input_file (str): Path to input .dat file
        output_file (str): Path to save noisy output file
        epsilon (float): Percentage of column mean to use as noise standard deviation
        num_constant_cols (int): Number of constant columns at the start (default 0)
    """
    # Load data
    data = np.loadtxt(input_file)
    
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)  # Handle single-column case
    
    # Process each non-constant column
    for col in range(num_constant_cols, data.shape[1]):
        col_mean = np.abs(np.mean(data[:, col]))
        noise_std = col_mean * epsilon
        
        # Generate Gaussian noise for this column
        noise = np.random.normal(0, noise_std, data.shape[0])
        
        # Apply noise to column
        data[:, col] += noise
    
    # Save to output file with full precision
    np.savetxt(output_file, data, fmt='%.18e')

def run_noisy_system_data_generation(input_txt_file, input_dat_file, epsilon_list=None):
    """
    Generates multiple noisy versions of a .dat file with different noise levels.
    
    Parameters:
        input_txt_file (str): Path to input .txt file
        input_dat_file (str): Path to input .dat file
        epsilon_list (list): List of noise levels (default: [1e-3, 1e-2, 5e-2, 1e-1])
    
    Returns:
        list: Paths to generated noisy files
    """
    if epsilon_list is None:
        epsilon_list = [1e-3, 1e-2, 5e-2, 1e-1]
    
    generated_files = []
    
    # Get number of constant columns from the system file
    num_constant_cols = 0
    if os.path.exists(input_txt_file):
        try:
            parsed_system = parse_system_data(input_txt_file)
            num_constant_cols = len(parsed_system['constants'])
            print(f"Detected {num_constant_cols} constant columns")
        except Exception as e:
            print(f"Couldn't parse system file, assuming no constant columns: {str(e)}")
    
    # Read the data file with header
    try:
        # First read just the header to get column names
        with open(input_dat_file, 'r') as f:
            header = f.readline().strip().split()
        
        # Then read the numerical data (skiprows=1 to skip header)
        data = np.loadtxt(input_dat_file, skiprows=1)
        
        if len(data.shape) == 1:
            data = data.reshape(-1, 1)  # Handle single-column case
        
    except Exception as e:
        print(f"Failed to read data file: {str(e)}")
        return False
    
    # Generate noisy versions
    for epsilon in epsilon_list:
        output_file = input_dat_file.replace('.dat', f'_noisy_{epsilon}.dat').replace('+', '')
        try:
            noisy_data = data.copy()
            
            # Add noise to non-constant columns
            for col in range(num_constant_cols, noisy_data.shape[1]):
                col_mean = np.abs(np.mean(noisy_data[:, col]))
                noise_std = col_mean * epsilon
                noise = np.random.normal(0, noise_std, noisy_data.shape[0])
                noisy_data[:, col] += noise
            
            # Save with original header
            with open(output_file, 'w') as f:
                # Write header
                f.write('\t'.join(header) + '\n')
                # Write data
                np.savetxt(f, noisy_data, fmt='%.18e', delimiter='\t')
            
            generated_files.append(output_file)
            print(f"Successfully created noisy file with ε={epsilon}: {output_file}")
            
        except Exception as e:
            print(f"Failed to create noisy file for ε={epsilon}: {str(e)}")
    
    return True