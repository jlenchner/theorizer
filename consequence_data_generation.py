import os
import re
import numpy as np
from sympy import symbols, sympify, Poly, solve
from collections import defaultdict

def evaluate_polynomial(target_polynomial, dataset, observed_constants, measured_derivatives, measured_variables):
    """
    Evaluates the polynomial for all rows in the dataset.
    
    Args:
        target_polynomial: String representation of the polynomial
        dataset: NumPy array containing the data
        observed_constants: List of constant variable names
        measured_derivatives: List of derivative variable names
        measured_variables: List of measured variable names
        
    Returns:
        NumPy array of evaluated polynomial values
    """
    evaluated_values = []
    
    # Combine all variable names in the correct order
    if '' not in measured_derivatives:
        variable_order = observed_constants + measured_derivatives + measured_variables
    else:
        variable_order = observed_constants + measured_variables
    
    print("Variables in evaluation: ", variable_order)
    
    # Convert polynomial string to sympy expression once (more efficient)
    try:
        expr = sympify(target_polynomial.replace('^', '**'))
    except Exception as e:
        print(f"Polynomial parsing error: {e}")
        return np.full(dataset.shape[0], np.nan)
    
    for i in range(dataset.shape[0]):
        # Create substitution dictionary for this row
        substitutions = dict(zip(symbols(' '.join(variable_order)), dataset[i]))
        
        try:
            # Substitute and evaluate
            evaluated = expr.subs(substitutions)
            evaluated_values.append(float(evaluated.evalf()))
        except Exception as e:
            evaluated_values.append(np.nan)
            print(f"Evaluation error at row {i}: {e}")
    
    return np.array(evaluated_values)

def extract_info_from_file(filepath):
    with open(filepath, 'r') as file:
        content = file.readlines()
    
    # Initialize variables
    target_polynomial = []
    in_polynomial = False
    measured_variables = []
    observed_constants = []
    measured_derivatives = []
    
    for line in content:
        line = line.strip()
        if line.startswith('Measured Variables:'):
            measured_variables = eval(line.split('Measured Variables:')[1].strip())
        elif line.startswith('Observed Constants:'):
            observed_constants = eval(line.split('Observed Constants:')[1].strip())
        elif line.startswith('Measured Derivatives:'):
            measured_derivatives = eval(line.split('Measured Derivatives:')[1].strip())
        elif line.startswith('Target Polynomial:'):
            in_polynomial = True
        elif in_polynomial:
            if line:  # Only add non-empty lines
                target_polynomial.append(line)
    
    # Join polynomial lines with spaces (removes newlines but preserves structure)
    target_polynomial = ' '.join(target_polynomial)
    
    return target_polynomial, measured_variables, observed_constants, measured_derivatives

def generate_dataset(target_polynomial, measured_variables, observed_constants, measured_derivatives, constant_data=True, derivative_data=True, region=[1,5]):
    # Initialize dataset as a 2D array with 1000 rows and 0 columns
    dataset = np.zeros((1000, 0))  # Start with an empty 2D array

    # Step 2a: Generate data for constants
    if constant_data and observed_constants:  # Only proceed if constants are provided
        constant_values = {const: np.random.uniform(0, 10) for const in observed_constants}  # Sample once
        if 'c' in observed_constants:
            constant_values['c'] = 1
        if 'G' in observed_constants:
            constant_values['G'] = 1
        constant_data_matrix = np.array([[constant_values[const]] * 1000 for const in observed_constants]).T  # Repeat 1000 times
        dataset = np.hstack((dataset, constant_data_matrix))  # Append constant data
    elif not constant_data and observed_constants:
        constant_values = {const: np.random.uniform(region[0], region[1], 1000) for const in observed_constants}  # Generate 1000 samples per constant
        constant_data_matrix = np.column_stack([constant_values[const] for const in observed_constants])  # Stack properly
        dataset = np.hstack((dataset, constant_data_matrix))  # Append constant data

    derivative_dependent_variable = []
    derivative_data_matrix = []
    
    # Step 2b: Generate data for derivatives only if needed
    if derivative_data:
        if 'd1' in measured_variables and ('dx1dt' in measured_derivatives or 'd2x1dt2' in measured_derivatives):
            d1_data = np.random.uniform(region[0], region[1], 1000)
            dx1dt_data = np.diff(d1_data, append=d1_data[0]) if 'dx1dt' in measured_derivatives else None
            if 'dx1dt' in measured_derivatives:
                d2x1dt2_data = np.diff(dx1dt_data, append=dx1dt_data[0]) if 'd2x1dt2' in measured_derivatives else None
            else: 
                d2x1dt2_data = np.random.uniform(region[0],region[1], 1000) if 'd2x1dt2' in measured_derivatives else None
            
            if dx1dt_data is not None and 'dx1dt' in measured_derivatives:
                derivative_data_matrix.append(dx1dt_data)
            if d2x1dt2_data is not None and 'd2x1dt2' in measured_derivatives:
                derivative_data_matrix.append(d2x1dt2_data)

        elif 'd1' not in measured_variables and ('dx1dt' in measured_derivatives or 'd2x1dt2' in measured_derivatives):
            dx1dt_data = np.random.uniform(region[0],region[1], 1000)
            if 'dx1dt' in measured_derivatives:
                d2x1dt2_data = np.diff(dx1dt_data, append=dx1dt_data[0]) if 'd2x1dt2' in measured_derivatives else None
            else: 
                d2x1dt2_data = np.random.uniform(region[0],region[1], 1000) if 'd2x1dt2' in measured_derivatives else None
            if dx1dt_data is not None and 'dx1dt' in measured_derivatives:
                derivative_data_matrix.append(dx1dt_data)
            if d2x1dt2_data is not None and 'd2x1dt2' in measured_derivatives:
                derivative_data_matrix.append(d2x1dt2_data)
        
        if 'd2' in measured_variables and ('dx2dt' in measured_derivatives or 'd2x2dt2' in measured_derivatives):
            d2_data = np.random.uniform(region[0], region[1], 1000)
            dx2dt_data = np.diff(d2_data, append=d2_data[0]) if 'dx2dt' in measured_derivatives else None
            if 'dx2dt' in measured_derivatives:
                d2x2dt2_data = np.diff(dx2dt_data, append=dx2dt_data[0]) if 'd2x2dt2' in measured_derivatives else None
            else: 
                d2x2dt2_data = np.random.uniform(region[0],region[1], 1000) if 'd2x2dt2' in measured_derivatives else None
            
            if dx2dt_data is not None and 'dx2dt' in measured_derivatives:
                derivative_data_matrix.append(dx2dt_data)
            if d2x2dt2_data is not None and 'd2x2dt2' in measured_derivatives:
                derivative_data_matrix.append(d2x2dt2_data)

        elif 'd2' not in measured_variables and ('dx2dt' in measured_derivatives or 'd2x2dt2' in measured_derivatives):
            dx2dt_data = np.random.uniform(region[0],region[1], 1000)
            if 'dx2dt' in measured_derivatives:
                d2x2dt2_data = np.diff(dx2dt_data, append=dx2dt_data[0]) if 'd2x2dt2' in measured_derivatives else None
            else: 
                d2x2dt2_data = np.random.uniform(region[0],region[1], 1000) if 'd2x2dt2' in measured_derivatives else None

            if dx2dt_data is not None and 'dx2dt' in measured_derivatives:
                derivative_data_matrix.append(dx2dt_data)
            if d2x2dt2_data is not None and 'd2x2dt2' in measured_derivatives:
                derivative_data_matrix.append(d2x2dt2_data)

        if 'd1' in measured_variables and ('dx1dt' in measured_derivatives or 'd2x1dt2' in measured_derivatives):
            derivative_data_matrix.append(d1_data)
            measured_variables.remove('d1')
            derivative_dependent_variable.append('d1')
        if 'd2' in measured_variables and ('dx2dt' in measured_derivatives or 'd2x2dt2' in measured_derivatives):
            derivative_data_matrix.append(d2_data)
            measured_variables.remove('d2')
            derivative_dependent_variable.append('d2')

    if derivative_data_matrix:
        dataset = np.hstack((dataset, np.column_stack(derivative_data_matrix)))
    
    # Step 2c: Sort measured variables by degree in target polynomial
    degree_dict = defaultdict(int)
    for var in measured_variables:
        matches = re.findall(rf'{var}\^(\d+)', target_polynomial)
        if matches:
            degree_dict[var] = sum(int(exp) for exp in matches)
        else:
            degree_dict[var] = 1
    sorted_variables = measured_variables #sorted(measured_variables, key=lambda x: degree_dict[x], reverse=True)
    print(sorted_variables)
    # Generate data for all but the last variable
    for var in sorted_variables[:-1]:
        var_data = np.random.uniform(region[0], region[1], 1000)
        dataset = np.hstack((dataset, var_data.reshape(-1, 1)))

    # Step 2d: Generate data for the last variable
    last_var = sorted_variables[-1]
    last_var_symbol = symbols(last_var)  # Convert to symbolic variable

    roots_column = np.zeros((dataset.shape[0], 1)) * np.nan  # Initialize with NaN

    complex_count = 0
    evaluation_fail_count = 0

    for i in range(1000):
        # Substitute all known values into the polynomial at once
        polynomial = target_polynomial

        for j, var in enumerate(observed_constants + measured_derivatives + derivative_dependent_variable + sorted_variables[:-1]):
            # Replace variables with placeholders to avoid invalid syntax
            polynomial = re.sub(rf'\b{var}\b', f'{var}_value', polynomial)

        # Replace '^' with '**' for exponentiation
        polynomial = polynomial.replace('^', '**')

        # Convert the polynomial into a symbolic expression
        try:
            expr = sympify(polynomial)
            
            # Substitute numeric values into the expression
            substitutions = {}
            for j, var in enumerate(observed_constants + measured_derivatives + derivative_dependent_variable + sorted_variables[:-1]):
                substitutions[symbols(f'{var}_value')] = dataset[i, j]
            
            expr = expr.subs(substitutions)
            
            # Extract coefficients of the polynomial in the last variable
            poly = Poly(expr, last_var_symbol)
            coefficients = poly.all_coeffs()

            # Solve for the roots
            roots = np.roots(coefficients)
            
            # Filter real roots
            real_roots = roots[np.isreal(roots)].real
            
            if len(real_roots) > 0:
                roots_column[i] = real_roots[0]  # Use the first real root
            else:
                complex_count += 1
        except Exception as e:
            evaluation_fail_count += 1
            print(f"An exception occurred: {e}")  # Skip if polynomial evaluation fails   

    # Append the roots column to the dataset
    dataset = np.hstack((dataset, roots_column))

    # Remove rows with NaN values
    dataset = dataset[~np.isnan(dataset).any(axis=1)] 
    print("Complex roots thrown out: ", complex_count)
    print("Number of failed evaluations: ", evaluation_fail_count)  
    
    return dataset

def generate_dataset_no_der(target_polynomial, measured_variables, observed_constants, constant_data=True, region=[1,5]):
    # Initialize dataset as a 2D array with 1000 rows and 0 columns
    dataset = np.zeros((1000, 0))  # Start with an empty 2D array

    # Step 2a: Generate data for constants
    if constant_data and observed_constants:  # Only proceed if constants are provided
        constant_values = {const: np.random.uniform(0, 10) for const in observed_constants}  # Sample once
        if 'c' in observed_constants:
            constant_values['c'] = 1
        if 'G' in observed_constants:
            constant_values['G'] = 1
        constant_data_matrix = np.array([[constant_values[const]] * 1000 for const in observed_constants]).T  # Repeat 1000 times
        dataset = np.hstack((dataset, constant_data_matrix))  # Append constant data
    elif not constant_data and observed_constants:
        constant_values = {const: np.random.uniform(region[0], region[1], 1000) for const in observed_constants}  # Generate 1000 samples per constant
        constant_data_matrix = np.column_stack([constant_values[const] for const in observed_constants])  # Stack properly
        dataset = np.hstack((dataset, constant_data_matrix))  # Append constant data
    
    # Step 2c: Sort measured variables by degree in target polynomial
    degree_dict = defaultdict(int)
    for var in measured_variables:
        matches = re.findall(rf'{var}\^(\d+)', target_polynomial)
        if matches:
            degree_dict[var] = sum(int(exp) for exp in matches)
        else:
            degree_dict[var] = 1
    sorted_variables = measured_variables #sorted(measured_variables, key=lambda x: degree_dict[x], reverse=True)

    # Generate data for all but the last variable
    for var in sorted_variables[:-1]:
        var_data = np.random.uniform(region[0], region[1], 1000)
        dataset = np.hstack((dataset, var_data.reshape(-1, 1)))

    # Step 2d: Generate data for the last variable
    last_var = sorted_variables[-1]
    last_var_symbol = symbols(last_var)  # Convert to symbolic variable

    roots_column = np.zeros((dataset.shape[0], 1)) * np.nan  # Initialize with NaN

    complex_count = 0
    evaluation_fail_count = 0

    for i in range(1000):
        # Substitute all known values into the polynomial at once
        polynomial = target_polynomial
        for j, var in enumerate(observed_constants + sorted_variables[:-1]):
            # Replace variables with placeholders to avoid invalid syntax
            polynomial = re.sub(rf'\b{var}\b', f'{var}_value', polynomial)

        # Replace '^' with '**' for exponentiation
        polynomial = polynomial.replace('^', '**')

        # Convert the polynomial into a symbolic expression
        try:
            expr = sympify(polynomial)
            
            # Substitute numeric values into the expression
            substitutions = {}
            for j, var in enumerate(observed_constants + sorted_variables[:-1]):
                substitutions[symbols(f'{var}_value')] = dataset[i, j]
            
            expr = expr.subs(substitutions)
            
            # Extract coefficients of the polynomial in the last variable
            poly = Poly(expr, last_var_symbol)
            coefficients = poly.all_coeffs()

            # Solve for the roots
            roots = np.roots(coefficients)
            
            # Filter real roots
            real_roots = roots[np.isreal(roots)].real
            
            if len(real_roots) > 0:
                roots_column[i] = real_roots[0]  # Use the first real root
            else:
                complex_count += 1
        except Exception as e:
            evaluation_fail_count += 1
            print(f"An exception occurred: {e}")  # Skip if polynomial evaluation fails   

    # Append the roots column to the dataset
    dataset = np.hstack((dataset, roots_column))

    # Remove rows with NaN values
    dataset = dataset[~np.isnan(dataset).any(axis=1)] 
    print("Complex roots thrown out: ", complex_count)
    print("Number of failed evaluations: ", evaluation_fail_count)  
    
    return dataset

def add_gaussian_noise(input_file, output_file, epsilon):
    # Load data
    data = np.loadtxt(input_file)
        
    # Generate Gaussian noise with standard deviation as 5% of the column's mean
    noise = np.random.normal(0, np.abs(np.mean(data[:,-1])) * epsilon, len(data[:,-1]))
    
    # Apply noise only to non-constant columns
    
    data[:, -1] += noise
    
    # Save to output file
    np.savetxt(output_file, data, fmt='%.18e')

def run_consequence_noiseless_data_generation(input_file, output_file, region=[1, 10]):
    np.random.seed(42)

    # Extract info from file
    target_polynomial, measured_variables, observed_constants, measured_derivatives = extract_info_from_file(input_file)
    # Create a copy of the original file content
    with open(input_file, 'r') as f:
        original_content = f.read()
    
    # Identify derivative-dependent variables
    derivative_dependent_vars = []
    if 'dx1dt' in measured_derivatives or 'd2x1dt2' in measured_derivatives:
        if 'd1' in measured_variables:
            derivative_dependent_vars.append('d1')
    if 'dx2dt' in measured_derivatives or 'd2x2dt2' in measured_derivatives:
        if 'd2' in measured_variables:
            derivative_dependent_vars.append('d2')
    
    # Reorder measured variables if needed
    if derivative_dependent_vars:
        new_measured_vars = derivative_dependent_vars.copy()
        for var in measured_variables:
            if var not in derivative_dependent_vars:
                new_measured_vars.append(var)
        
        # Update the content with reordered measured variables
        updated_content = original_content.replace(
            f"Measured Variables: {measured_variables}",
            f"Measured Variables: {new_measured_vars}"
        )
        
        with open(input_file, 'w') as f:
            f.write(updated_content)
        
        measured_variables = new_measured_vars.copy()
    
    # Handle empty derivatives case
    if not measured_derivatives:
        measured_derivatives = ['']
    
    if '' not in measured_derivatives:
        dataset = generate_dataset(target_polynomial, measured_variables.copy(), observed_constants, measured_derivatives, True, True, region)
    else:
        dataset = generate_dataset_no_der(target_polynomial, measured_variables.copy(), observed_constants, True, region)

    if len(dataset) == 0:
        print("Dataset Generation Failed. Skipping \n")
        return False
    
    # Save dataset to file
    np.savetxt(output_file, dataset, delimiter=' ')
    print("Reordered measured variables:", [v for v in derivative_dependent_vars if v in measured_variables] + 
          [v for v in measured_variables if v not in derivative_dependent_vars])
    print("Evaluation: ", sum(evaluate_polynomial(target_polynomial, dataset[0:5], observed_constants, measured_derivatives, measured_variables)))
    print(f"Generated dataset for system with shape {dataset.shape} \n")

    return True

def run_consequence_noisy_data_generation(input_file):
    epsilon_list = [1e-3,1e-2,1e-1,5e-2]

    for epsilon in epsilon_list:
        output_file = input_file.replace('.dat', f'_{epsilon}.dat')
        add_gaussian_noise(input_file, output_file,epsilon)