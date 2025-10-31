import os
import re
import numpy as np
from sympy import symbols, sympify, Poly, solve
from collections import defaultdict

# ---------------------------------------------------------------------
# Utility: evaluate a target polynomial on a set of rows (for sanity checks)
# ---------------------------------------------------------------------
def evaluate_polynomial(target_polynomial, dataset, observed_constants, measured_derivatives, measured_variables):
    """
    Evaluate the target polynomial on the first axis (rows) of `dataset`.

    Args:
        target_polynomial: The target polynomial as a string; '^' is allowed
                           and will be converted to Python '**'.
        dataset: 2D array where columns are ordered as
                 [constants] + [derivatives if any] + [bases used by derivatives (d1/d2) if any] + [measured vars].
                 (This matches how the dataset is constructed below.)
        observed_constants: Names of constants (in the order they appear in dataset).
        measured_derivatives: Names of derivatives (ordered). Use [''] if none.
        measured_variables: Names of measured variables (ordered).

    Returns:
        1D array of float evaluations (nan where evaluation failed).
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

# ---------------------------------------------------------------------
# Parsing of the consequence file emitted by your system
# ---------------------------------------------------------------------
def extract_info_from_file(filepath):
    """
    Extract target polynomial and variable lists from a consequence file.

    The file is expected to contain lines like:
      Measured Variables: [...]
      Observed Constants: [...]
      Measured Derivatives: [...]
      Target Polynomial:
      <one or more lines of the polynomial>

    Args:
        filepath: Path to the consequence file.

    Returns:
        (target_polynomial, measured_variables, observed_constants, measured_derivatives)
        where target_polynomial is a single-line string with '^' preserved.
    """
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


# ---------------------------------------------------------------------
# Dataset generators (with and without derivatives)
# ---------------------------------------------------------------------
def generate_dataset(target_polynomial,measured_variables,observed_constants,measured_derivatives,constant_data=True,derivative_data=True,region=[1, 5]):
    """
    Generate a dataset when there are no derivatives.

    Column order produced:
        [constants...] +
        [measured_variables except the last solved variable, in order] +
        [last variable solved from the polynomial]

    Args:
        target_polynomial: Target polynomial as a string; '^' allowed.
        measured_variables: Measured variable names (may include theta/trig/exponential).
        observed_constants: Constant names.
        constant_data: If True, constants are fixed across rows; else per-row sampling.
        region: [low, high] sampling interval for random draws.

    Returns:
        2D array where final column is the solved variable. Rows with complex/no solution are dropped.
    """
    N = 1000
    dataset = np.zeros((N, 0))  # start empty

    # ---- Constants ----
    const_data_cols = []
    if observed_constants:
        if constant_data:
            # One draw per constant, held fixed across N rows
            constant_values = {const: np.random.uniform(0, 10) for const in observed_constants}
            if 'c' in observed_constants:
                constant_values['c'] = 1
            if 'G' in observed_constants:
                constant_values['G'] = 1
            for const in observed_constants:
                const_data_cols.append(np.full(N, constant_values[const]))
        else:
            # Independent draws per row
            for const in observed_constants:
                const_data_cols.append(np.random.uniform(region[0], region[1], N))
    if const_data_cols:
        dataset = np.hstack([dataset] + [col.reshape(-1, 1) for col in const_data_cols])

    # ---- Derivatives (and their bases d1/d2 when needed) ----
    derivative_dependent_variable = []  # order of bases (subset of ['d1','d2']) used in derivative construction
    derivative_arrays = {}              

    # d1-branch
    if derivative_data and ('dx1dt' in measured_derivatives or 'd2x1dt2' in measured_derivatives):
        if 'd1' in measured_variables:
            d1_data      = np.random.uniform(region[0], region[1], N)
            dx1dt_data   = np.diff(d1_data, append=d1_data[0]) if 'dx1dt' in measured_derivatives else None
            if 'd2x1dt2' in measured_derivatives:
                if dx1dt_data is not None:
                    d2x1dt2_data = np.diff(dx1dt_data, append=dx1dt_data[0])
                else:
                    d2x1dt2_data = np.random.uniform(region[0], region[1], N)
            else:
                d2x1dt2_data = None

            if 'dx1dt'   in measured_derivatives and dx1dt_data is not None: derivative_arrays['dx1dt']   = dx1dt_data
            if 'd2x1dt2' in measured_derivatives and d2x1dt2_data is not None: derivative_arrays['d2x1dt2'] = d2x1dt2_data

            # mark d1 as derivative-dependent and remove from measured list (preserve order elsewhere)
            derivative_dependent_variable.append('d1')
            measured_variables = [v for v in measured_variables if v != 'd1']
        else:
            # No d1 measured; synthesize derivatives directly
            if 'dx1dt' in measured_derivatives:
                dx1dt_data = np.random.uniform(region[0], region[1], N)
                derivative_arrays['dx1dt'] = dx1dt_data
                if 'd2x1dt2' in measured_derivatives:
                    derivative_arrays['d2x1dt2'] = np.diff(dx1dt_data, append=dx1dt_data[0])
            elif 'd2x1dt2' in measured_derivatives:
                derivative_arrays['d2x1dt2'] = np.random.uniform(region[0], region[1], N)

    # d2-branch
    if derivative_data and ('dx2dt' in measured_derivatives or 'd2x2dt2' in measured_derivatives):
        if 'd2' in measured_variables:
            d2_data      = np.random.uniform(region[0], region[1], N)
            dx2dt_data   = np.diff(d2_data, append=d2_data[0]) if 'dx2dt' in measured_derivatives else None
            if 'd2x2dt2' in measured_derivatives:
                if dx2dt_data is not None:
                    d2x2dt2_data = np.diff(dx2dt_data, append=dx2dt_data[0])
                else:
                    d2x2dt2_data = np.random.uniform(region[0], region[1], N)
            else:
                d2x2dt2_data = None

            if 'dx2dt'   in measured_derivatives and dx2dt_data is not None: derivative_arrays['dx2dt']   = dx2dt_data
            if 'd2x2dt2' in measured_derivatives and d2x2dt2_data is not None: derivative_arrays['d2x2dt2'] = d2x2dt2_data

            # mark d1 as derivative-dependent and remove from measured list (preserve order elsewhere)
            derivative_dependent_variable.append('d2')
            measured_variables = [v for v in measured_variables if v != 'd2']
        else:
            # No d1 measured; synthesize derivatives directly
            if 'dx2dt' in measured_derivatives:
                dx2dt_data = np.random.uniform(region[0], region[1], N)
                derivative_arrays['dx2dt'] = dx2dt_data
                if 'd2x2dt2' in measured_derivatives:
                    derivative_arrays['d2x2dt2'] = np.diff(dx2dt_data, append=dx2dt_data[0])
            elif 'd2x2dt2' in measured_derivatives:
                derivative_arrays['d2x2dt2'] = np.random.uniform(region[0], region[1], N)

    # Append derivative columns in the declared order; then append their bases (if used)
    if measured_derivatives and '' not in measured_derivatives:
        for dname in measured_derivatives:
            if dname in derivative_arrays:
                dataset = np.hstack((dataset, derivative_arrays[dname].reshape(-1, 1)))

        # then append d1/d2 bases (if any) in that same order
        for base in derivative_dependent_variable:
            if base == 'd1':
                dataset = np.hstack((dataset, d1_data.reshape(-1, 1)))
            elif base == 'd2':
                dataset = np.hstack((dataset, d2_data.reshape(-1, 1)))

    # ---- theta / trig/exponential variables ----
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    trig_vars_present = [v for v in measured_variables if v in trig_funcs]
    need_theta_stream = ('theta' in measured_variables) or bool(trig_vars_present)

    theta_stream = None
    if need_theta_stream:
        theta_stream = np.random.uniform(region[0], region[1], N)

    trig_values = {}
    for v in trig_vars_present:
        trig_values[v] = trig_funcs[v](theta_stream)

     # ---- choose last var to solve for (exclude specials) ----
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    # Preserve measured order; pick last that isn't a special
    solve_candidates = [v for v in measured_variables if v not in specials]
    if len(solve_candidates) == 0:
        # Nothing left to solve for; signal failure
        return np.zeros((0, 0))

    last_var = solve_candidates[-1]
    known_measured_vars = [v for v in measured_variables if v != last_var]

    # ---- sample measured vars (except the solved one), in order ----
    measured_cols = []
    for var in known_measured_vars:
        if var == 'theta':
            measured_cols.append(theta_stream)
        elif var in trig_values:
            measured_cols.append(trig_values[var])
        else:
            measured_cols.append(np.random.uniform(region[0], region[1], N))

    if measured_cols:
        dataset = np.hstack([dataset] + [col.reshape(-1, 1) for col in measured_cols])

    # ---- solve for last_var per-row (drop rows with complex/no real solution) ----
    last_var_symbol = symbols(last_var)
    roots_column = np.full((N, 1), np.nan)

    complex_count = 0
    evaluation_fail_count = 0

    # Substitution order must mirror dataset columns:
    sub_order = []
    sub_order.extend(observed_constants)
    if measured_derivatives and '' not in measured_derivatives:
        sub_order.extend(measured_derivatives)
        sub_order.extend(derivative_dependent_variable) # bases were appended after derivatives
    sub_order.extend(known_measured_vars)

    for i in range(N):
        polynomial = target_polynomial

        # Replace variable names with *_value placeholders so we can substitute scalars
        for name in (observed_constants +
                     (measured_derivatives if (measured_derivatives and '' not in measured_derivatives) else []) +
                     derivative_dependent_variable +
                     known_measured_vars):
            polynomial = re.sub(rf'\b{re.escape(name)}\b', f'{name}_value', polynomial)

        polynomial = polynomial.replace('^', '**')

        try:
            expr = sympify(polynomial)
            substitutions = {}
            for j, var_name in enumerate(sub_order):
                substitutions[symbols(f'{var_name}_value')] = dataset[i, j]

            expr = expr.subs(substitutions)

            poly = Poly(expr, last_var_symbol)
            coeffs = poly.all_coeffs()
            roots = np.roots(coeffs)
            real_roots = roots[np.isreal(roots)].real
            if len(real_roots) > 0:
                roots_column[i, 0] = real_roots[0]
            else:
                complex_count += 1
        except Exception as e:
            evaluation_fail_count += 1
            # print(f"An exception occurred: {e}")
            continue

    # Append solved column and drop rows with NaN (failed/complex)
    dataset = np.hstack((dataset, roots_column))
    dataset = dataset[~np.isnan(dataset).any(axis=1)]

    print("Complex roots thrown out: ", complex_count)
    print("Number of failed evaluations: ", evaluation_fail_count)

    return dataset

def generate_dataset_no_der(target_polynomial,measured_variables,observed_constants,constant_data=True,region=[1, 5]):
    """
    Generate a dataset when there are no derivatives.

    Column order produced:
        [constants...] +
        [measured_variables except the last solved variable, in order] +
        [last variable solved from the polynomial]

    Args:
        target_polynomial: Target polynomial as a string; '^' allowed.
        measured_variables: Measured variable names (may include theta/trig/exponential).
        observed_constants: Constant names.
        constant_data: If True, constants are fixed across rows; else per-row sampling.
        region: [low, high] sampling interval for random draws.

    Returns:
        2D array where final column is the solved variable. Rows with complex/no solution are dropped.
    """
    N = 1000
    dataset = np.zeros((N, 0))

    # ---- Constants ----
    const_data_cols = []
    if observed_constants:
        if constant_data:
            constant_values = {const: np.random.uniform(0, 10) for const in observed_constants}
            if 'c' in observed_constants:
                constant_values['c'] = 1
            if 'G' in observed_constants:
                constant_values['G'] = 1
            for const in observed_constants:
                const_data_cols.append(np.full(N, constant_values[const]))
        else:
            for const in observed_constants:
                const_data_cols.append(np.random.uniform(region[0], region[1], N))
    if const_data_cols:
        dataset = np.hstack([dataset] + [col.reshape(-1, 1) for col in const_data_cols])

    # ---- theta / trig/exponential generated variables ----
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    trig_vars_present = [v for v in measured_variables if v in trig_funcs]
    need_theta_stream = ('theta' in measured_variables) or bool(trig_vars_present)

    theta_stream = None
    if need_theta_stream:
        theta_stream = np.random.uniform(region[0], region[1], N)

    trig_values = {}
    for v in trig_vars_present:
        trig_values[v] = trig_funcs[v](theta_stream)

    # ---- choose last var to solve (exclude specials) ----
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    solve_candidates = [v for v in measured_variables if v not in specials]
    if len(solve_candidates) == 0:
        return np.zeros((0, 0))

    last_var = solve_candidates[-1]
    known_measured_vars = [v for v in measured_variables if v != last_var]

    # ---- sample measured vars (except the solved one) ----
    measured_cols = []
    for var in known_measured_vars:
        if var == 'theta':
            measured_cols.append(theta_stream)
        elif var in trig_values:
            measured_cols.append(trig_values[var])
        else:
            measured_cols.append(np.random.uniform(region[0], region[1], N))
    if measured_cols:
        dataset = np.hstack([dataset] + [col.reshape(-1, 1) for col in measured_cols])

    # ---- solve for last_var per-row ----
    last_var_symbol = symbols(last_var)
    roots_column = np.full((N, 1), np.nan)

    complex_count = 0
    evaluation_fail_count = 0

    # Sub order here is constants + known measured
    sub_order = []
    sub_order.extend(observed_constants)
    sub_order.extend(known_measured_vars)

    for i in range(N):
        polynomial = target_polynomial
        for name in (observed_constants + known_measured_vars):
            polynomial = re.sub(rf'\b{re.escape(name)}\b', f'{name}_value', polynomial)
        polynomial = polynomial.replace('^', '**')

        try:
            expr = sympify(polynomial)
            substitutions = {}
            for j, var_name in enumerate(sub_order):
                substitutions[symbols(f'{var_name}_value')] = dataset[i, j]

            expr = expr.subs(substitutions)
            poly = Poly(expr, last_var_symbol)
            coeffs = poly.all_coeffs()
            roots = np.roots(coeffs)
            real_roots = roots[np.isreal(roots)].real
            if len(real_roots) > 0:
                roots_column[i, 0] = real_roots[0]
            else:
                complex_count += 1
        except Exception as e:
            evaluation_fail_count += 1
            # print(f"An exception occurred: {e}")
            continue

    dataset = np.hstack((dataset, roots_column))
    dataset = dataset[~np.isnan(dataset).any(axis=1)]

    print("Complex roots thrown out: ", complex_count)
    print("Number of failed evaluations: ", evaluation_fail_count)

    return dataset

# ---------------------------------------------------------------------
# Noise and top-level runners
# ---------------------------------------------------------------------
def add_gaussian_noise(input_file, output_file, epsilon):
    """
    Add Gaussian noise to the last column of a saved dataset.

    IMPORTANT: This matches your original behavior (noise on the solved column only).

    Args:
        input_file: Path to the base (noiseless) .dat file.
        output_file: Path to write the noisy .dat file.
        epsilon: Scale factor for the std dev; sigma = |mean(last_col)| * epsilon.
    """
    data = np.loadtxt(input_file)        
    noise = np.random.normal(0, np.abs(np.mean(data[:,-1])) * epsilon, len(data[:,-1]))
    data[:, -1] += noise
    np.savetxt(output_file, data, fmt='%.18e')

def run_consequence_noiseless_data_generation(input_file, output_file, region=[1, 10], seed=42):
    """
    Generate a noiseless dataset for a single consequence specification file.

    Steps:
      1) Parse file to get target polynomial and variable lists.
      2) If derivatives depend on d1/d2, reorder 'Measured Variables' so those bases come first.
         (This update is written back to the same file to keep metadata consistent.)
      3) Generate dataset using derivative-aware or no-derivative path.
      4) Save to `output_file`.

    Args:
        input_file: Path to the consequence.txt-like file with sections described above.
        output_file: Destination .dat path.
        region: [low, high] sampling interval for random draws.
        seed: int for data generation seed. 

    Returns:
        True if generation succeeded and wrote a non-empty dataset; False otherwise.
    """
    np.random.seed(seed) # determinism for sampling

    target_polynomial, measured_variables, observed_constants, measured_derivatives = extract_info_from_file(input_file)
    with open(input_file, 'r') as f:
        original_content = f.read()
    
    # Detect whether derivatives “use” d1/d2, so we can move d1/d2 up front
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
    
    np.savetxt(output_file, dataset, delimiter=' ')
    print("Reordered measured variables:", [v for v in derivative_dependent_vars if v in measured_variables] + 
          [v for v in measured_variables if v not in derivative_dependent_vars])
    print("Evaluation: ", sum(evaluate_polynomial(target_polynomial, dataset[0:5], observed_constants, measured_derivatives, measured_variables)))
    print(f"Generated dataset for system with shape {dataset.shape} \n")

    return True

def run_consequence_noisy_data_generation(input_file):
    """
    Emit multiple noisy variants of a saved dataset.

    Args:
        input_file: Path to the noiseless `.dat` file previously created.

    Side effects:
        Writes files alongside `input_file` with suffixes:
        _0.001.dat, _0.01.dat, _0.1.dat, _0.05.dat
        (Order preserved from your original list.)
    """
    epsilon_list = [1e-3,1e-2,1e-1,5e-2]
    for epsilon in epsilon_list:
        output_file = input_file.replace('.dat', f'_{epsilon}.dat')
        add_gaussian_noise(input_file, output_file,epsilon)
