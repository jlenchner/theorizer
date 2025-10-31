import os
import re
import ast
from typing import Dict, List, Tuple, Optional

import numpy as np
import sympy as sp
from sympy import symbols, sympify, Poly

from m2_functions import projection

# Special variables that should be derived from theta, not solved numerically
SPECIALS = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
TRIG_FUNCS = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}

def parse_system_data(file_path):
    """
    Parse a system description file with Variables/Constants/Derivatives/Equations
    and their Units of Measure.

    Args:
        file_path: Path to the system text file.

    Returns:
        Dict with keys:
          - system_number: Optional[int]
          - variables, constants, derivatives: List[str]
          - equations: List[str] (power '^' normalized to '**', whitespace removed)
          - var_units, const_units, deriv_units: List[str]
          - eqn_units: List[str] (one per equation)
    """
    with open(file_path, 'r') as file:
        content = file.read()

    system_data = {}
    
    # Optional system number
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
    
    return system_data

def save_dataset(dataset, filename, delimiter='\t'):
    """
    Save a columnar dataset with a header row.

    Args:
        dataset: Mapping of column name -> 1D numpy array (all same length).
        filename: Output path (.dat or .txt).
        delimiter: Column separator (default: tab).

    Raises:
        ValueError if dataset is empty or arrays have inconsistent lengths.
    """
    if not dataset:
        raise ValueError("Dataset is empty")
    
    lengths = [len(arr) for arr in dataset.values()]
    if len(set(lengths)) > 1:
        raise ValueError("Inconsistent array lengths in dataset")
    
    headers = list(dataset.keys())
    data = np.column_stack([dataset[var] for var in headers])
    
    # Write to file
    with open(filename, 'w') as f:
        f.write(delimiter.join(headers) + '\n')
        np.savetxt(f, data, delimiter=delimiter, fmt='%.8e')

def get_variables_from_equations(equations):
    """Return the set of variables appearing in each equation string."""
    return [extract_variables(eq) for eq in equations]

def extract_variables(equation):
    """
    Extract variable-like tokens from a single equation string.

    Returns:
        Unique variable names (order not guaranteed).
    """
    variables = re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', equation)
    return list(set(variables))

def get_base_variable(derivative):
    """
    Infer base stream for a derivative name like 'dx1dt' or 'd2x2dt2'.

    Returns:
        'd1', 'd2', ... if recognizable, else None.
    """
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
    Collect all variables by scanning equations in order, deduping while preserving first appearance,
    then reverse the list.

    This yields an order where variables from equation 1 end up on the right (last),
    those from later equations trend left (earlier). Helpful for elimination heuristics.

    Args:
        equations: List of equation strings.

    Returns:
        Reversed deduplicated variable list.
    """
    all_vars = []    
    equation_variables = get_variables_from_equations(equations)
    for vars_list in equation_variables:
        for var in vars_list:
            if var not in all_vars:
                all_vars.append(var)
    
    all_vars.reverse()    
    return all_vars

def seed_theta_family(variables, N, region):
    """
    If 'theta' and/or {'sinTheta','cosTheta','eTheta'} appear, generate a single theta stream and
    deterministically compute auxiliaries. Only keys present in `variables` are returned.

    Args:
        variables: Full variable list for the dataset.
        N: Number of samples.
        region: (low, high) for uniform sampling of theta.

    Returns:
        Dict mapping present specials -> NumPy arrays.
    """
    present_aux = [v for v in variables if v in TRIG_FUNCS]
    need_theta = ('theta' in variables) or bool(present_aux)
    out = {}

    if not need_theta:
        return out

    theta_stream = np.random.uniform(region[0], region[1], N)

    if 'theta' in variables:
        out['theta'] = theta_stream

    for v in present_aux:
        out[v] = TRIG_FUNCS[v](theta_stream)

    return out



def evaluate_polynomials(dataset, equations, variable_order, tolerance=1e-5):
    """
    Check how well a dataset satisfies polynomial equations == 0.

    Args:
        dataset: Mapping var -> values (1D arrays).
        equations: A single equation or list of equations (strings). Each should evaluate to ~0.
        variable_order: The argument order expected by each lambdified equation.
        tolerance: Max absolute error to consider a row "valid".

    Returns:
        Dict with:
          total_points, valid_points, pass_rate (%),
          max_absolute_error, valid_indices (list of row indices).
    """
    if isinstance(equations, str):
        equations = [equations]

    sym_vars = sp.symbols(' '.join(variable_order))
    eval_fns = []
    for eq in equations:
        try:
            expr = sp.sympify(eq.replace('^', '**'))
            fn = sp.lambdify(sym_vars, expr, modules="numpy")
            eval_fns.append(fn)
        except Exception as e:
            raise ValueError(f"Could not parse equation '{eq}': {str(e)}")
    
    try:
        data_matrix = np.column_stack([dataset[var] for var in variable_order])
    except KeyError as e:
        raise KeyError(f"Variable {e} is missing from the dataset. Check your variable_order list.")
    
    total_points = data_matrix.shape[0]
    valid_mask = np.ones(total_points, dtype=bool)
    max_error = 0.0
    
    for fn in eval_fns:
        try:
            results = fn(*data_matrix.T)
            results = np.array(results)
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
                                        constant_data=True, derivative_data=True, region=[1, 10],
                                        num_points=10000, seed=42, tol=1e-6):
    """
    Sample a dataset consistent with a (lex) Gröbner basis by solving one variable per polynomial
    row-wise, while sampling the others from a region.

    Args:
        grobner_basis: List of polynomials as strings.
        variables: All variables (defines column order for sympy symbols).
        observed_constants: Constants to include; if constant_data=True, draw a single value and repeat.
        measured_derivatives: Derivatives to include; attempt finite difference if a base stream exists.
        constant_data: If True, constants are fixed across rows; else sampled per row.
        derivative_data: If True, include derivative streams.
        region: Sampling range (low, high) for non-constant/random vars.
        num_points: Desired number of rows before filtering by solvability.
        seed: RNG seed.
        tol: Unused here (kept for symmetry; could be used for filtering).

    Returns:
        Dict var -> array (all same length) on success, or None if no valid rows were found.
    """

    np.random.seed(seed)
    data = {}
    N = num_points

    # Symbols in the exact order of `variables`
    sym_vars = symbols(' '.join(variables))

    # Solve sparser equations earlier
    sorted_basis = sorted(grobner_basis, key=lambda eq: len(extract_variables(eq)))
    num_valid_points = N

    # 0) Constants
    for const in observed_constants:
        if constant_data:
            data[const] = np.full(num_valid_points, np.random.uniform(1,10))
        else:
            data[const] = np.random.uniform(region[0], region[1], num_valid_points)

    # 1) Derivatives (forward diff if base exists)
    if derivative_data:
        for deriv in measured_derivatives:
            base_var = get_base_variable(deriv)  # e.g., 'd1' or 'd2'
            if base_var in variables and base_var in data:
                base_data = data[base_var]
                data[deriv] = np.diff(base_data, append=base_data[0])
            elif base_var in variables and base_var not in data:
                data[deriv] = np.random.uniform(region[0], region[1], num_valid_points)
            else:
                data[deriv] = np.random.uniform(region[0], region[1], num_valid_points)

    # 2) Theta & auxiliaries
    data.update(seed_theta_family(variables, num_valid_points, region))

    # 3) Solve one variable per polynomial (others sampled)
    for eq in sorted_basis:
        eq_expr = sympify(eq)
        eq_variables = extract_variables(eq)

        unknown_vars = [v for v in eq_variables if v not in data]

        if not unknown_vars:
            continue

        # Never solve for theta/sin/cos/e^theta; they are computed
        non_special_unknowns = [v for v in unknown_vars if v not in SPECIALS]

        if not non_special_unknowns:
            continue

        # Sample all but one unknown uniformly; solve for the last.
        for var in non_special_unknowns[:-1]:
            data[var] = np.random.uniform(region[0], region[1], num_valid_points)

        last_var = non_special_unknowns[-1]
        last_var_sym = sym_vars[variables.index(last_var)]
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
                if var == last_var:
                    continue
                if var in data:
                    subs_dict[sym_vars[variables.index(var)]] = data[var][i]
            try:
                coeffs_i = [c.subs(subs_dict) for c in coefficients]
                coeffs_i = [complex(c) for c in coeffs_i]  # numeric
                roots = np.roots(coeffs_i)
                real_roots = roots[np.isreal(roots)].real
                if real_roots.size > 0:
                    last_var_data[i] = real_roots[0]
                    valid_indices.append(i)
            except Exception:
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

    # 4) Fill any remaining variables by sampling
    for var in variables:
        if var not in data:
            data[var] = np.random.uniform(region[0], region[1], num_valid_points)

    return data



def generate_data_for_system(equations, observed_constants, measured_derivatives, num_points=1000, seed=42, 
                             constant_data=True, derivative_data=True, region=[1, 10]):
    """
    Build a Gröbner basis (via Macaulay2 projection of the identity map) and
    generate a dataset consistent with it.

    Args:
        equations: Axiom equations (strings in M2/SymPy-parseable form).
        observed_constants: Constant names to include.
        measured_derivatives: Derivative names to include.
        num_points, seed, constant_data, derivative_data, region: Generation controls.

    Returns:
        Dict var -> ndarray dataset.

    Raises:
        ValueError if no basis lines are found in the M2 output.
    """

    variables = collect_variables_in_sequence_and_reverse(equations)

    # Eliminate nothing (project onto the same variables) to get a GB for the ideal.
    projection(variables, equations, variables, [], filename='temp_grobner.txt')

    with open('temp_grobner.txt', 'r') as file:
        content = file.read()
        
        # Extract lines under the labeled section
        grobner_match = re.search(
            r'Polynomials of the Gröbner basis of the eliminated ideal:\s*\n(.*?)(?=\n\S+:|$)',
            content, 
            re.DOTALL
        )
        
        if grobner_match:
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

def run_noiseless_system_data_generation(input_file, output_file, region = [1,10]):
    """
    End-to-end: parse a system file, generate a clean dataset, and save it to disk.

    Args:
        input_file: Path to system definition file.
        output_file: Where to write the .dat file (header + data).
        region: Sampling range for non-constant variables.

    Returns:
        True on success, False if generation failed or was too small.
    """
    
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
            region=region
        )
        
        if len(dataset) == 0 or len(next(iter(dataset.values()))) < 100:  
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
    Add Gaussian noise to each non-constant column of a .dat file.

    Args:
        input_file: Path to input .dat (numeric, with or without header).
        output_file: Destination .dat.
        epsilon: Noise scale; std = epsilon * |column mean|.
        num_constant_cols: Number of leading constant columns to keep untouched.
    """
    data = np.loadtxt(input_file)
    
    if len(data.shape) == 1:
        data = data.reshape(-1, 1)  
    
    for col in range(num_constant_cols, data.shape[1]):
        col_mean = np.abs(np.mean(data[:, col]))
        noise_std = col_mean * epsilon        
        noise = np.random.normal(0, noise_std, data.shape[0])        
        data[:, col] += noise    
    np.savetxt(output_file, data, fmt='%.18e')

def run_noisy_system_data_generation(input_txt_file, input_dat_file, epsilon_list=None):
    """
    Produce multiple noisy variants of a saved dataset, preserving the header.

    Args:
        input_txt_file: System file (used to detect number of constants).
        input_dat_file: Clean data file (header in first line).
        epsilon_list: Noise std scales (default: [1e-3, 1e-2, 5e-2, 1e-1]).

    Returns:
        True if at least one noisy file was created; False if data read failed.
    """
    if epsilon_list is None:
        epsilon_list = [1e-3, 1e-2, 5e-2, 1e-1]
    
    generated_files = []    
    num_constant_cols = 0
    if os.path.exists(input_txt_file):
        try:
            parsed_system = parse_system_data(input_txt_file)
            num_constant_cols = len(parsed_system['constants'])
            print(f"Detected {num_constant_cols} constant columns")
        except Exception as e:
            print(f"Couldn't parse system file, assuming no constant columns: {str(e)}")
    
    try:
        with open(input_dat_file, 'r') as f:
            header = f.readline().strip().split()
        
        data = np.loadtxt(input_dat_file, skiprows=1)        
        if len(data.shape) == 1:
            data = data.reshape(-1, 1)  # Handle single-column case
    except Exception as e:
        print(f"Failed to read data file: {str(e)}")
        return False
    
    for epsilon in epsilon_list:
        output_file = input_dat_file.replace('.dat', f'_noisy_{epsilon}.dat').replace('+', '')
        try:
            noisy_data = data.copy()            
            for col in range(num_constant_cols, noisy_data.shape[1]):
                col_mean = np.abs(np.mean(noisy_data[:, col]))
                noise_std = col_mean * epsilon
                noise = np.random.normal(0, noise_std, noisy_data.shape[0])
                noisy_data[:, col] += noise
            
            with open(output_file, 'w') as f:
                f.write('\t'.join(header) + '\n')
                np.savetxt(f, noisy_data, fmt='%.18e', delimiter='\t')
            
            generated_files.append(output_file)
            print(f"Successfully created noisy file with ε={epsilon}: {output_file}")
            
        except Exception as e:
            print(f"Failed to create noisy file for ε={epsilon}: {str(e)}")
    
    return True
