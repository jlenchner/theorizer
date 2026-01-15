import os
import re
import numpy as np
from sympy import symbols, sympify, Poly, solve, sqrt, Symbol
from collections import defaultdict
from scipy.integrate import solve_ivp

# Mapping from derivative names to their base variables and order
DERIVATIVE_INFO = {
    'dx1dt': ('d1', 1),
    'd2x1dt2': ('d1', 2),
    'dx2dt': ('d2', 1),
    'd2x2dt2': ('d2', 2),
}


def generate_dataset_algebraic_derivative(
    target_polynomial,
    measured_variables,
    observed_constants,
    measured_derivatives,
    constant_data=True,
    region=[1, 5],
    N=1000
):
    """
    Generate dataset by solving algebraically for the derivative from the polynomial.
    
    This is a cleaner fallback than finite differences because the data will exactly
    satisfy the polynomial. Works even when base variable is not in measured variables.
    
    The workflow:
    1. Sample constants (fixed values)
    2. Sample all measured variables
    3. Solve the polynomial algebraically for the derivative value at each row
    
    Args:
        target_polynomial: String representation of the consequence polynomial
        measured_variables: List of measured variable names
        observed_constants: List of constant names
        measured_derivatives: List containing exactly one derivative name
        constant_data: If True, constants are fixed; if False, sampled per row
        region: [min, max] range for sampling
        N: Number of data points to generate
    
    Returns:
        NumPy array with columns in order: constants, measured_variables, derivative
        Returns empty array if generation fails
    """
    if len(measured_derivatives) != 1 or measured_derivatives[0] == '':
        print("Algebraic derivative method requires exactly one derivative")
        return np.zeros((0, 0))
    
    deriv_name = measured_derivatives[0]
    deriv_symbol = Symbol(deriv_name)
    
    # -------------------------
    # 0) Setup constants
    # -------------------------
    constant_values = {}
    if observed_constants:
        if constant_data:
            for const in observed_constants:
                if const == 'c':
                    constant_values[const] = 1
                elif const == 'G':
                    constant_values[const] = 1
                else:
                    constant_values[const] = np.random.uniform(1, 10)
        else:
            for const in observed_constants:
                constant_values[const] = np.random.uniform(region[0], region[1])
    
    # -------------------------
    # 1) Handle theta and trig functions
    # -------------------------
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    
    theta_stream = None
    need_theta_stream = ('theta' in measured_variables) or any(v in measured_variables for v in trig_funcs)
    if need_theta_stream:
        theta_stream = np.random.uniform(region[0], region[1], N)
    
    # -------------------------
    # 2) Generate all measured variables (not the derivative)
    # -------------------------
    var_data = {}
    for var in measured_variables:
        if var == 'theta':
            var_data[var] = theta_stream
        elif var in trig_funcs:
            var_data[var] = trig_funcs[var](theta_stream)
        else:
            var_data[var] = np.random.uniform(region[0], region[1], N)
    
    # -------------------------
    # 3) Parse polynomial and solve for the derivative
    # -------------------------
    poly_str = target_polynomial.replace('^', '**')
    try:
        poly_expr = sympify(poly_str)
    except Exception as e:
        print(f"Failed to parse polynomial: {e}")
        return np.zeros((0, 0))
    
    # Create symbolic variables for substitution
    all_symbols = {}
    for const in observed_constants:
        all_symbols[const] = Symbol(const)
    for var in measured_variables:
        all_symbols[var] = Symbol(var)
    
    # -------------------------
    # 4) Solve for derivative at each row
    # -------------------------
    deriv_values = np.full(N, np.nan)
    valid_indices = []
    complex_count = 0
    fail_count = 0
    
    for i in range(N):
        # Build substitution dictionary
        subs_dict = {}
        for const, val in constant_values.items():
            subs_dict[all_symbols[const]] = val
        for var in measured_variables:
            subs_dict[all_symbols[var]] = var_data[var][i]
        
        try:
            # Substitute known values
            expr_substituted = poly_expr.subs(subs_dict)
            
            # Try to get polynomial coefficients and use np.roots (faster)
            try:
                poly = Poly(expr_substituted, deriv_symbol)
                coeffs = [complex(c) for c in poly.all_coeffs()]
                roots = np.roots(coeffs)
                real_roots = roots[np.isreal(roots)].real
                
                if len(real_roots) > 0:
                    # Take the root that's closest to a reasonable physical value
                    valid_roots = real_roots[(real_roots >= region[0] * 0.1) & (real_roots <= region[1] * 100)]
                    if len(valid_roots) > 0:
                        deriv_values[i] = valid_roots[0]
                        valid_indices.append(i)
                    else:
                        # Take any real root
                        deriv_values[i] = real_roots[0]
                        valid_indices.append(i)
                else:
                    complex_count += 1
            except Exception:
                # Fallback to sympy solve
                solutions = solve(expr_substituted, deriv_symbol)
                if solutions:
                    sol = complex(solutions[0].evalf())
                    if np.isreal(sol):
                        deriv_values[i] = sol.real
                        valid_indices.append(i)
                    else:
                        complex_count += 1
                else:
                    fail_count += 1
        except Exception as e:
            fail_count += 1
            continue
    
    print(f"Algebraic derivative solver: {len(valid_indices)} valid, {complex_count} complex, {fail_count} failed")
    
    if len(valid_indices) == 0:
        return np.zeros((0, 0))
    
    # -------------------------
    # 5) Build dataset with valid rows only
    # -------------------------
    valid_indices = np.array(valid_indices)

    columns = []
    col_names = []

    # Constants first
    for const in observed_constants:
        columns.append(np.full(len(valid_indices), constant_values[const]))
        col_names.append(const)

    # Derivative SECOND
    columns.append(deriv_values[valid_indices])
    col_names.append(deriv_name)

    # Measured variables LAST
    for var in measured_variables:
        columns.append(var_data[var][valid_indices])
        col_names.append(var)

    dataset = np.column_stack(columns)
    
    print(f"Generated algebraic dataset with shape {dataset.shape}")
    print(f"Column order: {col_names}")
    
    return dataset


def generate_dataset_with_ode(
    target_polynomial,
    measured_variables,
    observed_constants,
    measured_derivatives,
    constant_data=True,
    region=[1, 5],
    N=1000
):
    """
    Generate dataset using ODE solver when exactly one derivative is present.
    
    The workflow:
    1. Sample constants (fixed values)
    2. Sample all measured variables EXCEPT the derivative and its base (d1 or d2)
    3. Express the polynomial as an ODE by solving for the derivative
    4. Use solve_ivp to integrate and get consistent (derivative, base_var) pairs
    
    Args:
        target_polynomial: String representation of the consequence polynomial
        measured_variables: List of measured variable names
        observed_constants: List of constant names
        measured_derivatives: List containing exactly one derivative name
        constant_data: If True, constants are fixed; if False, sampled per row
        region: [min, max] range for sampling
        N: Number of data points to generate
    
    Returns:
        NumPy array with columns in order: constants, derivative, base_var, other_vars
        Returns empty array if generation fails
    """
    if len(measured_derivatives) != 1 or measured_derivatives[0] == '':
        print("ODE method requires exactly one derivative")
        return np.zeros((0, 0))
    
    deriv_name = measured_derivatives[0]
    if deriv_name not in DERIVATIVE_INFO:
        print(f"Unknown derivative: {deriv_name}")
        return np.zeros((0, 0))
    
    base_var, deriv_order = DERIVATIVE_INFO[deriv_name]
    
    # Check if base variable is in measured variables
    if base_var not in measured_variables:
        print(f"Base variable {base_var} not in measured variables, falling back to random sampling")
        return np.zeros((0, 0))  # Signal to use fallback method
    
    # -------------------------
    # 0) Setup constants
    # -------------------------
    constant_values = {}
    if observed_constants:
        if constant_data:
            for const in observed_constants:
                if const == 'c':
                    constant_values[const] = 1
                elif const == 'G':
                    constant_values[const] = 1
                else:
                    constant_values[const] = np.random.uniform(1, 10)
        else:
            for const in observed_constants:
                constant_values[const] = np.random.uniform(region[0], region[1])
    
    # -------------------------
    # 1) Handle theta and trig functions
    # -------------------------
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    
    # -------------------------
    # 2) Identify other measured variables (not derivative, not base_var)
    # -------------------------
    other_vars = [v for v in measured_variables if v != base_var and v not in specials]
    
    # -------------------------
    # 3) Parse polynomial and solve for the derivative symbolically
    # -------------------------
    poly_str = target_polynomial.replace('^', '**')
    try:
        poly_expr = sympify(poly_str)
    except Exception as e:
        print(f"Failed to parse polynomial: {e}")
        return np.zeros((0, 0))
    
    deriv_sym = Symbol(deriv_name)
    base_sym = Symbol(base_var)
    
    # Try to solve for the derivative
    try:
        solutions = solve(poly_expr, deriv_sym)
        if not solutions:
            print(f"Could not solve polynomial for {deriv_name}")
            return np.zeros((0, 0))
        # Take the first solution (usually the positive root for physical quantities)
        deriv_expr = solutions[0]
    except Exception as e:
        print(f"Failed to solve for derivative: {e}")
        return np.zeros((0, 0))
    
    print(f"Solved: {deriv_name} = {deriv_expr}")
    
    # -------------------------
    # 4) Generate data using ODE solver
    # -------------------------
    results = []
    failed_count = 0
    
    # Need to generate more attempts than N because some may fail
    max_attempts = N * 5
    
    for attempt in range(max_attempts):
        if len(results) >= N:
            break
        
        # Sample other variables for this data point
        other_var_values = {}
        theta_value = None
        
        # Handle theta specially - sample once and compute trig values
        if 'theta' in measured_variables or any(v in measured_variables for v in trig_funcs):
            theta_value = np.random.uniform(region[0], region[1])
        
        for v in other_vars:
            if v in trig_funcs and theta_value is not None:
                other_var_values[v] = trig_funcs[v](theta_value)
            elif v == 'theta':
                other_var_values[v] = theta_value
            else:
                other_var_values[v] = np.random.uniform(region[0], region[1])
        
        # Create substitution dict (everything except deriv and base_var)
        subs_dict = {}
        for const, val in constant_values.items():
            subs_dict[Symbol(const)] = val
        for var, val in other_var_values.items():
            subs_dict[Symbol(var)] = val
        
        # Add theta family to substitution dict (they were excluded from other_vars but need to be substituted)
        if theta_value is not None:
            if 'theta' in measured_variables:
                subs_dict[Symbol('theta')] = theta_value
            for trig_name, trig_func in trig_funcs.items():
                if trig_name in measured_variables:
                    subs_dict[Symbol(trig_name)] = trig_func(theta_value)
        
        # Substitute into derivative expression
        try:
            deriv_expr_numeric = deriv_expr.subs(subs_dict)
        except Exception as e:
            failed_count += 1
            continue
        
        # Now deriv_expr_numeric should be a function of base_var only
        # Create a lambda for the ODE RHS
        try:
            if deriv_order == 1:
                # First order: dx/dt = f(x)
                # The derivative expression gives dx/dt as function of x
                deriv_func = deriv_expr_numeric
                
                def ode_rhs(t, y):
                    x_val = y[0]
                    try:
                        result = float(deriv_func.subs(base_sym, x_val).evalf())
                        if np.isnan(result) or np.isinf(result):
                            return [0.0]  # Will be caught by success check
                        return [result]
                    except:
                        return [0.0]
                
                # Random initial condition
                x0 = np.random.uniform(region[0], region[1])
                y0 = [x0]
                
                # Solve the ODE
                t_span = (0, 1)
                t_eval = np.linspace(0, 1, 50)
                
                sol = solve_ivp(ode_rhs, t_span, y0, t_eval=t_eval, method='RK45')
                
                if sol.success and len(sol.t) > 10:
                    # Pick a random point from the trajectory (avoid very start)
                    idx = np.random.randint(5, len(sol.t))
                    base_val = sol.y[0, idx]
                    
                    # Compute derivative at this point
                    deriv_val = float(deriv_func.subs(base_sym, base_val).evalf())
                    
                    # Check validity
                    if np.isfinite(base_val) and np.isfinite(deriv_val):
                        if region[0] <= base_val <= region[1] * 10:  # Allow some range expansion
                            results.append({
                                'constants': constant_values,
                                'derivative': (deriv_name, deriv_val),
                                'base_var': (base_var, base_val),
                                'other_vars': other_var_values,
                                'theta': theta_value
                            })
                else:
                    failed_count += 1
                    
            else:  # deriv_order == 2
                # Second order: dÂ²x/dtÂ² = f(x) or dÂ²x/dtÂ² = f(x, dx/dt)
                # Convert to system: dx/dt = v, dv/dt = f(x, v)
                deriv_func = deriv_expr_numeric
                
                # Check if dx1dt or dx2dt appears in the expression
                first_deriv_name = 'dx1dt' if base_var == 'd1' else 'dx2dt'
                first_deriv_sym = Symbol(first_deriv_name)
                has_first_deriv = first_deriv_sym in deriv_func.free_symbols
                
                def ode_rhs(t, y):
                    x_val, v_val = y[0], y[1]
                    try:
                        subs = {base_sym: x_val}
                        if has_first_deriv:
                            subs[first_deriv_sym] = v_val
                        accel = float(deriv_func.subs(subs).evalf())
                        if np.isnan(accel) or np.isinf(accel):
                            return [0.0, 0.0]
                        return [v_val, accel]
                    except:
                        return [0.0, 0.0]
                
                # Random initial conditions
                x0 = np.random.uniform(region[0], region[1])
                v0 = np.random.uniform(-1, 1)  # Smaller range for velocity
                y0 = [x0, v0]
                
                t_span = (0, 1)
                t_eval = np.linspace(0, 1, 50)
                
                sol = solve_ivp(ode_rhs, t_span, y0, t_eval=t_eval, method='RK45')
                
                if sol.success and len(sol.t) > 10:
                    idx = np.random.randint(5, len(sol.t))
                    base_val = sol.y[0, idx]
                    
                    # Compute second derivative at this point
                    subs = {base_sym: base_val}
                    if has_first_deriv:
                        subs[first_deriv_sym] = sol.y[1, idx]
                    deriv_val = float(deriv_func.subs(subs).evalf())
                    
                    if np.isfinite(base_val) and np.isfinite(deriv_val):
                        if region[0] <= base_val <= region[1] * 10:
                            results.append({
                                'constants': constant_values,
                                'derivative': (deriv_name, deriv_val),
                                'base_var': (base_var, base_val),
                                'other_vars': other_var_values,
                                'theta': theta_value
                            })
                else:
                    failed_count += 1
                    
        except Exception as e:
            failed_count += 1
            continue
    
    print(f"ODE generation: {len(results)} successful, {failed_count} failed")
    
    if len(results) == 0:
        return np.zeros((0, 0))
    
    # -------------------------
    # 5) Convert results to numpy array with correct column order
    # -------------------------
    # Column order: constants + derivative + base_var + other measured_vars (preserving order)
    
    # Build the dataset
    n_rows = len(results)
    columns = []
    col_names = []
    
    # Constants first
    for const in observed_constants:
        columns.append(np.full(n_rows, results[0]['constants'][const]))
        col_names.append(const)
    
    # Derivative column
    columns.append(np.array([r['derivative'][1] for r in results]))
    col_names.append(deriv_name)
    
    # Base variable column
    columns.append(np.array([r['base_var'][1] for r in results]))
    col_names.append(base_var)
    
    # Other measured variables in their original order (excluding base_var which is already added)
    for v in measured_variables:
        if v == base_var:
            continue  # Already added
        if v == 'theta':
            columns.append(np.array([r['theta'] if r['theta'] is not None else 0 for r in results]))
            col_names.append(v)
        elif v in trig_funcs:  # ADD THIS CONDITION FOR TRIG FUNCTIONS
            # Compute trig function from theta
            columns.append(np.array([trig_funcs[v](r['theta']) if r['theta'] is not None else 0 for r in results]))
            col_names.append(v)
        elif v in results[0]['other_vars']:
            columns.append(np.array([r['other_vars'][v] for r in results]))
            col_names.append(v)
    
    dataset = np.column_stack(columns)
    
    print(f"Generated ODE dataset with shape {dataset.shape}")
    print(f"Column order: {col_names}")
    
    return dataset

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

def generate_dataset(
    target_polynomial,
    measured_variables,
    observed_constants,
    measured_derivatives,
    constant_data=True,
    derivative_data=True,
    region=[1, 5]
):
    N = 1000
    dataset = np.zeros((N, 0))  # start empty

    # --------------------------
    # 0) Constants (same logic)
    # --------------------------
    const_data_cols = []
    if observed_constants:
        if constant_data:
            # one draw per constant, repeat N times
            constant_values = {const: np.random.uniform(0, 10) for const in observed_constants}
            if 'c' in observed_constants:
                constant_values['c'] = 1
            if 'G' in observed_constants:
                constant_values['G'] = 1
            for const in observed_constants:
                const_data_cols.append(np.full(N, constant_values[const]))
        else:
            # N samples per constant
            for const in observed_constants:
                const_data_cols.append(np.random.uniform(region[0], region[1], N))
    if const_data_cols:
        dataset = np.hstack([dataset] + [col.reshape(-1, 1) for col in const_data_cols])

    # ------------------------------------
    # 1) Derivatives (now strictly ordered)
    # ------------------------------------
    derivative_dependent_variable = []  # will contain 'd1' and/or 'd2' that were used to build derivatives
    derivative_arrays = {}              # map derivative name -> np.array

    # Prepare d1 branch (if needed)
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

    # Prepare d2 branch (if needed)
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

            derivative_dependent_variable.append('d2')
            measured_variables = [v for v in measured_variables if v != 'd2']
        else:
            if 'dx2dt' in measured_derivatives:
                dx2dt_data = np.random.uniform(region[0], region[1], N)
                derivative_arrays['dx2dt'] = dx2dt_data
                if 'd2x2dt2' in measured_derivatives:
                    derivative_arrays['d2x2dt2'] = np.diff(dx2dt_data, append=dx2dt_data[0])
            elif 'd2x2dt2' in measured_derivatives:
                derivative_arrays['d2x2dt2'] = np.random.uniform(region[0], region[1], N)

    # Append derivative columns strictly in measured_derivatives' order (ignore sentinel '')
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

    # ---------------------------------------------------------
    # 2) Theta & trig/exponential dependents (computed values)
    # ---------------------------------------------------------
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    trig_vars_present = [v for v in measured_variables if v in trig_funcs]
    need_theta_stream = ('theta' in measured_variables) or bool(trig_vars_present)

    theta_stream = None
    if need_theta_stream:
        # sample theta once for the whole dataset
        theta_stream = np.random.uniform(region[0], region[1], N)

    trig_values = {}
    for v in trig_vars_present:
        trig_values[v] = trig_funcs[v](theta_stream)

    # ---------------------------------------------------------
    # 3) Choose the last variable to solve for (exclude specials)
    # ---------------------------------------------------------
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    # Preserve measured order; pick last that isn't a special
    solve_candidates = [v for v in measured_variables if v not in specials]
    if len(solve_candidates) == 0:
        # Nothing left to solve for; signal failure
        return np.zeros((0, 0))

    last_var = solve_candidates[-1]
    known_measured_vars = [v for v in measured_variables if v != last_var]

    # ---------------------------------------------------------
    # 4) Generate measured vars (except last_var) in correct order
    # ---------------------------------------------------------
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

    # ---------------------------------------------------------
    # 5) Solve for last_var row-by-row
    # ---------------------------------------------------------
    last_var_symbol = symbols(last_var)
    roots_column = np.full((N, 1), np.nan)

    complex_count = 0
    evaluation_fail_count = 0

    # Build the substitution ordering that matches dataset columns:
    # constants + measured_derivatives + derivative_dependent_variable + known_measured_vars
    sub_order = []
    sub_order.extend(observed_constants)
    if measured_derivatives and '' not in measured_derivatives:
        sub_order.extend(measured_derivatives)
        sub_order.extend(derivative_dependent_variable)
    sub_order.extend(known_measured_vars)

    for i in range(N):
        polynomial = target_polynomial

        # Replace vars with placeholders so 'sympify' sees scalars
        for name in (observed_constants +
                     (measured_derivatives if (measured_derivatives and '' not in measured_derivatives) else []) +
                     derivative_dependent_variable +
                     known_measured_vars):
            polynomial = re.sub(rf'\b{re.escape(name)}\b', f'{name}_value', polynomial)

        polynomial = polynomial.replace('^', '**')

        try:
            expr = sympify(polynomial)
            substitutions = {}
            # Map the sub_order to the dataset columns (same order we appended)
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

def generate_dataset_no_der(
    target_polynomial,
    measured_variables,
    observed_constants,
    constant_data=True,
    region=[1, 5]
):
    N = 1000
    dataset = np.zeros((N, 0))

    # 0) Constants
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

    # 1) Theta + trig dependents
    trig_funcs = {'sinTheta': np.sin, 'cosTheta': np.cos, 'eTheta': np.exp}
    trig_vars_present = [v for v in measured_variables if v in trig_funcs]
    need_theta_stream = ('theta' in measured_variables) or bool(trig_vars_present)

    theta_stream = None
    if need_theta_stream:
        theta_stream = np.random.uniform(region[0], region[1], N)

    trig_values = {}
    for v in trig_vars_present:
        trig_values[v] = trig_funcs[v](theta_stream)

    # 2) Pick last var to solve (exclude specials)
    specials = {'theta', 'sinTheta', 'cosTheta', 'eTheta'}
    solve_candidates = [v for v in measured_variables if v not in specials]
    if len(solve_candidates) == 0:
        return np.zeros((0, 0))

    last_var = solve_candidates[-1]
    known_measured_vars = [v for v in measured_variables if v != last_var]

    # 3) Generate measured (except last) in correct order
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

    # 4) Solve for last_var
    last_var_symbol = symbols(last_var)
    roots_column = np.full((N, 1), np.nan)

    complex_count = 0
    evaluation_fail_count = 0

    # Substitution order here: constants + known_measured_vars
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
    
    # Handle empty derivatives case
    if not measured_derivatives:
        measured_derivatives = ['']
    
    # Check if we should use ODE-based method (exactly one derivative)
    use_ode_method = (
        len(measured_derivatives) == 1 and 
        measured_derivatives[0] != '' and 
        measured_derivatives[0] in DERIVATIVE_INFO
    )
    
    if use_ode_method:
        print(f"Using ODE-based data generation for derivative: {measured_derivatives[0]}")
        
        # For ODE method, we don't need to reorder - the function handles column ordering
        dataset = generate_dataset_with_ode(
            target_polynomial, 
            measured_variables.copy(), 
            observed_constants, 
            measured_derivatives,
            constant_data=True,
            region=region,
            N=1000
        )
        
        if len(dataset) == 0:
            print("ODE method failed, falling back to finite difference method")
            use_ode_method = False
        else:
            # Update the file with proper variable ordering for ODE output
            # ODE output column order: constants + derivative + base_var + other_vars
            deriv_name = measured_derivatives[0]
            base_var, _ = DERIVATIVE_INFO[deriv_name]
            
            # Reorder measured_variables to match ODE output: base_var first, then others
            new_measured_vars = [base_var] if base_var in measured_variables else []
            for var in measured_variables:
                if var != base_var:
                    new_measured_vars.append(var)
            
            if new_measured_vars != measured_variables:
                updated_content = original_content.replace(
                    f"Measured Variables: {measured_variables}",
                    f"Measured Variables: {new_measured_vars}"
                )
                with open(input_file, 'w') as f:
                    f.write(updated_content)
                measured_variables = new_measured_vars.copy()
    
    if not use_ode_method:
        # Try algebraic derivative solving first (cleaner than finite differences)
        if len(measured_derivatives) == 1 and measured_derivatives[0] != '' and measured_derivatives[0] in DERIVATIVE_INFO:
            print(f"Trying algebraic derivative solving for: {measured_derivatives[0]}")
            dataset = generate_dataset_algebraic_derivative(
                target_polynomial,
                measured_variables.copy(),
                observed_constants,
                measured_derivatives,
                constant_data=True,
                region=region,
                N=1000
            )
            
            if len(dataset) > 0:
                # Algebraic method succeeded
                print("Algebraic derivative method succeeded")
            else:
                # Both ODE and algebraic failed - give up on this system
                print("Algebraic method failed. Cannot generate data for this consequence.")
                return False
                
        elif not measured_derivatives or measured_derivatives[0] == '':
            # No derivatives case - solve normally
            print("No derivatives present, using standard algebraic solving.")
            dataset = generate_dataset_no_der(
                target_polynomial, 
                measured_variables.copy(), 
                observed_constants, 
                constant_data=True, 
                region=region
            )
        else:
            # Unknown derivative type or multiple derivatives - cannot handle
            print(f"Cannot handle derivative configuration: {measured_derivatives}")
            return False

    if len(dataset) == 0:
        print("Dataset Generation Failed. Skipping \n")
        return False
    
    # Save dataset to file
    np.savetxt(output_file, dataset, delimiter=' ')
    print(f"Generated dataset for consequence with shape {dataset.shape}")
    print("Evaluation: ", sum(evaluate_polynomial(target_polynomial, dataset[0:5], observed_constants, measured_derivatives, measured_variables)))

    return True

def run_consequence_noisy_data_generation(input_file):
    epsilon_list = [1e-3,1e-2,1e-1,5e-2]

    for epsilon in epsilon_list:
        output_file = input_file.replace('.dat', f'_{epsilon}.dat')
        add_gaussian_noise(input_file, output_file,epsilon)
