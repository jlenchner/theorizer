# SynPAT: Generating Synthetic Physical Theories with Data

This is a github repository for the synthetic theory and data generation system for the paper: (Insert paper link when made public)

## Prerequisites

### System Requirements
- Python 3.8+
- NumPy
- SymPy
- Macaulay2 (for projection analysis)

### Installing Dependencies

#### Python Dependencies
```bash
pip install numpy sympy
```

#### Macaulay2 Installation

   - macOS (using Homebrew): `brew install macaulay2`
   - Linux (Debian/Ubuntu): [Macaulay2's official website](https://macaulay2.com/Downloads/)
   - Windows: Download from [Macaulay2's official website](https://macaulay2.com/Downloads/)
   - Verify installation by running `M2` in your terminal

## Consequence Generation

### Projection Function (`projection`)
Performs algebraic projections using Macaulay2 to generate a consequence for a given axiom on a set of measured variables as described in the paper.

#### Usage Example
```python
result = projection(
    variables=['x', 'y', 'z'],
    axioms=['x^2 + y^2 - z^2'],
    measured_variables=['x', 'y'],
    non_measured_variables=['z'],
    filename='projection_results.txt'
)
```

### 2. Dataset Generation Functions

#### `generate_dataset`
Generates a synthetic dataset for a given consequence for the variables sampled from a uniform distribution of a specified region. 

#### Usage Example
```python
dataset = generate_dataset(
    target_polynomial='x^2 + a*x + b',
    axioms=['x^2 + y^2 - z^2'],
    observed_constants=['a', 'b'],
    observed_derivatives=['dx1dt'],
    region=[1, 5]
)
```

#### 'generate_data_for_system' 
Generates a synthetic dataset for a given system of axioms using the Groebner basis for the variables sampled from a uniform distribution of a specified region as described in the paper. 

#### Usage Example
```python
dataset = generate_dataset_for_system(
    measured_variables=['x', 'a', 'b'],
    observed_constants=['a', 'b'],
    observed_derivatives=['dx1dt'],
    region=[1, 5]
)
```

