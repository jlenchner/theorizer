# SynPAT: Generating Synthetic Physical Theories with Data

This is a GitHub repository for the synthetic theory and data generation system developed for the paper **SynPAT: Generating Synthetic Physical Theories with Data**. It generates symbolic physical systems and corresponding synthetic data to benchmark symbolic regression and scientific discovery algorithms.

---

## Prerequisites

### System Requirements
- Python 3.8+
- NumPy
- SymPy
- **Macaulay2** (required for consequence generation via algebraic projection)

### Installing Dependencies

#### Python Dependencies
Install with:
```bash
pip install numpy sympy
```

#### Macaulay2 Installation

- **macOS (Homebrew):**
  ```bash
  brew install macaulay2
  ```
- **Linux (Debian/Ubuntu):** Download and install from [Macaulay2's official website](https://macaulay2.com/Downloads/)
- **Windows:** Use the Windows installer from [Macaulay2's official website](https://macaulay2.com/Downloads/)

To verify the installation, run:
```bash
M2
```

---

## Usage: Generating Synthetic Systems and Data

The main script is `run_generation.py`. It generates symbolic systems, derivative structures, replacements, and consequences, and produces datasets with optional noise.

### Example Run

```bash
python run_generation.py \
    --vars "['m1', 'm2', 'd1', 'd2', 'W', 'p', 'Fc', 'Fg']" \
    --derivs "['dx1dt', 'dx2dt', 'd2x1dt2', 'd2x2dt2']" \
    --numVars "[6,7,8]" \
    --numDerivs "[2,3]" \
    --numEquations "[4,5]" \
    --numReplacements 5 \
    --numSystems 2 \
    --genReplacements \
    --genSysData \
    --genConsequence \
    --genConsequenceData \
    --numConstConseq 1 \
    --conseqDataRange "[1, 5]" \
    --sysDataRange "[2, 6]" \
    --timeout 3600 \
    --seed 123
```

### Argument Descriptions

| Argument               | Description                                                                 | Default |
|------------------------|-----------------------------------------------------------------------------|---------|
| `--vars`               | List of symbolic base variables used to construct the system.              | `['Fc', 'Fg', 'W', 'd1', 'd2', 'm1', 'm2', 'p', 'T', 'E']` |
| `--derivs`             | List of symbolic derivative names used as inputs.                          | `['dx1dt', 'd2x1dt2', 'dx2dt', 'd2x2dt2']` |
| `--numVars`            | List of integers specifying how many base variables to include.            | `[6, 7, 8, 9]` |
| `--numDerivs`          | List of integers specifying how many derivatives to include.               | `[2, 3, 4]` |
| `--numEquations`       | List of integers specifying how many equations should appear in each system. | `[4, 5, 6]` |
| `--numReplacements`    | Number of replacement axioms to generate per system.                       | `5` |
| `--numSystems`         | Number of independent systems to generate.                                 | `3` |
| `--numConstConseq`     | Maximum number of constants allowed in the consequence.                    | `1` |
| `--genReplacements`    | Enable generation of dimensionally consistent replacement axioms.          | `True` |
| `--genSysData`         | Generate synthetic data for each axiom system.                             | `True` |
| `--genConsequence`     | Compute an algebraic consequence using Macaulay2.                          | `True` |
| `--genConsequenceData` | Generate data for the derived consequence.                                 | `True` |
| `--conseqDataRange`    | Range `[start, end]` for consequence data sampling.                        | `[1, 10]` |
| `--sysDataRange`       | Range `[start, end]` for system data sampling.                             | `[1, 10]` |
| `--timeout`            | Max number of seconds to allow for generating a single system. `None` disables the timeout. | `3600` |
| `--seed`               | Random seed for reproducibility.                                           | `42` |

---

## Output Directory Structure

Generated systems and data are saved under:

```
dataset/{config}/System_{n}/
```

Example contents:

```
system.txt              # Symbolic equation system
system.dat              # Clean solution data
system_1e-1.dat         # Noisy data (ε = 1e-1)
system_1e-4.dat         # Noisy data (ε = 1e-4)
consequence.txt         # Derived symbolic consequence
consequence.dat         # Clean solution for consequence
consequence_1e-1.dat    # Noisy consequence data
replacement_i.txt       # ith system with replaced axiom 
```

Performance metrics (runtime and memory usage) are saved to:

```
Data_Gen_Statistics/{config}/System_{n}/performance.txt
```

---

## Dataset Download

The generated datasets are now hosted and available for download at [Huggingface SynPAT Dataset](https://huggingface.co/datasets/Karan0901/synpat-dataset).

You can clone or download the dataset directly from there as well. 

---

## Notes

- All ranges such as `--conseqDataRange` and `--sysDataRange` must be valid `[start, end]` pairs with `start < end`. Invalid formats will raise an error.
- The script automatically handles cleanup and only saves files relevant to the requested generation flags.
- All consequence generation using the elimination theorem is performed through Macaulay2 via the `run_consequence_generation()` function in `generate_dataset.py`.
- Noisy data variants are created for both system and consequence equations to support evaluation under data perturbations.

---

## License

This code is licensed under the MIT License.

---

For questions, issues, or contributions, please open an issue or pull request!
