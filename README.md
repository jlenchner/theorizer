# SynPAT: Generating Synthetic Physical Theories with Data

This is a GitHub repository for the synthetic theory and data generation system developed for the paper [**SynPAT: Generating Synthetic Physical Theories with Data**](https://www.arxiv.org/abs/2505.00878). It generates symbolic physical systems and corresponding synthetic data to benchmark symbolic regression and scientific discovery algorithms.

---

## Prerequisites

### System Requirements
- Python 3.9+
- NumPy
- SymPy
- Deprecated
- psutil
- **Macaulay2** (required for consequence generation via algebraic projection)

### Installing Dependencies

#### Python Dependencies
Install with:
```bash
pip install numpy sympy Deprecated psutil
```

Alternatively, all python dependencies can be installed using the included requirements file:

```bash
pip install -r requirements.txt
```

#### Macaulay2 Installation

- **macOS (Homebrew):**
  ```bash
  brew install Macaulay2/tap/M2
  ```
- **Linux (Debian/Ubuntu):** Download and install from [Macaulay2's official website](https://macaulay2.com/Downloads/)
- **Windows:** Use the Windows installer from [Macaulay2's official website](https://macaulay2.com/Downloads/)

To verify the installation, run:
```bash
M2
```

If installed correctly, you should see

```bash
Macaulay2, version <your version number>
Type "help" to see useful commands

i1 : 
```

To exit this, you can run

```bash
i1 : exit()
```

---

## Usage: Generating Synthetic Systems and Data

The main script is `run_generation.py`. It generates symbolic systems, derivative structures, replacements, and consequences, and produces datasets with optional noise.

### Example Run

```bash
python generate_dataset.py \
    --vars "['m1', 'm2', 'd1', 'd2', 'W', 'theta', 'sinTheta', 'p', 'Fc', 'Fg']" \
    --derivs "['dx1dt', 'dx2dt', 'd2x1dt2', 'd2x2dt2']" \
    --numVars "[6,7,8]" \
    --numDerivs "[2,3]" \
    --numEquations "[4,5]" \
    --numReplacements 5 \
    --numSystems 2 \
    --genReplacements True\
    --genSysData False\
    --genConsequence True\
    --genConsequenceData True\
    --numConstConseq 1 \
    --conseqDataRange "[1, 5]" \
    --sysDataRange "[2, 6]" \
    --maxConseqTerms "8" \
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
| `--numSystems`         | Number of independent systems to generate for each configuration.          | `3` |
| `--numConstConseq`     | Maximum number of constants allowed in the consequence.                    | `1` |
| `--genReplacements`    | Enable generation of dimensionally consistent replacement axioms.          | `True` |
| `--genSysData`         | Generate synthetic data for each axiom system.                             | `True` |
| `--genConsequence`     | Compute an algebraic consequence using Macaulay2.                          | `True` |
| `--genConsequenceData` | Generate data for the derived consequence.                                 | `True` |
| `--conseqDataRange`    | Range `[start, end]` for consequence data sampling.                        | `[1, 10]` |
| `--sysDataRange`       | Range `[start, end]` for system data sampling.                             | `[1, 10]` |
| `--maxConseqTerms`     | Upper bound on the number of monomials allowed in the consequence polynomial                             | `8` |
| `--timeout`            | Max number of seconds to allow for generating a single system. `None` disables the timeout. | `3600` |
| `--seed`               | Random seed for reproducibility.                                           | `42` |

---

## Output Directory Structure

Generated systems and data are saved under:

```
dataset/{config}/System_{n}/
```

Dataset contents:

```
- **Symbolic System**
  - `system.txt` — Base symbolic system with variable types, units, constants, and equations
  - `consequence.txt` — Derived consequence polynomial with metadata

- **Noiseless Data**
  - `system.dat` — Clean data for `system.txt`
  - `consequence.dat` — Clean data for `consequence.txt`

- **Noisy Data** (Gaussian noise with ε ∈ {0.001, 0.01, 0.05, 0.1}):
  - `system_ε.dat` — Noisy system data
  - `consequence_ε.dat` — Noisy consequence data

- **Replacement Systems**
  - `replacement_1.txt` to `replacement_5.txt` — `system.txt` with one axiom replaced (for ablation studies)

```

Performance metrics (runtime and memory usage) are saved to:

```
dataset_gen_statistics/{config}/System_{n}/performance.txt
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

# Third-Party Licenses

Our code is released under the MIT License (see `LICENSE`).

This project uses the following third-party software:

- **SymPy** — BSD 3-Clause  
  License: https://github.com/sympy/sympy/blob/master/LICENSE

- **Macaulay2** — GNU GPL (v2 or later / v3)  
  License: https://macaulay2.com/Downloads/Copyright/

---

## How to Cite

If you use **SynPAT** (code or datasets), please cite the paper and (optionally) the dataset.

**Paper**

> Jon Lenchner, Karan Srivastava, Joao Goncalves, Lior Horesh. *SynPAT: Generating Synthetic Physical Theories with Data*. arXiv:2505.00878.

- Preprint: https://www.arxiv.org/abs/2505.00878

```bibtex
@misc{synpat2025,
  title         = {SynPAT: Generating Synthetic Physical Theories with Data},
  author        = {Karan Srivastava},
  year          = {2025},
  eprint        = {2505.00878},
  archivePrefix = {arXiv},
  primaryClass  = {cs.LG},
  url           = {https://arxiv.org/abs/2505.00878}
}
