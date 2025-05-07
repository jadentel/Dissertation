
# Protein Folding with Genetic Algorithm (2D HP Model)

This project simulates protein folding using a Genetic Algorithm (GA) based on the 2D Hydrophobic-Polar (HP) lattice model. It supports various mutation strategies, tournament-based crossover, ASCII visualization, performance logging, and Bayesian Optimization for hyperparameter tuning.
![name image if not load](image name/dir)
---

## Project Structure

```
.
├── main.py                   # Entry point: runs GA on a test protein
├── Conformation.py           # Core folding logic, mutation, validation, fitness
├── Population.py             # Handles population initialization and evolution
├── Protein.py                # Protein sequence abstraction
├── bays.py                   # Bayesian Optimization of GA parameters
├── testing.py                # Performance benchmarking and CSV logging
├── protein_sequence.txt      # Example protein sequences with known optima
│
├── test/
│   ├── test.py               # Unit tests for all components
│   └── visual_utils.py       # Colorful visualization tools for mutations
```

---

## How to Run:

### Standard GA Run

```bash
python main.py
```

This runs the GA on the 24-length sequence `BBWWBWWBWWBWWBWWBWWBWWBB`. It prints:
- Final fitness
- ASCII diagram of the fold
- Conformation encoding string

---

## Unit Testing

Navigate to the `test` folder and run:

```bash
pytest unittest test.py
pytest unittest test.py -s # to print outputs
```

This tests:
- Protein sequence behavior
- Conformation validity, fitness, mutation correctness
- Population selection, crossover, insertion logic
- Visual difference between mutated conformations

---

## Performance Benchmarking

Run multiple silent GA experiments and analyze results:

```bash
python testing.py
```

Creates and saves a CSV file `ga_24seq_results.csv` with:
- Best fitness per run
- Evaluations to best
- Unique conformations found
- Generations and birth generation of best fold

Also prints a statistical summary (mean, stdev, min, max) to terminal.

---

## Bayesian Optimization

Use Bayesian Optimization to find the best GA parameters:

```bash
python bays.py
```

Tunes:
- `population_size` (500–2000)
- `mutation_probability` (0.01–0.4)
- `crossover_probability` (0.4–0.9)

Reports the best parameters and resulting GA performance.

---


## Protein Sequence Format

Proteins are sequences of:
- `'B'` – Hydrophobic (black)
- `'W'` – Polar (white)

Example from `protein_sequence.txt`:

```
BBWWBWWBWWBWWBWWBWWBWWBB  # optimal energy = -9
```

---

## Requirements

Install dependencies:

```bash
pip install bayesian-optimization termcolor
```

---

## Usage

![85](https://github.com/user-attachments/assets/6b61ae07-61f8-4bfd-94e1-fb5537ef007e)




