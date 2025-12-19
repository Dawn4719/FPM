## Fragmented Graph Pattern Matching on Large Graphs
This repository provides the implementation of algorithms for Fragmented Graph Pattern Matching problem.
It includes our proposed method (`FPM`) as well as baselines (`BsLn`, `BFraMatch`, `DFraMatch`, `SFraMatch`).

## Code Structure
The main implementations are located in the `matching` directory:

- `/matching/FPM.cpp` – Implementation of `FPM`, our proposed algorithm for solving the Fragmented Graph Pattern Matching problem.
- `/matching/FPMbase.cpp` – Implementations of the baseline methods used for comparison with `FPM`.
- `Dataset` – Example dataset.

---

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `FPM` in `build/matching` directory.

---

## Experiment Datasets
The datasets used in our paper can be downloaded [CC](https://networkrepository.com/MSRC-21C.php), [MC](https://networkrepository.com/MSRC-21.php), [DD](https://networkrepository.com/DD.php), [CL](https://networkrepository.com/CL-10M-1d8-L5.php).

---

## Run the procedure
### Example Commands
- **Run query on `CC` dataset (`-d CC`) with 3 queries (`-Q 3`) and distance threshold 1 (`-D 1`) on queryset 1 (`-q 1`) in `FPM` algorithm:**:
```zsh
./matching/FPM.out -q 1 -d CC -Q 3 -D 1
```
- **Run query on `CC` dataset (`-d CC`) with 3 queries (`-Q 3`) and distance threshold 1 (`-D 1`) on queryset 1 (`-q 1`) in `FPM` algorithm with 4 threads:**
```zsh
./matching/FPM.out -q 1 -d CC -Q 3 -D 1 -T 4
```
- **Run query on `CC` dataset (`-d CC`) with 3 queries (`-Q 3`) and distance threshold 1 (`-D 1`) on queryset 1 (`-q 1`) in `BsLn` algorithm:**:
```zsh
./matching/FPMBase.out -q 1 -d CC -Q 3 -D 1 -M 0
```

### Command Line Options

| Option | Required | Description                                                                                          |
|--------|----------|------------------------------------------------------------------------------------------------------|
| `-q`   | Yes      | Queryset name                                                                                        |
| `-d`   | Yes      | Dataset name                                                                                         |
| `-Q`   | Yes      | Query number in queryset                                                                             |
| `-D`   | Yes      | Distance threshold                                                                                   |
| `-M`   | Yes      | Method selection (`0`: BsLn, `1`: BFraMatch, `2`: DFraMatch, `3`: SFraMatch, **only for `FPMBase`**) |
| `-T`   | No       | Number of threads (**only for `FPM`**)                                                               |

---

This codebase is implemented on the framework provided by [SubgraphMatching](https://github.com/RapidsAtHKUST/SubgraphMatching)