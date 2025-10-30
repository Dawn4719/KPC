## Enumerating Cliques on $k$-partite Graphs
This repository provides the implementation of algorithms for **Enumerating Cliques on $k$-partite Graphs**. The main implementation is in KPC.cpp

---
## Compile the codes
When you have already downloaded the codes, run the following commands to compile our codes.
```
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `KPC` in `build` directory.

---
## Datasets
- **Default datasets** (e.g., `congress4P`, `congress5P`, `congress6P`) are included in the package.
- **Additional datasets** can be downloaded from [networkrepository](https://networkrepository.com/).

---
## Run the procedure
### Example Commands
- **Run $\lambda$-k-partite Clique query(5-4-partite Clique):**
```bash
./KPC 4 5 ./Dataset/congress4P 1
```
- **Run $(\alpha,\beta)$-k-partite Clique((1, 2)-4-partite Clique):**
```bash
./KPC 4 -1 ./Dataset/congress4P 2 1 2
```
- **Run $\delta$-balanced k-partite Clique(1-balanced 4-partite Clique):**
```bash
./KPC 4 1 ./Dataset/congress4P 3
```

### Command Line Options

| Argvs | Description |
|--------|-------------|
| `argv[1]`   | The number of partite sets |
| `argv[2]`   | The value of $\lambda$ (WHEN `argv[4]==1`) or $\delta$ (WHEN `argv[4]==3`)|
| `argv[3]`   | The path of dataset |
| `argv[4]`   | The type of $k$-partite clique (1 for $\lambda$-k-partite Clique$, 2 for $(\alpha,\beta)$-k-partite Clique$, 3 for $\delta$-balanced k-partite Clique$) |
| `argv[6]`   | The value of $\alpha$ |
| `argv[7]`   | The value of $\beta$ |

---

## External Dependencies
- MEBA: available at [C++ veresion](https://github.com/ddervs/MBEA) and [JAVA version](https://github.com/gaurabdg/maximal-biclique-enumeration)
---
