## Enumerating Cliques on $k$-partite Graphs
This repository provides the implementation of algorithms for **Enumerating Cliques on $k$-partite Graphs**. The main implementation is in KPC.cpp

---
## Compile the codes
When you already download the codes, run the following commands to compile our codes.
```
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `KPC` in `build` directory.

---
## Datasets
- **Default datasets** (e.g., `congress4KC`, `congress5KC`, `congress6KC`) are included in the package.
- **Additional datasets** can be downloaded from [networkrepository](https://networkrepository.com/).

---
## Run the procedure
### Example Commands  
- **Run $\lambda$-k-partite Clique query(5-4-partite Clique):**    
```bash
./KPC 4 5 ./Dataset/congress4KC -1 0
```
- **Run $(\alpha,\beta)$-k-partite Clique((1, 2)-4-partite Clique):**
```bash
./KPC 4 -1 ./Dataset/congress4KC 1 2
```
- **Run $\delta$-balanced k-partite Clique(1-balanced 4-partite Clique):**
```bash
./KPC 4 1 ./Dataset/congress4KC 0 -1
```

### Command Line Options

| Argvs | Description |
|--------|-------------|
| `argv[1]`   | The number of partite sets |
| `argv[2]`   | The value of $lambda$ (WHEN `argv[4]==-1 && argv[5]==0`) or $\delta$ (WHEN `argv[4]==0 && argv[5]==-1`)|
| `argv[3]`   | The path of dataset |
| `argv[4]`   | The value of $\alpha$ |
| `argv[5]`   | The value of $\beta$ |

---
