## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

```zsh
./matching/SubgraphMatching.out -q 1 -d CC -Q 3 -D 1 -M 0    
```

##
`StudyPerformance.cpp` in `./matching` is the code of FraMatch.
`StudyPerformance_old_sys_cal.cpp` in `./matching` is the code of BsLn, BFraMatch, DFraMatch.
`StudyPerformance_thread.cpp` in `./matching` is the muti-thread version of FraMatch.

## Experiment Datasets
The datasets used in our paper can be downloaded [CC](https://networkrepository.com/MSRC-21C.php), [MC](https://networkrepository.com/MSRC-21.php), [DD](https://networkrepository.com/DD.php), [CL](https://networkrepository.com/CL-10M-1d8-L5.php).
