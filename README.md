# Lennard-Jones Gas Simulation
This project aim is to reproduce the behavior of a Lennard-Jones gas by using the Monte Carlo method, 
and to study different phenomena, as phase changes, entropy and condensations.

## Compilation 
All code must be compiled using c++11. A standard compilation would be 
```bash
g++ -std=c++11 program_name_here.cpp
```
However, makefiles are included to save time and effort in this matter for some files, 
and diferent debugging, profiling and optimization checks are also automatized. Please do 
check the following table:

### Implemented Automatization
| Shell Command   | What it does oversimplified                                               |
|-----------------|---------------------------------------------------------------------------|
| make            | Compiles using O2 optimization, c++11 standard and runs program using time|
| make animation  | Same as make, redirects program output to gnuplot in real time            |
| make save_gif   | Same as animation, but saves animation in a gif                           |
| make debug      | Compiles with all warning flags, saves debug info using -g flag and runs  |
| make valgrind   | Compiles with -g flag and runs using valgrind                             |
| make cachegrind | Same as valgrind but using runs using cachegrind                          |
| make gprof      | Compiles with -pg flag, runs and display gprof info                       |
| make perf       | Same as gprof, but runs and displays perf info                            |