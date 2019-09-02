## Investigation of Finite Element Method Approaches on Nvidia GPUs

Program solves the 2D Laplace equation:  

<center> &#8711;<sup>2</sup> &phi;(<b>x</b>) = 0, <b>x</b> &#8712; &Omega; <br>
&phi;(<b>x</b>) = a, <b>x</b> &#8712; &Gamma;<sub>Da</sub>  <br>
&phi;(<b>x</b>) = b, <b>x</b> &#8712; &Gamma;<sub>Db</sub> <br>  
&part;&phi;(<b>x</b>) &frasl; &part;<b>n</b> = 0, <b>x</b> &#8712; &Gamma;<sub>N</sub>, <br>
</center>

using the finite element method.

There are three separate implementations within the code:
* Serial code implementing standard FEM in C++.
* CUDA version implementing standard FEM with linear system partitioned across threads.
* CUDA version implementing novel FEM-SES (Single Element Solution) approach.  

All code is stored in `/src/`

### Installing

#### Dependencies
The program `fem_solver` was built on `cuda01` with the following dependencies:
* C++11
* gcc 5.4.0
* CUDA 10.1
* Intel MKL 18.0.4  

To uphold these, make sure enviroment variables `PATH` and `LD_LIBRARY_PATH` are corrected, ammending them by:
```
export PATH=$PATH:/usr/local/cuda-10.1/bin:/home/support/apps/intel/18.0.4/bin/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/lib64:/home/support/apps/intel/18.0.4/mkl/lib/intel64/
```
#### Make
To make the program, simply run
```
cd src/
make
```

To remove all object files and executable, run,
```
make clean
```
or to remove all objects, timings, results and executable, run
```
make clean_all
```

### Usage
To get detailed usage information, run
```
./fem_solver -h
```
All results are stored in `/src/results/output_*_results.csv` and run times are stored in `/src/timings/*_timings.csv`.

Typical succesful run output looks as the screenshot below:

![Run Screen](run_screen.png?raw=true "Runtime Screenshot")

<b> Note: </b> When running this solver, there is no means of knowing a priori the amount of memory needed to be allocated by cuSOLVER library. Be conscious of DRAM on device when choose problem size.

An error log will be created and written to `/src/error.log` during runtime.

### Testing

To run test suite for timings, run
```
make timings
```
or to test for memory errors,
```
make mem_tests
```
The timings test script is lengthy - contains around 6,000 runs - so do not run if short on time. Output plots can be attained from the R code stored in `/src/plots/plots.R`

### Author
Alex Keating  
Completed as part of a research project for the attainment of an MSc. High Performance Computing.
