## Investigation of Finite Element Method Approaches on Nvidia GPUs

Program solves the 2D Laplace equation:  

<center> &#8711;<sup>2</sup> &phi;(<b>x</b>) = 0, <b>x</b> &#8712; &Omega; </center> <br>
&phi;(<b>x</b>) = a, <b>x</b> &#8712; &Gamma;<sub>Da</sub>  <br>
&phi;(<b>x</b>) = b, <b>x</b> &#8712; &Gamma;<sub>Db</sub> <br>  
&part;&phi;(<b>x</b>) &frasl; &part; n = a, <b>x</b> &#8712; &Gamma;<sub>N</sub>, <br>

using the finite element method.

There are three separate implementations within the code:
* Serial code implementing standard FEM in C++.
* CUDA version implementing standard FEM with linear system decomposed across threads.
* CUDA version implementing novel FEM-SES (Single Element Solution) approach.
All code is stored in `/src/`

### Installing

### Usage

### Testing

### Author
