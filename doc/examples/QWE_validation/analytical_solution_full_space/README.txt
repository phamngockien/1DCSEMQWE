This folder contains a C++ code that generates the analytical solution based on Ward and Hohmann (1988).
It can give the EM field of either x-oriented ED or z-oriented MD in a homogeneous fullspace.

You may want to change the code in analytical_solution.cc for a different model simulation.

To compile the code, you can use the command:
g++ analytical_solution.cc HomoBackground.cpp

Then you can compute the EM field by the command:
./a.out
