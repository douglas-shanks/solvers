###########################################################################
===========================================================================

1D and 2D sequential
--------------------

1D and 2D test programs to solve the heat equation
with Dirichlet boundary conditions and initial condition 

1D:

u_0(x,0) = sin(pi*x).

2D:

u_0(x,y,0) = sin(pi*x)sin(pi*y).

To run:

make clean
make

./test1D

or 

./test2D

This generates the exact solution, FTCS solution, BTCS solution
and CN solution. The 2 norm of the error is then computed for each and wall clock time computed.

2D Parallel
------------
We solve the heat equation in 2D using FTCS and BTCS (CN might be added) with 
Dirichlet boundary and chosen initial condition.
This version has a variable diffusion coefficient making the problem harder to
solve numerically and thus physically more interesting.

TODO:

Thread the parallel code using openMP

============================================================================
############################################################################

 
 
