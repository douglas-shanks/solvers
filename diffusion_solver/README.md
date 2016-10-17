###########################################################################
===========================================================================

1D and 2D
--------------------

1D and 2D test programs to solve the heat equation
with Dirichlet boundary conditions and given initial condition.
The initial conditions given below allow one to have a closed
form solution for the analytical solution allowing the error to be 
computed.

1D:

u_0(x,0) = sin(\pi * x).

2D:

u_0(x,y,0) = sin(\pi * x)*sin(\pi * y).

To run:

make clean
make

./test1D

or 

./test2D

This generates the exact solution, FTCS solution, BTCS solution
and CN solution. The 2-norm of the error is then computed for each and wall clock time computed.

============================================================================
############################################################################

 
 
