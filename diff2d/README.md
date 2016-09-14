#diff2d

Program to solve the linear diffusion PDE with variable diffusion 
coefficient.

- MPI parallelism (tick, tested on my 4 core device)
- Jacobi (tick)
- CG     (tick)
- PCG (Jacobi diag)  (tick)
- PCG (Jacobi block)  (to do)
- GMG, V-cycle iterative solver (to do)
- CG-GMG (CG preconditioned with one V-cycle) (to do)
- CG-AS with a coarse grid. Perhaps extend to some kind of multilevel method (i.e. >2 levels).

- Thread the whole problem using openMP
- 3D version?
