Just a test bed for different interative solvers for
the steady state version of the linear diffusion PDE.

- MPI parallelism (tick, tested on my 4 core device)
- Jacobi (tick)
- CG     (tick)
- PCG (Jacobi diag)  (tick)
- PCG (Jacobi block)  (to do)
- GMG, V-cycle iterative solver (to do)
- CG-GMG (CG preconditioned with one V-cycle) (to do)

- Thread the whole problem using openMP (I have up to 8 threads, 2 per core :( )
