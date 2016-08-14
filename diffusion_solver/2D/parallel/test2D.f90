
!=============================================================================
!
!	A test program to check that 2D FTCS, BTCS and CN converges
!		Parallel Version
!
!=============================================================================

program test2D
	
  use header
  include "mpif.h"
  
  real (kind=8) :: eps
  integer :: kmax, its,its_total
  parameter (eps = 1.0d-16, kmax = 1000)
  
  type(Matrix) ::   Acn, Acnrhs, Abtcs, P
  type(Vector) ::  x, xcn, y, Ucn, Ubtcs
  
  real(kind=8), dimension(:,:), allocatable :: Uprint
  real(kind=8) :: tmax, hx, ht, xx, yy, R1, pre_coeff
  real(kind=8) ::errorftcs, errorcn, errorbtcs

  integer:: N, m, Nt, flag
  integer i, j
  
  ! Needed by MPI
  integer       :: myid, numprocs, nrows, ibeg, iend, irow
  integer       :: ierr,position,size,size1,size2
  character*(100) :: buffer
  
  ! For Timing
  real(kind=8) :: t1, t2

!=============================================================================
!      Beginning of program - Initialisation of MPI context
!=============================================================================


  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)	
  call greetings(myid,numprocs) 

  
!=============================================================================

  if ( myid == 0 ) then
  
  	WRITE(*,*)
  	WRITE(*,*) '========================================================================='
  	WRITE(*,*)
  	WRITE(*,*) ' diff_solver solves the linear diffusion PDE '
  	WRITE(*,*)
  	WRITE(*,*) ' dU/dt = D (\alpha D) U in \Omega, '
  	WRITE(*,*) ' U = 0 on d\Omega '
  	WRITE(*,*) ' U(x,0) = sin(pi*x)sin(pi*y)'
  	WRITE(*,*)
  	WRITE(*,*) ' This is a test program where we solve the 2D heat equation '
  	WRITE(*,*) ' We compare the exact solution to that with the BTCS method '
  	WRITE(*,*) '  and Crank Nicholson method.'
  	WRITE(*,*)
  	WRITE(*,*) ' This version is written in parallel using MPI.'
  	WRITE(*,*)
  	WRITE(*,*) '========================================================================='
  	WRITE(*,*)
  
  ! Read the inputs
  	open(unit=2,file="input.dat")
  	read(2,*) tmax,m,Nt,flag
  	
  	WRITE(*,*) ' Inputs:'
  	WRITE(*,*)
  	WRITE(*,*) ' tmax. '
  	WRITE(*,*) 	tmax
  	WRITE(*,*) ' Nx,Ny '
  	WRITE(*,*) 	m
  	WRITE(*,*) ' Nt '
  	WRITE(*,*) 	Nt
  	WRITE(*,*) 
  	WRITE(*,*) ' Which solver for BTCS and CN. ' 
  	WRITE(*,*) ' 0 = CG , 1 = PCG (Jacobi), 2 = Jacobi (slow) '
  	WRITE(*,*) flag 
     
 ! Error check 
  	if (m<0 .or. Nt<0) then
 	 	print *, "N and Nt must be positive "
 	 	stop
  	endif
     
  ! Mesh spacing for time and space

  	hx = 1.0_8/real(m-1,8)
  	ht = tmax/real(Nt-1,8)
  	R1 = ht/(hx**2)
	pre_coeff = 1.0_8 / (1.0_8 + 4.0_8*R1)
	
  	print*, 'Value of CFL coefficient'
  	print*, 				R1
  	print*, ' '
  	
  	if (R1 >= 0.25_8) then
		 print *, "Cannot guarantee stability for FTCS "
  	endif

  end if 
         
!=============================================================================
!      Broadcast alpha, m, tmax and Nt to the other processes
!=============================================================================
 ! I should just pack this all in a buffer!
 
  call MPI_Bcast(tmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Nt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (m==0) stop


! Total grid points, not including Dirichlet boundary

  N = (m-1)**2
  
!=============================================================================
!      Calculate the start and end indices of the rows to be held locally
!=============================================================================

  nrows = N / numprocs
  ibeg  = myid * nrows + 1
  
  if (myid < numprocs-1) then
     iend = (myid+1) * nrows
  else
     iend  = N
     nrows = N - myid*nrows
  end if

! Allocate memory

  allocate(Uprint(m-1,m-1))	
  
  allocate(Abtcs%aa(5*nrows))
  allocate(Abtcs%jj(5*nrows))
  allocate(Abtcs%ii(N+1))
  
  allocate(Acn%aa(5*nrows))
  allocate(Acn%jj(5*nrows))
  allocate(Acn%ii(N+1))
  
  allocate(Acnrhs%aa(5*nrows))
  allocate(Acnrhs%jj(5*nrows))
  allocate(Acnrhs%ii(N+1))
  
  allocate(P%aa(5*nrows))
  allocate(P%jj(5*nrows))
  allocate(P%ii(N+1))
  
  allocate(x%xx(N))
  allocate(y%xx(N))
  allocate(xcn%xx(N))  

  allocate(Ubtcs%xx(N))
  allocate(Ucn%xx(N))
  
  Abtcs%n    = N
  Abtcs%ibeg = ibeg
  Abtcs%iend = iend
  
  Acn%n    = N
  Acn%ibeg = ibeg
  Acn%iend = iend
  
  Acnrhs%n    = N
  Acnrhs%ibeg = ibeg
  Acnrhs%iend = iend

  P%n    = N
  P%ibeg = ibeg
  P%iend = iend
    
  x%n    = N
  x%ibeg = ibeg
  x%iend = iend
  x%xx(x%ibeg:x%iend) = 0.0d0
  
  y%n    = N
  y%ibeg = ibeg
  y%iend = iend
  y%xx(y%ibeg:y%iend) = 0.0d0
  
  xcn%n  = N
  xcn%ibeg = ibeg
  xcn%iend = iend
  xcn%xx(xcn%ibeg:xcn%iend) = 0.0d0
  
  Ubtcs%n = N
  Ubtcs%ibeg = ibeg
  Ubtcs%iend = iend
!  Ubtcs%xx(Ubtcs%ibeg:Ubtcs%iend) = 0.0d0
  
  Ucn%n = N 
  Ucn%ibeg = ibeg
  Ucn%iend = iend
!  Ucn%xx(Ucn%ibeg:Ucn%iend) = 0.0d0
  
!=============================================================================
!						BTCS
!=============================================================================
 
! Initial solution, u_0
print*, 'jere'
  alpha = 1.0_8
  call exactsoln2D(Ubtcs,alpha,0.0_8,m)
  print*, 'jere'
! Make BTCS matrix  
  call BTCSmatrix(Abtcs,alpha,tmax,m,Nt,ibeg,iend) 
  
   its_total = 0   
	! now solve and produce output information 
	 if (flag == 0) then
	 
	 	call cpu_time(t1)
  		do j=2,Nt
  		
     	    call cg(Abtcs,y,Ubtcs,eps,kmax,its)
     	    Ubtcs = y
     	    its_total = its_total + its 
     	    
  		end do
  		call cpu_time(t2)
  		
  	elseif (flag == 1) then
	    
	    ! Get the preconditioning matrix, problem with MPI
print*, 'jere'
	    call BTCSmatrix_PRE(P,alpha,tmax,m,Nt,ibeg,iend)
print*, 'jere'
	 	call cpu_time(t1)
  		do j=2,Nt

     	    call pcg(Abtcs,y,Ubtcs,eps,kmax,its,P,pre_coeff)
     	    Ubtcs = y
     	    its_total = its_total + its 

  		end do
  		call cpu_time(t2)
  		
     elseif (flag == 2) then	
	 
	 	call cpu_time(t1)
  		do j=2,Nt
  		
     	    call Jacobi(Abtcs,y,Ubtcs,eps,kmax,its,myid)
     	    Ubtcs = y
    		its_total = its_total + its 
    		
  		end do
  		call cpu_time(t2)
  		
 	 endif

!  print out the solution
  Uprint = 0.0_8    
  call FDsolution(Uprint,Ubtcs,m)
  
  open(unit=2, file='ubtcs.txt', ACTION="write", STATUS="replace")
  write(2,*) m
  write(2,*) 
  write(2, *)( Uprint)
  close(2) 
! Compute the error
  !call exactsoln2D(Uexact,alpha,tmax,m)
  !call error(Uexact,Ubtcs,errorbtcs)
  
   if (myid == 0) then 
  
  	write(*,*)
  	write(*,*) '========================================================================='
  	print*,    ' BTCS Solution '
  	write(*,*) '========================================================================='
  	write(*,*)
  
  	write(*,*) ' Error '
 	 write(*,*)
  	print*, errorbtcs	 
  	write(*,*)
  	print*, ' time for solve'   
  	print'(f12.6)', t2-t1
  	write(*,*)
    print*, ' Solver iterations total'   
    print*, its_total
    print*, ' Average Solver iterations'
    print*, ((its_total/(Nt)))
  	
  	WRITE(*,*)
  	WRITE(*,*) '=========================================================================' 
 	WRITE(*,*) ' Program finished '
 	WRITE(*,*) '========================================================================='
 	
  end if

! Deallocate memory

  deallocate(Uprint)
  
  deallocate(Abtcs%aa)
  deallocate(Abtcs%jj)
  deallocate(Abtcs%ii)
  
  deallocate(Acn%aa)
  deallocate(Acn%jj)
  deallocate(Acn%ii)
  
  deallocate(Acnrhs%aa)
  deallocate(Acnrhs%jj)
  deallocate(Acnrhs%ii)
  
  deallocate(x%xx)
  deallocate(y%xx)
  deallocate(xcn%xx)
  
  deallocate(Ubtcs%xx)
  deallocate(Ucn%xx)

  call MPI_Finalize(ierr) 

end program test2D
