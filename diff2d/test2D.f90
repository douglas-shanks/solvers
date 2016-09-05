
!=============================================================================
!
!	A test program to check that 2D FTCS, BTCS and CN converges
!		Parallel Version
!
!=============================================================================

program test2D
	
  use header
  include "mpif.h"
  
  real (kind=8) :: eps, norm
  integer :: kmax, its
  parameter (eps = 1.0d-15, kmax = 10000)
  
  type(Matrix) ::  Abtcs, P
  type(Vector) ::  x, y, Uinitial, Ubtcs
  
  real(kind=8), dimension(:,:), allocatable :: Uprint
  real(kind=8) :: alpha, tmax, hx, ht, xx, yy, R1, pre_coeff
  real(kind=8) :: errorbtcs

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
  	WRITE(*,*) ' diff2d solves the linear diffusion PDE '
  	WRITE(*,*)
  	WRITE(*,*) ' dU/dt = grad(alpha grad U) in \Omega, '
  	WRITE(*,*) ' U = 0 on d\Omega '
  	WRITE(*,*) ' U(x,0) = g(x,y)'
  	WRITE(*,*)
  	WRITE(*,*) ' This is a test program where we solve the 2D heat equation '
  	WRITE(*,*) ' with the BTCS method using different linear solvers.'
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
  	WRITE(*,*)
  	WRITE(*,*) ' 0 = CG , 1 = PCG (Jacobi), 2 = Jacobi (slow) '
  	WRITE(*,*) flag
  	WRITE(*,*) '========================================================================='
  	WRITE(*,*) 	
  	IF (flag==0) then
  		write(*,*) 'We are using CG'
  	ELSEIF (flag==1) then
  		write(*,*) 'We are using PCG'
  	ELSE
  		write(*,*) 'We are using Jacobi'
  	ENDIF
  	WRITE(*,*)		
    WRITE(*,*) '=========================================================================' 
    WRITE(*,*)
      	WRITE(*,*) flag
 ! Error check 
  	if (m<0 .or. Nt<0) then
 	 	print *, "N and Nt must be positive "
 	 	stop
  	endif
     
  ! Mesh spacing for time and space

  	hx = 1.0_8/real(m-1,8)
  	ht = tmax/real(Nt-1,8)

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
  
  allocate(P%aa(5*nrows))
  allocate(P%jj(5*nrows))
  allocate(P%ii(N+1))
  
  allocate(x%xx(N))
  allocate(y%xx(N)) 

  allocate(Uinitial%xx(N))
  allocate(Ubtcs%xx(N))

  
  Abtcs%n    = N
  Abtcs%ibeg = ibeg
  Abtcs%iend = iend

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
  
  Uinitial%n = N
  Uinitial%ibeg = ibeg
  Uinitial%iend = iend
  Uinitial%xx(Uinitial%ibeg:Uinitial%iend) = 0.0d0
  
  Ubtcs%n = N
  Ubtcs%ibeg = ibeg
  Ubtcs%iend = iend
!  Ubtcs%xx(Ubtcs%ibeg:Ubtcs%iend) = 0.0d0
 
! Mesh spacing for time and space

  hx = 1.0_8/real(m-1,8)
  ht = tmax/real(Nt-1,8)
  R1 = alpha*ht/(hx**2)

!=============================================================================
!					initial condition
!=============================================================================

  !if (myid == 0) then
  	call initial2D(Uinitial,1.0_8,0.0_8,m)
  	call FDsolution(Uprint,Uinitial,m)

  	open(unit=2, file='Uinitial.txt', ACTION="write", STATUS="replace")
  	write(2,*)  m
  	write(2,*) 
  	write(2, *) Uprint
  	close(2)
  !end if
  
!=============================================================================
!						BTCS
!=============================================================================
 
! Initial solution, u_0

  call initial2D(Ubtcs,1.0_8,0.0_8,m)
  
! Make BTCS matrix  
  call BTCSmatrix(Abtcs,tmax,m,Nt,ibeg,iend) 
  
    its_total = 0   
	! now solve and produce output information 
	 if (flag == 0) then
	 
	 	call cpu_time(t1)
  		do j=2,Nt
  		
     	    call cg(Abtcs,y,Ubtcs,eps,kmax,its,norm)
     	    Ubtcs = y
     	    its_total = its_total + its
     	    write(*,*) 'At iterate'
     	    write(*,*) j
     	    write(*,*) 'Iterations taken'
     	    write(*,*) its
     	    write(*,*) ' 2-norm of the residual'
     	    write(*,*) norm 
     	    
  		end do
  		call cpu_time(t2)
  		
  	elseif (flag == 1) then
	    
	    ! Get the preconditioning matrix

	    call BTCSmatrix_PRE(P,tmax,m,Nt,ibeg,iend)

	 	call cpu_time(t1)
  		do j=2,Nt

     	    call pcg(Abtcs,y,Ubtcs,eps,kmax,its,P,norm)
     	    Ubtcs = y
     	    its_total = its_total + its
     	    write(*,*) 'At iterate'
     	    write(*,*) j
     	    write(*,*) 'Iterations taken'
     	    write(*,*) its
     	    write(*,*) ' 2-norm of the residual'
     	    write(*,*) norm      	     

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
!  call initial2D(Uinitial,1.0_8,tmax,m)
!  call error(Uinitial,Ubtcs,errorbtcs)
  
   if (myid == 0) then 
  
  	write(*,*)
  	write(*,*) '========================================================================='
  	print*,    ' BTCS Solution '
  	write(*,*) '========================================================================='
  	write(*,*)
  
  	write(*,*) ' residual '
 	write(*,*)
    write(*,*) ' 2-norm of the residual'
    write(*,*)	norm
  	write(*,*)
  	print*, ' time for solve'   
  	print'(f12.6)', t2-t1
  	write(*,*)
    print*, ' Solver iterations total'   
    print*, its_total
    print*, ' Average Solver iterations'
    print*, NINT(real(its_total/(Nt-1)))
  	
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
  
  deallocate(x%xx)
  deallocate(y%xx)
  
  deallocate(Uinitial%xx)
  deallocate(Ubtcs%xx)

  call MPI_Finalize(ierr) 

end program test2D
