
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
  integer :: kmax, its
  parameter (eps = 1.0d-16, kmax = 1000)
  
  type(Matrix) ::   Aftcs, Acn, Acnrhs, Abtcs, P
  type(Vector) ::  x, xcn, y, Uexact, Uftcs, Uftcs1, Ucn, Ubtcs
  
  real(kind=8), dimension(:,:), allocatable :: Uprint
  real(kind=8) :: alpha, tmax, hx, ht, xx, yy, R1, pre_coeff
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
 
  call MPI_Pack_size(2,MPI_INTEGER, &
                     & MPI_COMM_WORLD,size1,ierr);
  call MPI_Pack_size(2,MPI_DOUBLE_PRECISION, &
                     & MPI_COMM_WORLD,size2,ierr);

  size = size1 + size2
  position = 0
  
!=============================================================================

  if ( myid == 0 ) then
  
  	WRITE(*,*)
  	WRITE(*,*) '========================================================================='
  	WRITE(*,*)
  	WRITE(*,*) ' diff_solver solves the linear diffusion PDE '
  	WRITE(*,*)
  	WRITE(*,*) ' dU/dt = \alpha D^2 U in \Omega, '
  	WRITE(*,*) ' U = 0 on d\Omega '
  	WRITE(*,*) ' U(x,0) = sin(pi*x)sin(pi*y)'
  	WRITE(*,*)
  	WRITE(*,*) ' This is a test program where we solve the 2D heat equation '
  	WRITE(*,*) ' We compare the exact solution to that with the FTCS method '
  	WRITE(*,*) ' BTCS method and Crank Nicholson method.'
  	WRITE(*,*)
  	WRITE(*,*) ' This version is written in parallel using MPI.'
  	WRITE(*,*)
  	WRITE(*,*) '========================================================================='
  	WRITE(*,*)
  
  ! Read the inputs
  	open(unit=2,file="input.dat")
  	read(2,*) alpha,tmax,m,Nt,flag
  	
  	WRITE(*,*) ' Inputs:'
  	WRITE(*,*)
  	WRITE(*,*) ' alpha. '
  	WRITE(*,*) 	alpha
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
  	R1 = alpha*ht/(hx**2)
	pre_coeff = 1.0_8 / (1.0_8 + 4.0_8*R1)
	
  	print*, 'Value of CFL coefficient'
  	print*, 				R1
  	print*, ' '
  	
  	if (R1 >= 0.25_8) then
		 print *, "Cannot guarantee stability for FTCS "
  	endif
  
  	call MPI_Pack(alpha,1,MPI_DOUBLE_PRECISION, &
    	      & buffer,100,position,MPI_COMM_WORLD,ierr);
  	call MPI_Pack(m,1,MPI_INTEGER, &
    	      & buffer,100,position,MPI_COMM_WORLD,ierr);
  	call MPI_Pack(tmax,1,MPI_DOUBLE_PRECISION, &
    	      & buffer,100,position,MPI_COMM_WORLD,ierr);
  	call MPI_Pack(Nt,1,MPI_INTEGER, &
    	      & buffer,100,position,MPI_COMM_WORLD,ierr);
  end if 
         
!=============================================================================
!      Broadcast alpha, m, tmax and Nt to the other processes
!=============================================================================

    call MPI_Bcast(buffer,size,MPI_PACKED,0,MPI_COMM_WORLD,ierr);
    
!=============================================================================
!         Unpack to slaves
!=============================================================================

  if (myid > 0) then
     call MPI_Unpack(buffer,100,position,alpha,1, &
          & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr);
     call MPI_Unpack(buffer,100,position,m,1, &
          & MPI_INTEGER,MPI_COMM_WORLD,ierr);
     call MPI_Unpack(buffer,100,position,tmax,1, &
          & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr);
     call MPI_Unpack(buffer,100,position,Nt,1, &
          & MPI_INTEGER,MPI_COMM_WORLD,ierr);
  end if

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
  
  allocate(Aftcs%aa(5*nrows))
  allocate(Aftcs%jj(5*nrows))
  allocate(Aftcs%ii(N+1))
  
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

  allocate(Uexact%xx(N))
  allocate(Uftcs%xx(N))
  allocate(Uftcs1%xx(nrows))
  allocate(Ubtcs%xx(N))
  allocate(Ucn%xx(N))
  
  Aftcs%n    = N
  Aftcs%ibeg = ibeg
  Aftcs%iend = iend
  
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
  
  Uexact%n = N
  Uexact%ibeg = ibeg
  Uexact%iend = iend
  Uexact%xx(Uexact%ibeg:Uexact%iend) = 0.0d0
  
  Uftcs%n = N
  Uftcs%ibeg = ibeg
  Uftcs%iend = iend
  
  Uftcs1%n = nrows
  Uftcs1%ibeg = ibeg
  Uftcs1%iend = iend
  Uftcs%xx(Uftcs%ibeg:Uftcs%iend) = 0.0d0
  
  Ubtcs%n = N
  Ubtcs%ibeg = ibeg
  Ubtcs%iend = iend
!  Ubtcs%xx(Ubtcs%ibeg:Ubtcs%iend) = 0.0d0
  
  Ucn%n = N 
  Ucn%ibeg = ibeg
  Ucn%iend = iend
!  Ucn%xx(Ucn%ibeg:Ucn%iend) = 0.0d0
 
! Mesh spacing for time and space

  hx = 1.0_8/real(m-1,8)
  ht = tmax/real(Nt-1,8)
  R1 = alpha*ht/(hx**2)

!=============================================================================
!					EXACT SOLUTION
!=============================================================================

  !if (myid == 0) then
  	call exactsoln2D(Uexact,alpha,tmax,m)
  	call FDsolution(Uprint,Uexact,m)

  	open(unit=2, file='uexact.txt', ACTION="write", STATUS="replace")
  	write(2,*)  m
  	write(2,*) 
  	write(2, *) Uprint
  	close(2)
  !end if
  
!=============================================================================
!						FTCS
!=============================================================================
  
! Initial solution, u_0

  call exactsoln2D(Uftcs,alpha,0.0_8,m)

  call FTCSmatrix(Aftcs,alpha,tmax,m,Nt,ibeg,iend) 

! Then it should just be a case of U^{j} = A*U{j-1}
  ierr = 0
  call cpu_time(t1)
  do j=2,Nt

    x = Uftcs
	call Mat_Mult(Aftcs,x,y)
	Uftcs = y
	
  end do
  call cpu_time(t2)
  !  print out the solution
   Uprint = 0.0_8
  
   call FDsolution(Uprint,Uftcs,m)
  
   open(unit=2, file='uftcs.txt', ACTION="write", STATUS="replace")
   write(2,*) m
   write(2,*) 
   write(2, *) Uprint
   close(2)
      
! Compute the error

  call error(Uftcs,Uexact,errorftcs)

  if (myid == 0) then 

    write(*,*)
    write(*,*) '========================================================================='
    print*,    ' FTCS Solution '
    write(*,*) '========================================================================='

    write(*,*) 
    write(*,*) ' Error '
    write(*,*)
    print*, errorftcs	 
    write(*,*)
    print*, ' time for solve'   
    print'(f12.6)', t2-t1
  	
  end if 
  
!=============================================================================
!						BTCS
!=============================================================================
 
! Initial solution, u_0

  call exactsoln2D(Ubtcs,alpha,0.0_8,m)
  
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

	    call BTCSmatrix_PRE(P,alpha,tmax,m,Nt,ibeg,iend)

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
  call exactsoln2D(Uexact,alpha,tmax,m)
  call error(Uexact,Ubtcs,errorbtcs)
  
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
    print*, NINT(real(its_total/(Nt-1)))
  	
  	WRITE(*,*)
  	WRITE(*,*) '=========================================================================' 
 	WRITE(*,*) ' Program finished '
 	WRITE(*,*) '========================================================================='
 	
  end if

!=============================================================================
!				Crank Nicholson
!=============================================================================

! Initial solution, u_0
!  call exactsoln2D(Ucn,alpha,0.0_8,m)

! Make CN matrix and RHS 
!  call CNmatrix(Acn,alpha,tmax,m,Nt,ibeg,iend) 
  
! make the rhs
!  call CNrhs(Acnrhs,alpha,tmax,m,Nt,ibeg,iend) 
    
!    its_total = 0
	! now solve and produce output information
!     if (flag == 0 .Or. flag ==1) then
     
!     	call cpu_time(t1)
!  		do j=2,Nt
    
!        x = Ucn
!	    call Mat_Mult(Acnrhs,x,y)
!     	call cg(Acn,x,y,eps,kmax,its)
!     	Ucn = x
!		its_total = its_total + its 
		
!        end do
!        call cpu_time(t2)
     !elseif (flag ==1) then
     
      !  print*, 'Error not implemented'
        
!     elseif (flag == 2) then	
!        call cpu_time(t1)
!  		do j=2,Nt
    
!        x = Ucn
!	    call Mat_Mult(Acnrhs,x,y)
!     	call Jacobi(Acn,x,y,eps,kmax,its,myid)
!     	Ucn = x
!		its_total = its_total + its 
		
!        end do
!        call cpu_time(t2)
    	
! 	 endif

!  print out the solution

!  Uprint = 0.0_8    
!  call FDsolution(Uprint,Ucn,m)
  
!  open(unit=2, file='ucn.txt', ACTION="write", STATUS="replace")
!  write(2,*) m
!  write(2,*) 
!  write(2, *)( Uprint)
!  close(2) 
  
! Compute the error
!  call exactsoln2D(Uexact,alpha,tmax,m)
!  call error(Uexact,Ucn,errorcn)
  
!   if (myid == 0) then 
  
!  	write(*,*)
 ! 	write(*,*) '========================================================================='
!  	print*,    ' Crank Nicholson Solution '
!  	write(*,*) ' Does not currently have a preconditioner' 
!  	write(*,*)'========================================================================='
!  	write(*,*)
!  	write(*,*) ' Error '
! 	write(*,*)
!  	print*, errorcn	 
!  	write(*,*)
!  	print*, ' time for solve'   
!  	print'(f12.6)', t2-t1
!  	write(*,*)
!    print*, ' Solver iterations total'   
!    print*, its_total
!    print*, ' Average Solver iterations'
!    print*, NINT(real(its_total/(Nt-1)))
!  end if
  	

 


! Deallocate memory

  deallocate(Uprint)
  
  deallocate(Aftcs%aa)
  deallocate(Aftcs%jj)
  deallocate(Aftcs%ii)
  
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
  
  deallocate(Uexact%xx)
  deallocate(Uftcs%xx)
  deallocate(Ubtcs%xx)
  deallocate(Ucn%xx)

  call MPI_Finalize(ierr) 

end program test2D
