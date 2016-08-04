
!---------------------------------------------------------------------
!
!	A test program to check that 2D FTCS, BTCS and CN converges
!
!---------------------------------------------------------------------
program test2D
	
  use header
  
  real (kind=8) :: eps
  integer :: kmax, its, its_total
  parameter (eps = 1.0d-16, kmax = 1000)
  
  type(Matrix) ::   Aftcs, Acn, Acnrhs, Abtcs, Mbtcs, Mcn
  type(Vector) ::  x, xcn, y, Uexact, Uftcs, Ucn, Ubtcs
  
  real(kind=8), dimension(:,:), allocatable :: Uprint
  real(kind=8) :: alpha, tmax, hx, ht, xx, yy, R1
  real(kind=8) ::errorftcs, errorcn, errorbtcs

  integer:: N, m, Nt, flag
  integer i, j
  real(kind=8) :: t1, t2

WRITE(*,*)
  WRITE(*,*) '========================================================================='
  WRITE(*,*)
  WRITE(*,*) ' dU/dt = \alpha D^2 U in \Omega, U = 0 on d\Omega '
  WRITE(*,*)
  WRITE(*,*) ' This is a test program where we solve the 2D heat equation '
  WRITE(*,*) ' We compare the exact solution to that with the FTCS method '
  WRITE(*,*) ' BTCS method and Crank Nicholson method.'
  WRITE(*,*) 'This version is sequential.'
  WRITE(*,*)
  WRITE(*,*) '========================================================================='
  WRITE(*,*)
  
  ! Read the inputs
  open(unit=2,file="input.dat")
  read(2,*) alpha,tmax,m,Nt,flag
  WRITE(*,*) 
  WRITE(*,*) ' Inputs:'
  WRITE(*,*)
  WRITE(*,*) ' alpha. '
  WRITE(*,*) 	alpha
  WRITE(*,*)   
  WRITE(*,*) ' tmax. '
  WRITE(*,*) 	tmax
  WRITE(*,*) 
  WRITE(*,*) ' Nx,Ny '
  WRITE(*,*) 	m
  WRITE(*,*) 
  WRITE(*,*) ' Nt '
  WRITE(*,*) 	Nt
  WRITE(*,*) 
  WRITE(*,*) ' Which solver for BTCS and CN. ' 
  WRITE(*,*) ' 0 = CG , 1 = PCG (Jacobi), 2 = Jacobi (slow) '
  WRITE(*,*) flag 
  WRITE(*,*) 
       
! Total grid points, not including Dirichlet boundary
  N = (m-1)**2
  
 ! Error check 
  if (m<0 .or. Nt<0) then
 	 print *, "N and Nt must be positive "
 	 stop
  endif

! Allocate memory

  allocate(Uprint(m-1,m-1))	
  
  allocate(Aftcs%aa(5*n))
  allocate(Aftcs%jj(5*N))
  allocate(Aftcs%ii(N+1))
  
  allocate(Abtcs%aa(5*N))
  allocate(Abtcs%jj(5*N))
  allocate(Abtcs%ii(N+1))
    
  allocate(Mbtcs%aa(5*N))
  allocate(Mbtcs%jj(5*N))
  allocate(Mbtcs%ii(N+1))
  
  allocate(Acn%aa(5*N))
  allocate(Acn%jj(5*N))
  allocate(Acn%ii(N+1))
  
  allocate(Mcn%aa(5*N))
  allocate(Mcn%jj(5*N))
  allocate(Mcn%ii(N+1))
  
  allocate(Acnrhs%aa(5*N))
  allocate(Acnrhs%jj(5*N))
  allocate(Acnrhs%ii(N+1))
  
  allocate(x%xx(N))
  allocate(y%xx(N))
  allocate(xcn%xx(N))  

  allocate(Uexact%xx(N))
  allocate(Uftcs%xx(N))
  allocate(Ubtcs%xx(N))
  allocate(Ucn%xx(N))
  
  Aftcs%n    = N
  Abtcs%n    = N
  Mbtcs%n    = N
  Acn%n    = N
  Mcn%n    = N
  Acnrhs%n    = N
  x%n    = N
  y%n    = N
  xcn%n  = N
  Uexact%n = N
  Uftcs%n = N
  Ubtcs%n = N
  Ucn%n = N 
 
 
! Mesh spacing for time and space

  hx = 1.0_8/real(m-1,8)
  ht = tmax/real(Nt-1,8)
  R1 = alpha*ht/(hx**2)

  print*, 'Value of CFL coefficient'
  print*, 				R1
  print*, ' '
  if (R1 >= 0.25_8) then
	 print *, "Cannot guarantee stability for FTCS "
  endif

!-----------------------------------------------------------
!					EXACT SOLUTION
!-----------------------------------------------------------

  call exactsoln2D(Uexact,alpha,tmax,m)

  call FDsolution(Uprint,Uexact,m)

  open(unit=2, file='uexact.txt', ACTION="write", STATUS="replace")
  write(2,*)  m
  write(2,*) 
  write(2, *) Uprint
  close(2)
  
!-----------------------------------------------------------
!						FTCS
!-----------------------------------------------------------

  write(*,*)
  write(*,*)  '=================================================================='
  print*,     'FTCS Solution'
  write(*,*)  '=================================================================='
  write(*,*) 
  
! Initial solution, u_0

  call exactsoln2D(Uftcs,alpha,0.0_8,m)

  call FTCSmatrix(Aftcs,alpha,tmax,m,Nt)
  
! Then it should just be a case of U^{j} = A*U{j-1}

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

  write(*,*) ' Error '
  write(*,*)
  print*, errorftcs	 
  write(*,*)
  print*, ' time for solve'   
  print'(f12.6)', t2-t1

!-----------------------------------------------------------
!						BTCS
!-----------------------------------------------------------

  write(*,*)
  write(*,*)  '=================================================================='
  print*,     'BTCS Solution'
  write(*,*)  '=================================================================='
  write(*,*)
 
! Initial solution, u_0

  call exactsoln2D(Ubtcs,alpha,0.0_8,m)
  
! Make BTCS matrix

  call BTCSmatrix(Abtcs,alpha,tmax,m,Nt)
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
  		
  	else if (flag == 1) then

  		! Make the preconditioner
		call BTCSmatrix_PRE(Mbtcs,alpha,tmax,m,Nt)

	 	call cpu_time(t1)
  		do j=2,Nt
  		
     	    call pcg(Abtcs,y,Ubtcs,eps,kmax,its,Mbtcs)
     	    Ubtcs = y
     	    its_total = its_total + its 
  		end do
  		call cpu_time(t2)
  		
     elseif (flag == 2) then	
	 
	 	call cpu_time(t1)
  		do j=2,Nt
  		
     	    call Jacobi(Abtcs,y,Ubtcs,eps,kmax,its)
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

  call error(Uexact,Ubtcs,errorbtcs)

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


!-----------------------------------------------------------
!				Crank Nicholson
!-----------------------------------------------------------

! Initial solution, u_0
 
  call exactsoln2D(Ucn,alpha,0.0_8,m)
! Make CN matrix and RHS 

  call CNmatrix(Acn,alpha,tmax,m,Nt)
! make the rhs

  call CNrhs(Acnrhs,alpha,tmax,m,Nt)
    its_total = 0
	! now solve and produce output information
     if (flag == 0) then
     
     	call cpu_time(t1)
  		do j=2,Nt
    
        x = Ucn
	    call Mat_Mult(Acnrhs,x,y)
     	call cg(Acn,x,y,eps,kmax,its)
     	Ucn = x
		its_total = its_total + its  
        end do
        call cpu_time(t2)
        
  	else if (flag == 1) then

  		! Make the preconditioner
		call CNmatrix_PRE(Mcn,alpha,tmax,m,Nt)

	 	call cpu_time(t1)
  		do j=2,Nt
  		    
  		    x = Ucn
	   		call Mat_Mult(Acnrhs,x,y)
     	    call pcg(Acn,x,y,eps,kmax,its,Mcn)
     	    Ucn = x
     	    its_total = its_total + its 
  		end do
  		call cpu_time(t2)
  		        
     elseif (flag == 2) then	
        call cpu_time(t1)
  		do j=2,Nt
    
        x = Ucn
	    call Mat_Mult(Acnrhs,x,y)
     	call Jacobi(Acn,x,y,eps,kmax,its)
     	Ucn = x
		its_total = its_total + its 
        end do
        call cpu_time(t2)
     	
 	 endif

  write(*,*)
  write(*,*)  '=================================================================='
  print*,     'Crank Nicholson Solution'
  write(*,*)  '=================================================================='
  write(*,*)

!  print out the solution

  Uprint = 0.0_8
  
  call FDsolution(Uprint,Ucn,m)
  
  open(unit=2, file='ucn.txt', ACTION="write", STATUS="replace")
  write(2,*) m
  write(2,*) 
  write(2, *)( Uprint)
  close(2)
  
! Compute the error

  call error(Uexact,Ucn,errorcn)

  write(*,*) ' Error '
  write(*,*)
  print*, errorcn	 
  write(*,*)
  print*, ' time for solve'   
  print'(f12.6)', t2-t1
  write(*,*)
  print*, ' Solver iterations total'   
  print*, its_total
  print*, ' Average Solver iterations'
  print*, NINT(real(its_total/(Nt-1)))

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


 WRITE(*,*)  '=================================================================='
 WRITE(*,*)  'Program finished'
 WRITE(*,*)  '=================================================================='

end program test2D
