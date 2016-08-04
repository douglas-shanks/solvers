!=============================================================
!
! Driver for the solution of  he discretied PDE
! 		-\grad a(x,y) \grad U + sigma U = g,  in D
!									  U = 0,  on D
! 
! and assortment of linear solvers are used
!=============================================================

program test_steady

  use header
  include "mpif.h"

  real (kind=8) :: eps
  integer :: kmax
  parameter (eps = 1.0d-8, kmax = 9999)
  real(kind=8) t1, t2; 

  type(Matrix)  :: A, P
  type(Vector)  :: u, u_ex, b
  
  real (kind=8) :: norm
  real(kind=8), dimension(:,:), allocatable :: Uprint
  
  integer       :: m
  integer       :: myid, nprocs, nrows, ibeg, iend
  integer       :: ierr
  integer		:: position,size,size1,size2
  character*(100) :: buffer
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Beginning of program - Initialisation of MPI context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

 
  call MPI_Pack_size(1,MPI_INTEGER, &
                     & MPI_COMM_WORLD,size1,ierr);
  call MPI_Pack_size(1,MPI_DOUBLE_PRECISION, &
                     & MPI_COMM_WORLD,size2,ierr);

  size = size1 + size2
  position = 0
  
  !call greetings(myid,nprocs)

  if ( myid == 0 ) then
     print*,'Number of discretisation points per coordinate direction:'
     open(unit=2,file="input.dat")
     read(2,*) m, sigma
     print*, 'Value of m =',  m
	 print*, ' Value of sigma=', sigma
!     if (mod((m-1)*(m-1),nprocs) /= 0) then
!        print*,'(m-1)*(m-1) has to be a multiple of nprocs =',nprocs
!        m = 0
!     end if

     close(2)
     
     call MPI_Pack(sigma,1,MPI_DOUBLE_PRECISION, &
          & buffer,100,position,MPI_COMM_WORLD,ierr);
          
     call MPI_Pack(m,1,MPI_INTEGER, &
          & buffer,100,position,MPI_COMM_WORLD,ierr);
          
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Broadcast m to the other processes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   call MPI_Bcast(buffer,size,MPI_PACKED,0,MPI_COMM_WORLD,ierr);
   
  !call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  !call MPI_Bcast(sigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  !if (m==0) stop
  
  if (myid > 0) then
  
     call MPI_Unpack(buffer,100,position,sigma,1, &
          & MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr);
     call MPI_Unpack(buffer,100,position,m,1, &
          & MPI_INTEGER,MPI_COMM_WORLD,ierr);
          
  end if
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate the start and end indices of the rows to be held locally
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  n     = (m-1)*(m-1)

  nrows = n / nprocs
  ibeg  = myid * nrows + 1
  iend  = (myid+1) * nrows

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for A, u, u_ex and b and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(Uprint(m-1,m-1))	
  
  allocate(A%aa(5*nrows))
  allocate(A%jj(5*nrows))
  allocate(A%ii(n+1))

  A%n    = n
  A%ibeg = ibeg
  A%iend = iend

  allocate(P%aa(5*nrows))
  allocate(P%jj(5*nrows))
  allocate(P%ii(n+1))

  P%n    = n
  P%ibeg = ibeg
  P%iend = iend

  allocate(u%xx(n))
  allocate(b%xx(n))
  allocate(u_ex%xx(n))

  b%n    = n
  b%ibeg = ibeg
  b%iend = iend

  u%n    = n
  u%ibeg = ibeg
  u%iend = iend

  u_ex%n    = n
  u_ex%ibeg = ibeg
  u_ex%iend = iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the matrix A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call sysmatrix(A,m,ibeg,iend,sigma) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Set exact solution u_ex to a random vector and calculate rhs b = A*u
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call random_number(u_ex%xx(ibeg:iend))
  
  !call rhs(b,m)
  
  call Mat_Mult(A,u_ex,b)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Apply Jacobi's method to solve the system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 t1 = MPI_Wtime()
 
 call Jacobi(A,u,b,eps,kmax,its,myid)
 
 t2 = MPI_Wtime()

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Check the error
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Error(u,u_ex,norm)
  
  if (myid == 0) then
     if (its > kmax) then
        write(*,100) kmax,norm
     else
        write(*,110) its,norm
     endif

100  format('Maximum number of iterations (',i5,          &
    &       ') reached, norm of the error is ',e10.4)
110  format('After ',i5,' Jacobi iterations the norm of the error is ',e10.4)
	 
	 write(6,*)  "Elapsed time (s) is ", t2 - t1 
  end if
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Apply CG to solve the system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

 t1 = MPI_Wtime()
 
 call cg(A,u,b,eps,kmax,its)
 
 t2 = MPI_Wtime()

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Check the error for CG
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Error(u,u_ex,norm)
  
  if (myid == 0) then
     if (its > kmax) then
        write(*,200) kmax,norm
     else
        write(*,210) its,norm
     endif

200  format('Maximum number of iterations (',i5,          &
    &       ') reached, norm of the error is ',e10.4)
210  format('After ',i5,' CG iterations the norm of the error is ',e10.4)
	 
	 write(6,*)  "Elapsed time (s) is ", t2 - t1 
  end if
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Apply PCG to solve the system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
 call sysmatrix_pre(P,m,ibeg,iend,sigma,P)
 
 t1 = MPI_Wtime()
 
 call pcg(A,u,b,eps,kmax,its,P)
 
 t2 = MPI_Wtime()

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Check the error for PCG
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Error(u,u_ex,norm)
  
  if (myid == 0) then
     if (its > kmax) then
        write(*,300) kmax,norm
     else
        write(*,310) its,norm
     endif

300  format('Maximum number of iterations (',i5,          &
    &       ') reached, norm of the error is ',e10.4)
310  format('After ',i5,' PCG iterations the norm of the error is ',e10.4)
	 
	 write(6,*)  "Elapsed time (s) is ", t2 - t1 
  end if 	
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Write the solution to a file for postprocessing in Matlab
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  if (nprocs == 1) then
  	  	call FDsolution(Uprint,u,m)

  		open(unit=2, file='solution.txt', ACTION="write", STATUS="replace")
  		write(2,*)  m
  		write(2,*) 
  		write(2, *) Uprint
  		close(2)
  end if		

  deallocate(Uprint)
  deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(P%aa)
  deallocate(P%jj)
  deallocate(P%ii)
  deallocate(u_ex%xx)
  deallocate(u%xx)
  deallocate(b%xx)

  call MPI_Finalize(ierr)

end program test_steady




