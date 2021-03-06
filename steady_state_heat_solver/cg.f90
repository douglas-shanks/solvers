!
!  Description: Subroutine implementing the CG method (parallel version)
!
! -----------------------------------------------------------------------

subroutine cg(A,u,b,tol,maxits,its,rnorm)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Arguments:
!      
!       Input:
!
!         A          local part of system matrix in compressed row format
!         b          right hand side vector in distributed form
!         tol        tolerance for iterative method
!         maxits     maximum number of iterations
!        
!       Output:
!
!         u          solution vector in distibuted form 
!         its        number of iterations
!		  rnorm 	norm of residual
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  use header

  
  implicit none

  type(Matrix), intent(in)    :: A
  type(Vector), intent(inout) :: u
  real(kind=8), intent(out) :: rnorm
  type(Vector), intent(in)    :: b
  real(kind=8), intent(in) :: tol
  integer, intent(in)  :: maxits
  integer, intent(out) :: its

  real(kind=8) :: one, zero
  parameter (one = 1.0_8, zero = 0.0_8)

  type(Vector) :: p
  type(Vector) :: q
  type(Vector) :: r
  real(kind=8) :: rtr,alpha,beta,gamma,delta,norm,norm0
  integer      :: n_loc, ierr

  n_loc = b%iend - b%ibeg + 1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for additional vectors p, q and r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(p%xx(b%n))
  allocate(q%xx(b%n))
  allocate(r%xx(b%n))

  p%n    = b%n
  p%ibeg = b%ibeg
  p%iend = b%iend

  q%n    = b%n
  q%ibeg = b%ibeg
  q%iend = b%iend

  r%n    = b%n
  r%ibeg = b%ibeg
  r%iend = b%iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program - Initialise solution vector and other vectors
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx(u%ibeg:u%iend) = zero

  call dcopy(n_loc,b%xx(b%ibeg),1,r%xx(r%ibeg),1)
  call dcopy(n_loc,b%xx(b%ibeg),1,p%xx(p%ibeg),1)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate initial residual norm and stop if small enough
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Vec_Dot(r,r,rtr)

  norm0 = sqrt(rtr)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Iterate - up to kmax iterations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  do its=1,maxits

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Implementation of one CG iteration
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     call Mat_Mult(A,p,q)

     call Vec_Dot(p,q,gamma)

     alpha = rtr / gamma

     call daxpy(n_loc,alpha,p%xx(p%ibeg),1,u%xx(u%ibeg),1)
     call daxpy(n_loc,-alpha,q%xx(q%ibeg),1,r%xx(r%ibeg),1)

     call Vec_Dot(r,r,delta)

     beta = delta / rtr
     rtr  = delta

     norm = sqrt(rtr)
     rnorm = norm
     if (norm/norm0 < tol) exit

     call dscal(n_loc,beta,p%xx(p%ibeg),1)
     call daxpy(n_loc,one,r%xx(r%ibeg),1,p%xx(p%ibeg),1)

  end do

  deallocate(p%xx)
  deallocate(q%xx)
  deallocate(r%xx)

end subroutine cg

