!
!  Description: Subroutine implementing the CG method (sequential version)
!
! -----------------------------------------------------------------------

subroutine cg(A,u,b,tol,maxits,its)
  
  use header

  implicit none

  type(Matrix), intent(in) :: A
  type(Vector), intent(inout) :: u
  type(Vector), intent(in) :: b
  real(kind=8), intent(in) :: tol
  integer, intent(in) :: maxits
  integer, intent(out) :: its

  real(kind=8) :: one, zero
  parameter (one = 1.0_8, zero = 0.0_8)

  type(Vector) :: p
  type(Vector) :: q, r
  real(kind=8) :: rtr,alpha,beta,gamma,delta,norm,norm0
  integer      :: n

  n = b%n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for additional vectors p and q
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(p%xx(b%n))
  allocate(q%xx(b%n))
  allocate(r%xx(b%n))

  p%n = b%n
  q%n = b%n
  r%n = b%n
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program - Initialise solution vector and other vectors
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx = zero

  call dcopy(n,b%xx,1,r%xx,1)
  call dcopy(n,b%xx,1,p%xx,1)
  ! for preconditioner do p = inv(M)zzp
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate initial residual norm and stop if small enough
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  ! this changes if preconditioned (p,r,rtr) ?
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
     
     call daxpy(n,alpha,p%xx,1,u%xx,1)
     call daxpy(n,-alpha,q%xx,1,r%xx,1)
     
     ! we want inv(M)*r
     
     call Vec_Dot(r,r,delta)
     
     beta = delta / rtr
     rtr  = delta
     
     norm = sqrt(rtr)
     
     if (norm/norm0 < tol) exit
     
     call dscal(n,beta,p%xx,1)
     call daxpy(n,one,r%xx,1,p%xx,1)
     
  end do
 
  deallocate(p%xx)
  deallocate(q%xx)
  deallocate(r%xx)
  
end subroutine cg

