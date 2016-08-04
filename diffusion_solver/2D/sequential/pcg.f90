!
!  Description: Subroutine implementing the CG method (sequential version)
!
! -----------------------------------------------------------------------

subroutine pcg(A,u,b,tol,maxits,its,M)
  
  use header

  implicit none

  type(Matrix), intent(in) :: A,M
  type(Vector), intent(inout) :: u
  type(Vector), intent(in) :: b
  real(kind=8), intent(in) :: tol
  integer, intent(in) :: maxits
  integer, intent(out) :: its

  real(kind=8) :: one, zero
  parameter (one = 1.0_8, zero = 0.0_8)

  type(Vector) :: p
  type(Vector) :: q, r, z
  real(kind=8) :: rtr,alpha,beta,gamma,delta,norm,norm0
  integer      :: n

  n = b%n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for additional vectors p and q
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(p%xx(b%n))
  allocate(q%xx(b%n))
  allocate(r%xx(b%n))
  allocate(z%xx(b%n))

  p%n = b%n
  q%n = b%n
  r%n = b%n
  z%n = b%n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program - Initialise solution vector and other vectors
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx = zero
  z%xx = zero
  
  call dcopy(n,b%xx,1,r%xx,1)
  
  ! z = M^{-1}r
  call Mat_Mult(M,r,z)
  
  ! p = z
  call dcopy(n,z%xx,1,p%xx,1)
  
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate initial residual norm and stop if small enough
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  ! rtr = r'*z
  call Vec_Dot(r,z,rtr)

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
     
     ! z = M^{-1}r
     call Mat_Mult(M,r,z)
     
     ! delta = z'*r
     call Vec_Dot(z,r,delta)
     
     beta = delta / rtr
     rtr  = delta
     
     norm = sqrt(rtr)
     
     if (norm/norm0 < tol) exit
     
     call dscal(n,beta,p%xx,1)
     
     ! p_{k+1} = z_{k+1} + b_{k}*p_{k}
     call daxpy(n,one,z%xx,1,p%xx,1)
     
  end do
 
  deallocate(p%xx)
  deallocate(q%xx)
  deallocate(r%xx)
  deallocate(z%xx)
  
end subroutine pcg

