!
!  Description: Subroutine implementing Jacobi's method 
!               (parallel version)
!
! -----------------------------------------------------------------------

subroutine Jacobi(A,u,b,eps,kmax,its,myid)

  use header
  include "mpif.h"

  type(Matrix), intent(inout) :: A
  type(Vector), intent(inout) :: u
  type(Vector), intent(inout) :: b
  real(kind=8), intent(in) :: eps
  integer, intent(in) :: kmax,myid
  integer, intent(out) :: its

  type(Vector) :: u_old
  real(kind=8) :: norm
  integer :: i,i_j,j

  allocate(u_old%xx(u%n))
  u_old%n    = u%n
  u_old%ibeg = u%ibeg
  u_old%iend = u%iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program - Initialise the solution vector
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx(u%ibeg:u%iend) = 0.0_8

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Iterate - up to kmax iterations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  do its=1,kmax

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Gather the entire vector u_k in u_old on each processor and
!     start calculating u_{k+1} by setting it to b (locally).
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     u_old%xx(u%ibeg:u%iend) = u%xx(u%ibeg:u%iend)
     u%xx(u%ibeg:u%iend)     = b%xx(u%ibeg:u%iend)

     call sparsegather(u_old,A%bw)

     do i=u%ibeg,u%iend
        do i_j = A%ii(i)+1, A%ii(i+1)-1

           j = A%jj(i_j)
           u%xx(i) = u%xx(i) - A%aa(i_j) * u_old%xx(j)
           
        end do
        u%xx(i) = u%xx(i) / A%aa(A%ii(i))      
     end do
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate || u_{k+1} - u_k ||_2 and check the stopping criterion   
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     call Error_jac(u,u_old,norm)
     !if (myid == 0) write(*,*) its, norm
     if (norm < eps) exit

  end do

end subroutine Jacobi








