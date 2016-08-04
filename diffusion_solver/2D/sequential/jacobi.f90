!
!  Description: Subroutine implementing Jacobi's method 
!               (sequential version)
!
! -----------------------------------------------------------------------

subroutine Jacobi(A,u,b,eps,kmax,its)

  use header

  type(Matrix) :: A
  type(Vector) :: u
  type(Vector) :: b
  real(kind=8), intent(in) :: eps
  integer, intent(in) :: kmax
  integer, intent(out) :: its

  type(Vector) :: u_old
  real(kind=8) :: norm
  integer :: i,i_j,j

  allocate(u_old%xx(u%n))
  u_old%n = u%n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program - Initialise the solution vector
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx = 0.0_8

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Iterate - up to kmax iterations
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  do its=1,kmax

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Store u_k in u_old and start calculating u_{k+1} by setting it to b   
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     u_old%xx = u%xx
     u%xx     = b%xx

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     For each component of u_{k+1} do the Jacobi iteration Example 7.4(3)
!
!       Note that in the compressed row storage format the nonzero 
!       entries of row i are stored in 
!
!         A%aa(A%ii(i)), A%aa(A%ii(i)+1), ..., A%aa(A%ii(i+1)-1)
!
!       the according (global) column numbers are stored in
!
!         A%jj(A%ii(i)), A%jj(A%ii(i)+1), ..., A%jj(A%ii(i+1)-1)   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     do i=1,u%n
        do i_j = A%ii(i)+1, A%ii(i+1)-1

           j = A%jj(i_j)
           u%xx(i) = u%xx(i) - A%aa(i_j) * u_old%xx(j)
           
        end do
        u%xx(i) = u%xx(i) / A%aa(A%ii(i))      
     end do
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate || u_{k+1} - u_k ||_2 and check the stopping criterion   
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     call Error(u,u_old,norm)
     !write(*,*) its, norm
     if (norm < eps) exit

  end do

end subroutine Jacobi



