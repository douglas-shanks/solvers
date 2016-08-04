!
!  Subroutine to multiply a matrix and a vector in compressed row storage 
!  format , i.e to calculate b = A*u.
!
!------------------------------------------------------------------------

subroutine Mat_Mult(A,u,b)

  use header

  type(Matrix) :: A
  type(Vector) :: u, b
  integer :: i,i_j,j

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Calculate each component of b by taking the scalar product of the
!     i-th row of A and the vector u.
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

  do i = 1,b%n

     b%xx(i) = 0.0_8
     do i_j = A%ii(i), A%ii(i+1)-1

        j = A%jj(i_j)
        b%xx(i) = b%xx(i) + A%aa(i_j) * u%xx(j)

     end do

  end do

end subroutine Mat_Mult
