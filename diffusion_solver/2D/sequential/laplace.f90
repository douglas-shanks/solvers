!
!  Subroutine to set matrix elements for the 2-D, five-point stencil 
!  in compressed row storage format.
!
!  Note: this uses the standard ordering of the unknowns as presented on
!        the handout in Week 5; ordering the indices (i,j) in the order 
!
!    (1,1),(2,1),(3,1),...,(m-1,1),(1,2),(2,2),...,(m-1,2),...,(m-1,m-1).
!
!------------------------------------------------------------------------

subroutine Laplace(A,m)

  use header

  type(Matrix) :: A
  integer, intent(in) :: m

  integer :: row,col,i,j,inz,c
  real(kind=8) :: v

  c = m*m

  inz = 0
  do row=1,(m-1)*(m-1)

! Set the diagonal entry in row

     v   = 4.0_8 * c
     inz = inz + 1

     A%aa(inz) = v
     A%ii(row) = inz
     A%jj(inz) = row

! Set the off-diagonal entries in row. To do this we calculate 
! the indices (i,j) in the cartesian numbering of the unknowns
! and then set the off-diagonal entries accordingly  

     v = -1.0_8 * c
     j = (row-1)/(m-1) + 1
     i = row - (j-1)*(m-1)  

! If i = 1 then there is no (geometric) left neighbour, thus no entry

     if (i > 1) then
        col = row - 1  
        inz = inz + 1
     
        A%aa(inz) = v
        A%jj(inz) = col
     endif

! If i = m-1 then there is no (geometric) right neighbour, thus no entry

     if (i < m-1) then
        col = row + 1  
        inz = inz + 1
     
        A%aa(inz) = v
        A%jj(inz) = col
     endif

! If j = 1 then there is no (geometric) lower neighbour, thus no entry

     if (j > 1) then
        col = row - m + 1 
        inz = inz + 1
     
        A%aa(inz) = v
        A%jj(inz) = col
     endif

! If j = m-1 then there is no (geometric) upper neighbour, thus no entry

     if (j < m-1) then
        col = row + m - 1 
        inz = inz + 1
     
        A%aa(inz) = v
        A%jj(inz) = col
     endif
  end do
 
! Set A%n (the dimension of A), A%nnz (the number of nonzero entries) and
! A%ii(A%n+1) which is needed to address the nonzero entries in the last
! row (in order to know where the arrays A%aa and A%jj end)

  A%nnz       = inz
  A%ii(A%n+1) = inz + 1

end subroutine Laplace







