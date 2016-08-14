!==========================================================================
!
!  Subroutine to set the RHS for the 2-D, five-point stencil 
!  of a general 2nd-order elliptic PDE in compressed row storage format 
!  (parallel version).
!  The discretied PDE
! 		-\grad a(x,y) \grad U + sigma U = g,  in D
!									  U = 0,  on D
! Also options to add non zero boundary conditions in rhs.f90 function 
!==========================================================================

subroutine sysmatrix(A,m,ibeg,iend,sigma)

  use header

  type(Matrix), intent(inout) :: A
  integer, intent(in) :: m,ibeg,iend

  integer :: row,col,i,j,inz,c
  real(kind=8) :: xm,xp,ym,yp, sigma

  interface
    function alpha(x,y) result (val)

    use header

    real(kind=8), intent(in) :: x,y
    real(kind=8) :: val

    end function alpha
  end interface

  c = m*m

  inz = 0
  do row=ibeg,iend

! Set the diagonal entry and the off-diagonal entries in each row. 
! To do this we calculate the indices (i,j) in the cartesian numbering 
! of the unknowns and then set the entries accordingly
  
     j = (row-1)/(m-1) + 1
     i = row - (j-1)*(m-1)  

! Set the diagonal entry in row, i.e. the sum of the coefficient in the
! four mesh cells adjacent to the given node 

     xm = (i-0.5d0)/m
     xp = (i+0.5d0)/m
     ym = (j-0.5d0)/m
     yp = (j+0.5d0)/m

     inz = inz + 1

     A%aa(inz) = c * (alpha(xm,ym) + alpha(xp,ym) + alpha(xm,yp) + alpha(xp,yp)) + sigma
     A%ii(row) = inz
     A%jj(inz) = row

! If i = 1 then there is no (geometric) left neighbour, thus no entry

     if (i > 1) then
        col = row - 1  
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * c * (alpha(xm,ym) + alpha(xm,yp))
        A%jj(inz) = col
     endif

! If i = m-1 then there is no (geometric) right neighbour, thus no entry

     if (i < m-1) then
        col = row + 1  
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * c * (alpha(xp,ym) + alpha(xp,yp))
        A%jj(inz) = col
     endif

! If j = 1 then there is no (geometric) lower neighbour, thus no entry

     if (j > 1) then
        col = row - m + 1 
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * c * (alpha(xm,ym) + alpha(xp,ym))
        A%jj(inz) = col
     endif

! If j = m-1 then there is no (geometric) upper neighbour, thus no entry

     if (j < m-1) then
        col = row + m - 1
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * c * (alpha(xm,yp) + alpha(xp,yp))
        A%jj(inz) = col
     endif
  end do
 
! Set A%n (the dimension of A), A%nnz (the number of nonzero entries) and
! A%ii(A%n+1) which is needed to address the nonzero entries in the last
! row (in order to know where the arrays A%aa and A%jj end)

  A%bw   = m-1
  A%nnz  = inz
  A%ii(iend+1) = A%nnz + 1

end subroutine sysmatrix
