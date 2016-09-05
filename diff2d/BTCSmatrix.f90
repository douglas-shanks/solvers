! Subroutine to form matrix arriving from 
! BTCS for the heat equation
! CRS version

subroutine BTCSmatrix(A,tmax,m,Nt,ibeg,iend)

  use header

  type(Matrix) :: A
  integer, intent(in) :: m, Nt, ibeg, iend
  real(kind=8), intent(in) :: tmax
  
  real(kind=8):: hx, ht, R
  integer :: row,col,i,j,inz
  real(kind=8) :: v
  
  real(kind=8) :: xm,xp,ym,yp

  interface
    function alpha(x,y) result (val)

    use header

    real(kind=8), intent(in) :: x,y
    real(kind=8) :: val

    end function alpha
  end interface
  
! Mesh spacing for time and space

  hx = 1.0_8/real(m,8)
  ht = tmax/real(Nt-1,8)
  R = ht/(hx**2)

  A%bw = m+1

  inz = 0
  do row=ibeg,iend

! Set the diagonal entry and off-diagonal entries in row. 

     j = (row-1)/(m-1) + 1
     i = row - (j-1)*(m-1)  
     
     xm = (i-0.5d0)/m
     xp = (i+0.5d0)/m
     ym = (j-0.5d0)/m
     yp = (j+0.5d0)/m
     
     inz = inz + 1

     A%aa(inz) = 1.0_8 + R*(alpha(xm,ym) + alpha(xp,ym) + alpha(xm,yp) + alpha(xp,yp))
     A%ii(row) = inz
     A%jj(inz) = row



! If i = 1 then there is no (geometric) left neighbour, thus no entry

     if (i > 1) then
        col = row - 1  
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * R * (alpha(xm,ym) + alpha(xm,yp))
        A%jj(inz) = col
     endif

! If i = m-1 then there is no (geometric) right neighbour, thus no entry

     if (i < m-1) then
        col = row + 1  
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * R * (alpha(xp,ym) + alpha(xp,yp))
        A%jj(inz) = col
     endif

! If j = 1 then there is no (geometric) lower neighbour, thus no entry

     if (j > 1) then
        col = row - m + 1 
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * R * (alpha(xm,ym) + alpha(xp,ym))
        A%jj(inz) = col
     endif

! If j = m-1 then there is no (geometric) upper neighbour, thus no entry

     if (j < m-1) then
        col = row + m - 1 
        inz = inz + 1
     
        A%aa(inz) = -0.5d0 * R * (alpha(xm,yp) + alpha(xp,yp))
        A%jj(inz) = col
     endif
  end do
 
! Set A%n (the dimension of A), A%nnz (the number of nonzero entries) and
! A%ii(A%n+1) which is needed to address the nonzero entries in the last
! row (in order to know where the arrays A%aa and A%jj end)

  A%bw   = m-1
  A%nnz  = inz
  A%ii(A%iend+1) = A%nnz + 1
  
end subroutine BTCSmatrix

