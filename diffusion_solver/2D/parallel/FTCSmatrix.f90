subroutine FTCSmatrix(A,alpha,tmax,m,Nt,ibeg,iend)

! Subroutine to form matrix arriving from 
! FTCS for the heat equation
! CRS version

  use header

  type(Matrix) :: A
  integer, intent(in) :: m, Nt,ibeg,iend
  real(kind=8), intent(in) :: alpha, tmax
  
  real(kind=8):: hx, ht, R1
  integer :: row,col,i,j,inz
  real(kind=8) :: v
  
  ! interface for the diffusion coefficient
!  interface
!    function alpha(x,y) result (val)

!    use header

!    real(kind=8), intent(in) :: x,y
!    real(kind=8) :: val

!    end function alpha
!  end interface
  
! Mesh spacing for time and space

  hx = 1.0_8/real(m,8)
  ht = tmax/real(Nt-1,8)
  R1 = alpha*ht/(hx**2)

  inz = 0
  do row=ibeg,iend

! Set the diagonal entry and off-diagonal entries in row. To do this we 
    
     j = (row-1)/(m-1) + 1
     i = row - (j-1)*(m-1) 
     inz = inz + 1

     A%aa(inz) = 1.0_8 - 4.0_8*R1
     A%ii(row) = inz
     A%jj(inz) = row

! If i = 1 then there is no (geometric) left neighbour, thus no entry

     if (i > 1) then
        col = row - 1  
        inz = inz + 1
     
        A%aa(inz) =R1
        A%jj(inz) = col
     endif

! If i = m-1 then there is no (geometric) right neighbour, thus no entry

     if (i < m-1) then
        col = row + 1  
        inz = inz + 1
     
        A%aa(inz) = R1
        A%jj(inz) = col
     endif

! If j = 1 then there is no (geometric) lower neighbour, thus no entry

     if (j > 1) then
        col = row - m + 1 
        inz = inz + 1
     
        A%aa(inz) = R1
        A%jj(inz) = col
     endif

! If j = m-1 then there is no (geometric) upper neighbour, thus no entry

     if (j < m-1) then
        col = row + m - 1 
        inz = inz + 1
     
        A%aa(inz) = R1
        A%jj(inz) = col
     endif
  end do
 
! Set A%n (the dimension of A), A%nnz (the number of nonzero entries) and
! A%ii(A%n+1) which is needed to address the nonzero entries in the last
! row (in order to know where the arrays A%aa and A%jj end)

  A%bw   = m-1
  A%nnz  = inz
  A%ii(A%iend+1) = A%nnz + 1
  
end subroutine FTCSmatrix

