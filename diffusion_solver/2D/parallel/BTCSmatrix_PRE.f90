! Subroutine to form jacobi preconditioner matrix arriving from 
! BTCS for the heat equation
! CRS version

subroutine BTCSmatrix_PRE(A,alpha,tmax,m,Nt,ibeg,iend)

  use header

  type(Matrix) :: A
  integer, intent(in) :: m, Nt, ibeg, iend
  real(kind=8), intent(in) :: alpha, tmax
  
  real(kind=8):: hx, ht, R1
  integer :: row,col,i,j,inz
  real(kind=8) :: v
  
! Mesh spacing for time and space

  hx = 1.0_8/real(m-1,8)
  ht = tmax/real(Nt-1,8)
  R1 = alpha*ht/(hx**2)

  A%bw = 1

  inz = 0
  do row=ibeg,iend

! Set the diagonal entry in row
	
     v   = 1.0_8 / (1.0_8 + 4.0_8*R1)
     inz = inz + 1

     A%aa(inz) = v
     A%ii(row) = inz
     A%jj(inz) = row

  end do
 
! Set A%n (the dimension of A), A%nnz (the number of nonzero entries) and
! A%ii(A%n+1) which is needed to address the nonzero entries in the last
! row (in order to know where the arrays A%aa and A%jj end)

  A%bw   = m-1
  A%nnz  = inz
  A%ii(A%iend+1) = A%nnz + 1
  
end subroutine BTCSmatrix_PRE

