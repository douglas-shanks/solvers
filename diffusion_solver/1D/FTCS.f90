subroutine FTCS(U,alpha,tmax,Nx,Nt)

! subroutine FTCS
! computes the exact solution of the heat equation with Dirichlet
! boundary and initial condition u_{0}(x) = sin(\p ix)
! using FTCS.
! found error, always start i from 1!

integer, intent(in) :: Nx, Nt
real(kind=8), intent(in) :: alpha,tmax
real(kind=8), intent(out), dimension(Nt,Nx) :: U

integer :: i,j
real (kind=8) :: xx,hx,xt,ht,R1,R2
real (kind=8) :: PI = 4.0_8*ATAN(1.0_8)
 
hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)

R1 = alpha*ht/(hx**2)
R2 = 1.0 - 2.0*R1

do i = 1,Nx
	xx = (i-1)*hx
	U(1,i) = SIN(PI*real(xx))
end do

do j=2,Nt
	do i=2,Nx-1
		U(j,i) = R1*U(j-1,i-1) + R2*U(j-1,i) + R1*U(j-1,i+1)
	end do	
end do

end subroutine FTCS
