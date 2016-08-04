subroutine exactsoln(U,alpha,tmax,Nx)

! subroutine exactsoln
! computes the exact solution of the heat equation with Dirichlet
! boundary and initial condition u_{0}(x) = sin(\p ix)
! At the final time Tf
! found error, always start i from 1!

integer, intent(in) :: Nx
real(kind=8), intent(in) :: alpha,tmax
real(kind=8), intent(out), dimension(Nx) :: U

integer :: i
real (kind=8) :: xx,hx
real (kind=8) :: PI = 4.0_8*ATAN(1.0_8)
 
hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)

do i = 1,Nx+1
	xx = (i-1)*hx
	U(i) = SIN(PI*real(xx))*EXP(-tmax*alpha*PI**2)
end do



end subroutine exactsoln
