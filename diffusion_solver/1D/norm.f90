subroutine norm(error,U1,U2,Nt,Nx)

! This actually computes the error not the norm

integer, intent(in) :: Nx, Nt 
real(kind=8), intent(in), dimension(Nx) :: U1
real(kind=8), intent(in), dimension(Nt,Nx) :: U2
real(kind=8), intent(out) :: error

error = 0
do i=1,Nx
	error = error + SQRT((U1(i)- U2(Nt,i))**2)/SQRT(real(Nx,8) 	  )
enddo

end subroutine norm
