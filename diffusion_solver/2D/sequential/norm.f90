subroutine norm(error,U1,U2,m)

! Computes the norm

 use header
 
 type(Vector) :: U1, U2
 integer, intent(in) :: m
 real(kind=8), intent(out) :: error

 error = 0
 
 do i=1,(m-1)**2
 
 	error = error + SQRT((U1%xx(i)- U2%xx(i))**2) / SQRT(real((m-1)**2,8))
 
 enddo

end subroutine norm
