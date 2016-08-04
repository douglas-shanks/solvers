subroutine CN_rhs(phi,b,Nx,Nt,tmax)

! finds the right-hand side vector in the approximation of Poisson's 
! equation on the unit square

integer, intent(in) :: Nx,Nt
real(kind=8), intent(in) :: tmax
real(kind=8), intent(in), dimension(Nx) :: b
real(kind=8), intent(out), dimension(Nx) :: phi
integer :: i,j
real (kind=8) :: xx,hx,xt,ht,R

hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)

R = alpha*ht/(hx*hx)

phi(1) = b(1)
do  i = 2,Nx-1
    phi(i) =  R*b(i+1) + 2.0_8*(1.0_8 - R)*b(i) + R*b(i-1)
end do
phi(Nx) = b(Nx)

end subroutine CN_rhs
