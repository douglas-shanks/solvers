subroutine BTCS_rhs(b,Nx,Nt,tmax)

! finds the right-hand side vector in the approximation of Poisson's 
! equation on the unit square

integer, intent(in) :: Nx,Nt
real(kind=8), intent(in) :: tmax
real(kind=8), intent(out), dimension(Nx) :: b
integer :: i,j
real (kind=8) :: xx,xt,hx,ht

hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)

do  i = 1,Nx
    b(i) = (1.0/ht)*b(i)
end do
 
end subroutine BTCS_rhs
