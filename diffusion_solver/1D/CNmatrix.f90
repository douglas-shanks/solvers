subroutine CNmatrix(A,alpha,tmax,Nx,Nt)

! Subroutine to form matrix arriving from 
! BTCS for the heat equation
! Non CRS version

real(kind=8), intent(out), dimension(Nx,Nx) :: A
integer, intent(in) :: Nx, Nt
real(kind=8), intent(in) :: alpha, tmax
real (kind=8) :: xx,hx,xt,ht,R

hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)

R = alpha*ht/(hx*hx)

A = 0.0_8
do i = 1,Nx		
	A(i,i) = 2.0_8 + 2.0_8*R
	! Off diag entries
	if (i>1) then  
     	A(i,i-1) = -R
    end if
    if (i<Nx) then 
     	A(i,i+1) = -R
    end if
end do 

A(1,1) = 1.0_8
A(1,2) = 0.0_8
A(Nx,Nx) = 1.0_8
A(Nx,Nx-1) = 0.0_8

end subroutine CNmatrix
