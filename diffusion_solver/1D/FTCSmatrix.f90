subroutine FTCSmatrix(A,alpha,tmax,Nx,Nt)

! Subroutine to form matrix arriving from 
! FTCS for the heat equation
! Non CRS version

real(kind=8), intent(out), dimension(Nx,Nx) :: A
integer, intent(in) :: Nx, Nt
real(kind=8), intent(in) :: alpha, tmax
real (kind=8) :: xx,hx,xt,ht,R1,R2

hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)
R1 = alpha*ht/(hx**2)
R2 = 1.0 - 2.0*R1

A = 0.0
do i = 1,Nx		
	A(i,i) = R2
	! Off diag entries
	if (i>2) then  
     	A(i,i-1) = R1
    end if
    if (i<Nx) then 
     	A(i,i+1) = R1
    end if
    A(1,1) = 1.0_8
	A(1,2) = 0.0_8
	A(Nx,Nx) = 1.0_8
	A(Nx,Nx-1) = 0.0_8
end do 

end subroutine FTCSmatrix
