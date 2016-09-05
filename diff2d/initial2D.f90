subroutine initial2D(U,alpha,tmax,m)

! subroutine initial2D
! computes the intial condition for the heat equation with Dirichlet
! boundary 
! For example the initial condition u_{0}(x,y) = sin(\pi x)sin(\pi y)

  use header
  
  ! Needs to be inout here as memory for U is read in
  type(Vector), intent(inout) :: U
  integer, intent(in) :: m
  real(kind=8), intent(in) :: alpha,tmax


  integer :: i,j
  real (kind=8) :: xx1,yy1,hx,hy,lambda
  real (kind=8) :: PI = 4.0_8*ATAN(1.0_8)
   
  hx = 1.0_8/real(m,8)
  ht = tmax/real(Nt-1,8)

! Smooth sine centered

!  do i = 1,m-1
!	 do j = 1,m-1
	 
!		 xx1 = (i)*hx
!		 yy1 = (j)*hx
!		 lambda = SQRT(alpha)*PI*SQRT(2.0_8)
!		 U%xx(i+(j-1)*(m-1)) = SIN(PI*xx1)*SIN(PI*yy1)*EXP(-lambda**2.0*tmax)
		 
!	 end do
!  end do

! Gaussian pulse
  do i = 1,m-1
	 do j = 1,m-1
	 
		 xx1 = (i)*hx
		 yy1 = (j)*hx
         if ( i == m/2 .and. j == m/2 ) then
		 	U%xx(i+(j-1)*(m-1)) = m/4
		 else
		 	U%xx(i+(j-1)*(m-1)) = 0.0_8
		 endif
	 end do
  end do

end subroutine initial2D
