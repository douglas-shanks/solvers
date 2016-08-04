program test1D

!---------------------------------------------------------------------
!
!	A test program to check that 1D FTCS, BTCS and CN converges
!
!---------------------------------------------------------------------

implicit none

integer:: Nx, Nt
integer i,j
real(kind=8), dimension(:), allocatable :: Uexact,xftcs,xbtcs,yftcs,xcn,ycn
real(kind=8), dimension(:,:), allocatable :: Uftcs,Ubtcs,Ucn,Aftcs,Abtcs,Acn
real(kind=8) :: xx,hx,xt,ht,R1,R2
real(kind=8):: alpha,tmax
real(kind=8):: errorftcs,errorbtcs,errorcn,alpha1,beta1

! extra storage needed by dgesv, dpbsv:
integer, dimension(:), allocatable :: ipiv1, ipiv2
integer :: info1, info2

! Define Pi for use
real (kind=8) :: PI = 4.0_8*ATAN(1.0_8)

! Constants needed for dgmev
alpha1=1.d0; beta1=0.d0

WRITE(*,*)
WRITE(*,*) '------------------------------------------------------------------'
WRITE(*,*) '----This is a test program----------------------------------------'
WRITE(*,*) '----We solve the 1D heat equation---------------------------------'
WRITE(*,*) '----Comparing the exact solution to that with the FTCS method-----'
WRITE(*,*) '----BTCS method and Crank Nicholson method -----------------------'
WRITE(*,*) '------------------------------------------------------------------'
WRITE(*,*)
WRITE(*,*) '---------------------What is alpha -------------------------------'
READ(*,*) alpha
WRITE(*,*) '---------------------What is tmax --------------------------------'
READ(*,*) tmax
WRITE(*,*)
WRITE(*,*) '---------------------What is Nx-----------------------------------'
READ(*,*) Nx
WRITE(*,*) '---------------------What is Nt-----------------------------------'
READ(*,*) Nt
WRITE(*,*)
if (Nx<1 .or. Nt<1) then
	print *, "Nx and Nt must be positive "
    stop
endif

! Mesh spacing for time and space
hx = 1.0_8/real(Nx-1,8)
ht = tmax/real(Nt-1,8)
R1 = alpha*ht/(hx**2)

print*, '-------------------------value of r ------------------------------'
print*,               R1

! Allocate storage
allocate(yftcs(Nx))
allocate(xftcs(Nx))
allocate(xbtcs(Nx))
allocate(xcn(Nx))
allocate(ycn(Nx))
allocate(Uexact(Nx))
allocate(Uftcs(Nt,Nx))
allocate(Ubtcs(Nt,Nx))
allocate(Ucn(Nt,Nx))
allocate(Aftcs(Nx,Nx))
allocate(Abtcs(Nx,Nx))
allocate(Acn(Nx,Nx))
allocate(ipiv1(Nx))
allocate(ipiv2(Nx))

!-----------------------------------------------------------
!					EXACT SOLUTION
!-----------------------------------------------------------

call exactsoln(Uexact,alpha,tmax,Nx)

!-----------------------------------------------------------
!						FTCS
!-----------------------------------------------------------

!call FTCS(Ua,alpha,tmax,Nx,Nt)
! Now do the matrix version and check solution agrees
call FTCSmatrix(Aftcs,alpha,tmax,Nx,Nt)

! Initial solution, u_0
do i = 1,Nx
	xx = (i-1)*hx
	Uftcs(1,i) = SIN(PI*real(xx))
end do

! Then it should just be a case of U^{j} = A*U{j-1}
do j=2,Nt

	xftcs = Uftcs(j-1,:)
  	call dgemv('Nx',Nx,Nx,alpha1,Aftcs,Nx,xftcs,1,beta1,yftcs,1)
	Uftcs(j,:) = yftcs
	
end do

! enforce boundary conditions
Uftcs(:,1) = 0.0
Uftcs(:,Nx)= 0.0

! Compute the error
call norm(errorftcs,Uexact,Uftcs,Nt,Nx)

write(*,*)
write(*,*) '------------------------------------------------------------------'
print*,    '----------------------- FTCS Solution ----------------------------'
write(*,*) '------------------------------------------------------------------'
write(*,*)
write(*,*) '-------------------- Error ---------------------------------------'
write(*,*)
print*, errorftcs	 
write(*,*) '------------------------------------------------------------------'
write(*,*)

!-----------------------------------------------------------
!						BTCS
!-----------------------------------------------------------

! Initial solution, u_0
do i = 1,Nx
	xx = (i-1)*hx
	Ubtcs(1,i) = SIN(PI*real(xx))
end do

! Then it should just be a case of A U^{j} = (1/dt)U{j-1}
do j=2,Nt
	! Make BTCS matrix and RHS (This has to be done at each
	! iterate as apparetnly dgesv changes the matrix!
	call BTCSmatrix(Abtcs,alpha,tmax,Nx,Nt)
	
	do  i = 1,Nx
    	xbtcs(i) = Ubtcs(j-1,i)
	end do
	!xbtcs = Ubtcs(j-1,:)
	!call BTCS_rhs(xbtcs,Nx,Nt,tmax)

	! now solve, and produce output information
  	call dgesv(Nx,1,Abtcs,Nx,ipiv1,xbtcs,Nx,info1)    
	Ubtcs(j,:) = xbtcs	
	
end do

! enforce boundary conditions
Ubtcs(:,1) = 0.0
Ubtcs(:,Nx)= 0.0

! Compute the error
call norm(errorbtcs,Uexact,Ubtcs,Nt,Nx)

write(*,*) '------------------------------------------------------------------'
print*,    '------------------------ BTCS Solution ---------------------------'
write(*,*) '------------------------------------------------------------------'
write(*,*)
write(*,*) '--------------------- Error --------------------------------------'
print*,
print*, errorbtcs	 
write(*,*) '------------------------------------------------------------------'
write(*,*)

!-----------------------------------------------------------
!				Crank Nicholson
!-----------------------------------------------------------

! Initial solution, u_0
do i = 1,Nx
	xx = (i-1)*hx
	Ucn(1,i) = SIN(PI*real(xx))
end do

do j=2,Nt

	! Make CN matrix and RHS (This has to be done at each
	! iterate as apparetnly dgesv changes the matrix!
	call CNmatrix(Acn,alpha,tmax,Nx,Nt)
	
	! make the rhs
	xcn(1) = Ucn(j-1,1)
	do  i = 2,Nx-1
    	xcn(i) =  R1*Ucn(j-1,i+1) + 2.0_8*(1.0_8 - R1)*Ucn(j-1,i) + R1*Ucn(j-1,i-1)
	end do
	xcn(Nx) = Ucn(j-1,Nx)
	!xcn = Ucn(j-1,:)
	!call CN_rhs(xcn,Ucn(j-1,Nx),Nx,Nt,tmax)

	! now solve, and produce output information
  	call dgesv(Nx,1,Acn,Nx,ipiv2,xcn,Nx,info2)  
	Ucn(j,:) = xcn
 
end do

! enforce boundary conditions
!Ucn(:,1) = 0.0
!Ucn(:,Nx)= 0.0	

! Compute the error
call norm(errorcn,Uexact,Ucn,Nt,Nx)

write(*,*) '------------------------------------------------------------------'
print*,     '------------------------ CN Solution -----------------------------'
write(*,*) '------------------------------------------------------------------'
write(*,*)
write(*,*) '--------------------- Error --------------------------------------'
print*,
print*, errorcn	 
write(*,*) '------------------------------------------------------------------'
write(*,*)
!print*, Uexact
!print*,Uftcs(Nt,:)
!print*,Ubtcs(Nt,:)
!print*,Ucn(Nt,:)

!---------------------------------------------------------------
! deallocate storage
!---------------------------------------------------------------
deallocate(Uexact,Uftcs,Ubtcs,Ucn)
deallocate(Aftcs,Abtcs,Acn)
deallocate(xftcs,yftcs,ycn)
deallocate(xbtcs,xcn)
deallocate(ipiv1,ipiv2)

WRITE(*,*) '------------------------------------------------------------------' 
WRITE(*,*) '--------------- Program finished ---------------------------------'
WRITE(*,*) '------------------------------------------------------------------'

end program test1D
