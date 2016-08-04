!==================================================================
!
!   Function evaluating the material parameter at the point (x,y) in 
! 	the PDE
! 		-\grad a(x,y) \grad U + sigma U = g,  in D
!									  U = 0,  on D
!   (At the moment alpha(x,y) = 1.0 everywhere, i.e. the Laplacian)
!
!==================================================================

function alpha(x,y) result (val)

use header

real(kind=8), intent(in) :: x,y
real(kind=8) :: val

!val = 1.0d0
!val = SIN(x)*COS(y)

!diffusion coefficient different in each half (big improvement from PCG)

!if((0.25_8<=x .AND. x<=0.75_8).AND.(0.25_8<=y .AND. y<=0.75_8)) then
!   val=1.0d-6
!else
!   val = 1.0d0
!end if

!diffusion coefficient different in each quadrant (big improvement from PCG)

if((0.5_8<=x .AND. x<=0.75_8).AND.(0.0_8<=y .AND. y<=0.25_8))then
   val=1.0d2
elseif((0.25_8<=x .AND. x<=0.5_8).AND.(0.25_8<=y .AND. y<=0.5_8))then
   val=1.0d2
elseif((0.5_8<=x .AND. x<=0.75_8).AND.(0.5_8<=y .AND. y<=0.75_8))then
   val=1.0d2
elseif((0.25_8<=x .AND. x<=0.5_8).AND.(0.75_8<=y .AND. y<=1.0_8))then
   val=1.0d2
else
   val = 1.0d0
end if

end function alpha
