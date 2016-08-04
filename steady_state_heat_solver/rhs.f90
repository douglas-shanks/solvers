!==========================================================================
!
!  Subroutine to set the RHS for the 2-D, five-point stencil 
!  of a general 2nd-order elliptic PDE in compressed row storage format 
!  (parallel version).
!  The discretied PDE
! 		-\grad a(x,y) \grad U + sigma U = g,  in D
!									  U = 0,  on D
! g is set such that u = y when i = 1
!					   = 0 when j = 1
!					   = y when i = m-1
!					   = 1 when i = m-1
!==========================================================================

subroutine RHS(b,m)

  use header

  type(Vector), intent(inout) :: b
  integer, intent(in) :: m

  integer :: row,i,j,c
  real(kind=8) :: xm,xp,ym,yp

  interface
    function alpha(x,y) result (val)

    use header

    real(kind=8), intent(in) :: x, y
    real(kind=8) :: val

    end function alpha
  end interface

  c = m*m

  b%xx(b%ibeg:b%iend) = 0.0d0

  do row=b%ibeg,b%iend

! To set the entries in RHS we calculate the indices (i,j) in the Cartesian 
! numbering of the unknowns and then set the entries accordingly. There are 
! only entries in the RHS for unknowns that are "connected" to the boundary
  
     j = (row-1)/(m-1) + 1
     i = row - (j-1)*(m-1)  

     xm = (i-0.5d0)/m
     xp = (i+0.5d0)/m
     ym = (j-0.5d0)/m
     yp = (j+0.5d0)/m


! If i = 1 then the unknown is next to the left boundary where u=y

     if (i == 1) then
        b%xx(row) = b%xx(row) + 0.5d0*c*(alpha(xm,ym) + alpha(xm,yp)) * dble(j)/m
     endif

! If i = m-1 then the unknown is next to the right boundary where u=y also

     if (i == m-1) then
        b%xx(row) = b%xx(row) + 0.5d0*c*(alpha(xp,ym) + alpha(xp,yp)) * dble(j)/m
     endif

! If j = m-1 then the unknown is next to the top boundary where u=1

     if (j == m-1) then
        b%xx(row) = b%xx(row) + 0.5d0*c*(alpha(xm,yp) + alpha(xp,yp))
     endif

  end do
 
end subroutine RHS
