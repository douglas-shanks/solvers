! -----------------------------------------------------------------------
subroutine Error(u,v,norm)
! -----------------------------------------------------------------------
!
!  Description: Subroutine to calculate || u - v ||_2 
!
! -----------------------------------------------------------------------

  use header

  type(Vector) :: u, v
  real(kind=8), intent(out) :: norm

  integer :: i

  norm = 0.0_8
  
  do i=1,u%n
     norm = norm + (u%xx(i) - v%xx(i)) * (u%xx(i) - v%xx(i))
  end do

  norm = sqrt(norm/u%n)

end subroutine Error
