! Creates a FD solution on a grid that can be written to file for
! post processing
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

subroutine FDsolution(uu,u,m)

  use header

  implicit none

  integer, intent(in)      :: m
  type(Vector), intent(in) :: u

  real(kind=8), intent(inout), dimension(m-1,m-1):: uu
  integer :: ierr,j,mdelta


  do j=1,m-1
     uu(j,:) = u%xx((j-1)*(m-1)+1:j*(m-1))
  end do


end subroutine FDsolution
