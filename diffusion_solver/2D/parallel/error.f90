! -----------------------------------------------------------------------
subroutine Error(u,v,norm)
! -----------------------------------------------------------------------
!
!  Description: Subroutine to calculate || u - v ||_2 
!
! -----------------------------------------------------------------------

  use header
  include "mpif.h"

  type(Vector), intent(in) :: u, v
  real(kind=8), intent(out) :: norm

  real(kind=8) :: mynrm
  integer :: i

  mynrm = 0.0_8
  
  ! need this ?call sparsegather(v,u%n)
  
  do i=u%ibeg,u%iend 
     mynrm = mynrm + (u%xx(i) - v%xx(i)) * (u%xx(i) - v%xx(i))
  end do

  ! Taking this out solved the problem .... why? Because its just computing hte error on one processor. It's not bringing al the results back.
  !call MPI_Allreduce(mynrm,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
  !                   MPI_COMM_WORLD,ierr)
  norm = sqrt(mynrm/u%n)

end subroutine Error
