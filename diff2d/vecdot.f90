!
!  Description: Subroutine to calculate the dot product u.v
!	Parallel version
!
! -----------------------------------------------------------------------

subroutine Vec_Dot(x,y,sigma)

  use header
  include "mpif.h"

  type(Vector), intent(in) :: x,y
  real (kind=8), intent(out)  :: sigma
  real(kind=8) :: ddot
  real(kind=8) :: mysigma

! Use a level 1 BLAS call to caluclate the dot product

  external ddot
  mysigma = ddot(x%iend-x%ibeg+1,x%xx(x%ibeg),1,y%xx(y%ibeg),1)

! Combine product from all processes and distribute the result to all processes
 
  call MPI_ALLREDUCE(mysigma,sigma,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                    MPI_COMM_WORLD,ierr)

end subroutine Vec_Dot
