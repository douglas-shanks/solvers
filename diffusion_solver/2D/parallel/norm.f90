subroutine norm(error,U1,U2)

! Computes the norm (parallel)

 use header
 include "mpif.h"
 
 type(Vector), intent(in) :: U1, U2
 real(kind=8), intent(out) :: error
 real(kind=8) :: myerror
 error = 0
 
 do i=U2%ibeg,U2%iend
 
 	myerror = myerror + SQRT((U1%xx(i)- U2%xx(i))**2) / SQRT(real((m-1)**2,8))
 
 enddo
call MPI_Allreduce(myerror,error,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                     MPI_COMM_WORLD,ierr)
end subroutine norm
