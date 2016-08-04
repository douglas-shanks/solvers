
subroutine sparsegather(x,BW)
 
 use header
 include "mpif.h"

 type(Vector), intent(inout) :: x
 integer, intent(in) :: BW

 integer :: myid,numprocs,stat(MPI_STATUS_SIZE),ierr
 
 call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
 call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)

 if(myid<numprocs-1) then
    call MPI_Send(x%xx(x%iend-BW+1:x%iend),BW,&
         MPI_DOUBLE_PRECISION,myid+1,myid+1, MPI_COMM_WORLD,ierr)
 endif

 if(myid>0) then
    call MPI_RECV(x%xx(x%ibeg-BW:x%ibeg-1),BW,&
         MPI_DOUBLE_PRECISION,myid-1, &           
         MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
 endif

 if(myid>0) then
    call MPI_Send(x%xx(x%ibeg:x%ibeg+BW-1),BW,&
         MPI_DOUBLE_PRECISION,myid-1,myid-1, MPI_COMM_WORLD,ierr)
 endif

 if(myid<numprocs-1) then
    call MPI_RECV(x%xx(x%iend+1:x%iend+BW),BW,&
         MPI_DOUBLE_PRECISION,myid+1, &           
         MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
 endif


end subroutine sparsegather




