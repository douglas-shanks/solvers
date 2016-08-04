
subroutine greetings(myid,numprocs)

 include "mpif.h"

 integer, intent(in) :: myid, numprocs

 integer :: len
 parameter (len = MPI_MAX_PROCESSOR_NAME+1)
 character*(len) name,chbuf
 integer :: stat(MPI_STATUS_SIZE),ierr,i,dest,reslen

 call MPI_GET_PROCESSOR_NAME(name,reslen,ierr)

 if (myid == 0) then
  
    print*,' There are ',numprocs,'processes running'
    print*,' Master process ',myid,' runs on  ',name

    do i=1,numprocs-1
       call MPI_RECV(chbuf,len,MPI_CHARACTER,MPI_ANY_SOURCE, &           
            MPI_ANY_TAG,MPI_COMM_WORLD, stat,ierr)
       print*,' Slave process ',stat(MPI_SOURCE), ' runs on  ',chbuf
    end do

 else

    dest = 0
    call MPI_SEND(name,len,MPI_CHARACTER,dest,myid, &
         MPI_COMM_WORLD,ierr)
 endif

end subroutine greetings




