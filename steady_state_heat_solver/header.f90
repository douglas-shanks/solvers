module header

  implicit none

  type Matrix
     real (kind=8), allocatable, dimension(:) :: aa
     integer, allocatable, dimension(:) :: ii
     integer, allocatable, dimension(:) :: jj
     integer :: n
     integer :: nnz
     integer :: ibeg
     integer :: iend
     integer :: bw
  end type Matrix

  type Vector
     real (kind=8), allocatable, dimension(:) :: xx
     integer :: n
     integer :: ibeg
     integer :: iend
  end type Vector

end module 
