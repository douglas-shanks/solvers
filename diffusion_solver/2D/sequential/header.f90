module header

  implicit none

  type Matrix
     real (kind=8), allocatable, dimension(:) :: aa
     integer, allocatable, dimension(:) :: ii
     integer, allocatable, dimension(:) :: jj
     integer :: n
     integer :: nnz
  end type Matrix

  type Vector
     real (kind=8), allocatable, dimension(:) :: xx
     integer :: n
  end type Vector

end module 
