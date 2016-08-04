!
!  Description: Subroutine to calculate the dot product u.v
!
! -----------------------------------------------------------------------

subroutine Vec_Dot(u,v,sigma)

  use header

  type(Vector), intent(in) :: u, v
  real(kind=8), intent(out) :: sigma
  real(kind=8) :: ddot

  external ddot  ! Uses BLAS Level 1

  sigma = ddot(u%n,u%xx,1,v%xx,1)

end subroutine Vec_Dot
