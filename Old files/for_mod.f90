!module forces

!contains

subroutine calc_forces(N,pos,forces)

!implicit none
integer, intent(in)  :: N
real(8), intent(in)  :: pos(N,3)
real(8), intent(inout)  :: forces(N,3)
real(8) :: dx,dy,dz,rc,sig,ep,dr(N,N)
integer :: i,j


!calculating the radial distance between each pair of particles
do j= 1,N
   do i=1,N
      dx = pos(i,1) - pos(j,1)
      dy = pos(i,2) - pos(j,2)
      dz = pos(i,3) - pos(j,3)
      dr(i,j) = sqrt(dx**2 + dy**2 + dz**2)
      print *, dr(i,j)
      if (dr(i,j) < rc) then
        dr(i,j) = 0
      end if
   enddo
enddo

end subroutine calc_forces

!end module forces
