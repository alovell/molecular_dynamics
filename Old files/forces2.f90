      module forces

      contains

      subroutine calc_forces(N,pos,forces)

      !implicit none
      integer, intent(in)  :: N
      real(8), intent(in)  :: pos(N,3)
      real(8), intent(inout)  :: forces(N,3)
      real(8) :: dx(N),dy(N),dz(N),rc,sig,ep,dr(N,N),mass
      integer :: i,j
      real(8)  :: dt !time step
      
      rc=1 !define cut-off radius, place holder for now
      mass=18 !place holder mass for argon


      !calculating the radial distance between each pair of particles
      do j= 1,N
        do i=1,N
          dx(i) = pos(i,1) - pos(j,1)
          dy(i) = pos(i,2) - pos(j,2)
          dz(i) = pos(i,3) - pos(j,3)
          dr(i,j) = sqrt(dx(i)**2 + dy(i)**2 + dz(i)**2)
          print *, dr(i,j)
          if (dr(i,j) > rc) then
            dr(i,j) = 0
          end if
        enddo
      enddo
      
      !Leonnard-Jones Potential
      !VLJ = 4*ep*[(sig/dr)^12 - (sig/dr)^6
      !(rm,ep) is the minimum of the potential well 
      !with rm = 2^(1/6)*sig
      
      ep=1 !place holder for epsilon
      sig=2/(2**(0.16666667)) !2 is a place holder for rm
      !FLJ from -grad VLJ
      !Fi = 24*ep[2*sig**12/r**14 - sig**6/r**8]*ri
      
      !calculate the force on each particle
      do j= 1,N
        do i=1,N
	  !if distance is greater than critical distance
	  !or the particles are the same force is zero
	  if (dr(i,j)==0) then 
	    forces(i,1)=0
	    forces(i,2)=0
	    forces(i,3)=0
	  else
	    !else caluclate the force from VLJ
            forces(i,1) = 24*ep*dx(i)*(2*sig**12/(r(i,j)**14) 
     & - sig**6/(r(i,j)**8))
	    forces(i,2) = 24*ep*dy(i)*(2*sig**12/(r(i,j)**14) 
     & - sig**6/(r(i,j)**8))
	    forces(i,3) = 24*ep*dz(i)*(2*sig**12/(r(i,j)**14) 
     & - sig**6/(r(i,j)**8))
	  end if
        enddo
      enddo
      
      !print forces to see if there are zeros where we expect
      do j=1,N
        do i=1,N
          print *, forces(i,j)
	enddo 
      enddo
      
      !this section is commented out - need to read in mom vector
      !also need to read out new pos and mom vectors
      !calculates the updated momenta
      !p(k+1) = p(k) + F(i)*t
      
      !here is the actual code for the momenta
      !do i=1,N
        !mom(i,1) = mom(i,1)+forces(i,1)*dt
        !mom(i,2) = mom(i,2)+forces(i,2)*dt
        !mom(i,3) = mom(i,3)+forces(i,3)*dt
      !enddo 
      
      !calculates the updated positions
      !do i=1,N
        !pos(i,1) = pos(i,1) + mom(i,1)*dt/mass + (0.5d0*forces(i,1)*dt**2)/mass
	!pos(i,2) = pos(i,2) + mom(i,2)*dt/mass + (0.5d0*forces(i,2)*dt**2)/mass
	!pos(i,1) = pos(i,3) + mom(i,3)*dt/mass + (0.5d0*forces(i,3)*dt**2)/mass
      !enddo 

      end subroutine calc_forces

      end module forces
