 
      subroutine calc_forces(N,pos,forces,l)

      !implicit none
      integer, intent(in)  :: N
      real(8), intent(in)  :: pos(N,3),l
      real(8), intent(inout)  :: forces(N,3)
!f2py intent (in,out) :: forces
      real(8), parameter :: rc =3.2
      real(8) :: dr_vec(3),partialforce(3),sig,ep,dr2,F
      integer :: i,j
      !real(8)  :: dt !time step
      
      !rc=3.2 !define cut-off radius, place holder for now
      !mass=1._8 !place holder mass for argon

      !Leonnard-Jones Potential
      !VLJ = 4*ep*[(sig/dr)^12 - (sig/dr)^6
      !(rm,ep) is the minimum of the potential well 
      !with rm = 2^(1/6)*sig
      
      ep=1 !place holder for epsilon
      sig=1 !recommended by jos
      !FLJ from -grad VLJ
      !Fi = 24*ep[2*sig**12/r**14 - sig**6/r**8]*ri
      

      forces(:,:) = 0.0 
      print '(3F12.6)', pos(1:3,:)
      print *, l
      !calculating the radial distance between each pair of particles
      !Afterwards, do forces
      do j = 1, N
        do i = j+1, N
          dr_vec = pos(i,:) - pos(j,:)
          dr_vec = dr_vec - nint(dr_vec/l)*l
!          if dr_vec(:)
          dr2 = dot_product(dr_vec,dr_vec)
          if (dr2 < rc*rc) then
            dr2 = 1/dr2
            f= 2*dr2**7 - dr2**4
!            print *, dr2
            if (f>1.D3) then 
              print *, i, j, dr2
              print '(3F12.6)', pos(i,:), pos(j,:)
              pause
            end if
            partialforce(:) = 24*dr_vec*f
            forces(j,:) = forces(j,:) + partialforce(:)
            forces(i,:) = forces(i,:) - partialforce(:)
          end if
        enddo
      enddo
 !     forces(:,:) = 0.0 
      end subroutine calc_forces
      
     
      !print forces to see if there are zeros where we expect
      !do j=1,N
      !  do i=1,N
      !    print *, forces(i,j)
      !	enddo 
      !enddo
      
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

      


