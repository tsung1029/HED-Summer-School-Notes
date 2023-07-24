!-----------------------------------------------------------------------
      program f2pydemo2
      use init1
      use rpush1
      implicit none
!
! Inputs:
! indx = exponent which determines grid points in x direction:
! nx = 2**indx.
      integer, parameter :: indx =   9
! npx = number of electrons distributed in x direction.
      integer, parameter :: npx = 18432
! nloop = number of time steps in the simulation
      integer, parameter :: nloop = 100
! dt = time interval between successive calculations.
      real, parameter :: dt = 0.1
! vtx = thermal velocity of electrons in x direction
! vx0 = drift velocity of electrons in x direction.
      real, parameter :: vtx = 1.0, vx0 = 0.0
! ci = reciprocal of velocity of light.
      real, parameter :: ci = 0.1
! idimp = number of particle coordinates = 2
      integer, parameter :: idimp = 2
! wke = particle kinetic energy
      real :: wke = 0.0
!
! declare scalars for standard code
      integer :: np, nx, n
!
! declare arrays for standard code:
! part = particle array
      real, dimension(:,:), allocatable :: part
!
! rpart = helper object for relativistic particles
      type (part1d) :: rpart
!
! initialize scalars for standard code
! nx = number of grid points in x direction
! np = total number of particles in simulation
      nx = 2**indx; np = npx; 
!
! construct helper object for electrons: update rpart
      call new_part1d(rpart,0.0,0.0,idimp,np)
!
! allocate data for standard code
      allocate(part(idimp,np))
!
! initialize uniform plasma and maxwellian velocity: updates part
      call DISTR1(part,vtx,vx0,npx,nx,0)
!
      write (*,*) 'initial:part(:,1)=',part(:,1)
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop 
!
! push free-streaming relativistic particles: updates part, wke
         wke = 0.0
         call rpush1zf(rpart,part,dt,ci,wke,nx)
!
      enddo
!
! * * * end main iteration loop * * *
!
      write (*,*) 'final:part(:,1)=',part(:,1)
!
      stop
      end program
