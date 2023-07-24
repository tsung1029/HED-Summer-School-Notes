!-----------------------------------------------------------------------
!
      module rpush1
      implicit none
!
      type part1d
! qm = charge of particle, in units of e
! qbm = particle charge/mass ratio
! idimp = size of phase space = 2
! nop = number of particles in array
         real :: qm, qbm
! idimp = size of phase space = 2
! nop = number of particles in array
         integer :: idimp, nop
      end type
!
      contains
!
!-----------------------------------------------------------------------
      subroutine new_part1d(this,qm,qbm,idimp,nop)
! this subroutine creates a descriptor for 1d particles
! this = part1d descriptor of particle data
      implicit none
      type (part1d) :: this
      real, intent(in) :: qm, qbm
      integer, intent(in) :: idimp, nop
! create descriptor
      this%qm = qm; this%qbm = qbm
      this%idimp = idimp; this%nop = nop
      end subroutine
!
!-----------------------------------------------------------------------
      subroutine rpush1zf(this,part,dt,ci,ek,nx)
! for 1d code, this subroutine updates particle co-ordinate for
! particles with fixed velocities and periodic boundary conditions
! 9 flops, 2 divides, 1 sqrt/particle, 2 loads, 1 store
! input: all, output: part, ek
! equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! dtg = dt/sqrt(1.+(px(t+dt/2)*px(t+dt/2)*ci*ci)
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = sum((px(t-dt/2))**2)/(1. + gamma)
! where gamma = sqrt(1.+((px(t+dt/2)*ci)**2)
! idimp = size of phase space = 2
! nop = number of particles
! nx = system length in x direction
      implicit none
      type (part1d), intent(in) :: this
      integer, intent(in) :: nx
      real, intent(in) :: ci, dt
      real, intent(inout) :: ek
      real, dimension(:,:), intent(inout) :: part
! local data
      integer :: i
      real :: ci2, anx, dx, px, p2, gam, dtg
      double precision sum1
      anx = real(nx)
      ci2 = ci*ci
      sum1 = 0.0d0
      do i = 1, this%nop
         dx = part(1,i)
         px = part(2,i)
! average kinetic energy
         p2 = px*px
         gam = sqrt(1.0 + p2*ci2)
         sum1 = sum1 + p2/(1.0 + gam)
! update inverse gamma
         dtg = dt/gam
! new position
         dx = dx + px*dtg
! periodic boundary conditions
         if (dx.lt.0.0) dx = dx + anx
         if (dx.ge.anx) dx = dx - anx
! set new position
         part(1,i) = dx
      enddo
! normalize kinetic energy
      ek = ek + real(sum1)
      end subroutine
!
      end module
