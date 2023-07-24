!-----------------------------------------------------------------------
      subroutine PUSH1ZF(part,dt,ek,idimp,nop,nx)
! for 1d code, this subroutine updates particle co-ordinate for
! particles with fixed velocities and periodic boundary conditions
! 6 flops/particle, 2 loads, 1 store
! input: all, output: part, ek
! equations used are:
! x(t+dt) = x(t) + vx(t+dt/2)*dt
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! dt = time interval between successive calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2)
! idimp = size of phase space = 2
! nop = number of particles
! nx = system length in x direction
      implicit none
      integer nop, idimp, nx
      real part, dt, ek
      dimension part(idimp,nop)
! local data
      integer i
      real anx, dx, vx
      double precision sum1
      anx = real(nx)
      sum1 = 0.0d0
      do 10 i = 1, nop
      dx = part(1,i)
      vx = part(2,i)
! average kinetic energy
      sum1 = sum1 + vx*vx
! new position
      dx = dx + vx*dt
! periodic boundary conditions
      if (dx.lt.0.0) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
! set new position
      part(1,i) = dx
   10 continue
! normalize kinetic energy
      ek = ek + .125*real(sum1)
      return
      end
