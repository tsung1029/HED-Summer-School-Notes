!-----------------------------------------------------------------------
      module init1
!
      contains  
!
!-----------------------------------------------------------------------
      subroutine DISTR1(part,vtx,vdx,npx,nx,ipbc)
! for 1d code, this subroutine calculates initial particle co-ordinate
! and velocity, with uniform density and maxwellian velocity with drift
! part(1,n) = position x of particle n
! part(2,n) = velocity vx of particle n
! vtx = thermal velocity of particles in x direction
! vdx = drift velocity of particles x direction
! npx = number of particles distributed in x direction
! nx = system length in x direction
! ipbc = particle boundary condition = (0,1,2) =
! (none,2d periodic,2d reflecting)
! ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer, intent(in) :: npx, nx, ipbc
      real, intent(in) :: vtx, vdx
      real, dimension(:,:), intent(inout) :: part
! local data
      integer j, idimp, nop
      real edgelx, at1, sum1
      double precision dsum1
      double precision ranorm
      idimp = size(part,1)
      nop = size(part,2)
! set boundary values
      edgelx = 0.0
      at1 = real(nx)/real(npx)
      if (ipbc.eq.2) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile and maxwellian velocity distribution
      do j = 1, npx
         part(1,j) = edgelx + at1*(real(j) - .5)
         part(2,j) = vtx*ranorm()
      enddo
! add correct drift
      dsum1 = 0.0d0
      do j = 1, npx
         dsum1 = dsum1 + dble(part(2,j))
      enddo
      sum1 = real(dsum1)/real(npx) - vdx
      part(2,:) = part(2,:) - sum1
      return
      end
!
      end module
!
!-----------------------------------------------------------------------
      function ranorm()
! this program calculates a random number y from a gaussian distribution
! with zero mean and unit variance, according to the method of
! mueller and box:
!    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
!    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
! where x is a random number uniformly distributed on (0,1).
! written for the ibm by viktor k. decyk, ucla
      implicit none
      integer iflg,isc,i1,r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end

