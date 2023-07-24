! Fortran Library for Skeleton 2-1/2D Electromagnetic PIC Code
! written by Viktor K. Decyk, UCLA
!-----------------------------------------------------------------------
      subroutine DISTR2H(part,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,idimp,nop,&
     &nx,ny,ipbc)
! for 2-1/2d code, this subroutine calculates initial particle
! co-ordinates and velocities with uniform density and maxwellian
! velocity with drift
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
! npx/npy = initial number of particles distributed in x/y direction
! idimp = size of phase space = 5
! nop = number of particles
! nx/ny = system length in x/y direction
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
! ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, idimp, nop, nx, ny, ipbc
      real vtx, vty, vtz, vdx, vdy, vdz
      real part
      dimension part(idimp,nop)
! local data
      integer j, k, k1, npxy
      real edgelx, edgely, at1, at2, at3, sum1, sum2, sum3
      double precision dsum1, dsum2, dsum3
      double precision ranorm
      npxy = npx*npy
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      at1 = real(nx)/real(npx)
      at2 = real(ny)/real(npy)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         at1 = real(nx-2)/real(npx)
         at2 = real(ny-2)/real(npy)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         at1 = real(nx-2)/real(npx)
      endif
! uniform density profile
      do 20 k = 1, npy
      k1 = npx*(k - 1)
      at3 = edgely + at2*(real(k) - 0.5)
      do 10 j = 1, npx
      part(1,j+k1) = edgelx + at1*(real(j) - 0.5)
      part(2,j+k1) = at3
   10 continue
   20 continue
! maxwellian velocity distribution
      do 30 j = 1, npxy
      part(3,j) = vtx*ranorm()
      part(4,j) = vty*ranorm()
      part(5,j) = vtz*ranorm()
   30 continue
! add correct drift
      dsum1 = 0.0d0
      dsum2 = 0.0d0
      dsum3 = 0.0d0
      do 40 j = 1, npxy
      dsum1 = dsum1 + part(3,j)
      dsum2 = dsum2 + part(4,j)
      dsum3 = dsum3 + part(5,j)
   40 continue
      sum1 = dsum1
      sum2 = dsum2
      sum3 = dsum3
      at1 = 1./real(npxy)
      sum1 = at1*sum1 - vdx
      sum2 = at1*sum2 - vdy
      sum3 = at1*sum3 - vdz
      do 50 j = 1, npxy
      part(3,j) = part(3,j) - sum1
      part(4,j) = part(4,j) - sum2
      part(5,j) = part(5,j) - sum3
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine GBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny,  &
     &nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. Using the Boris Mover.
! scalar version using guard cells
! 119 flops/particle, 1 divide, 29 loads, 5 stores
! input: all, output: part, ek
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t)), omy = (q/m)*by(x(t),y(t)), and
! omz = (q/m)*bz(x(t),y(t)).
! position equations used are:
! x(t+dt)=x(t) + vx(t+dt/2)*dt
! y(t+dt)=y(t) + vy(t+dt/2)*dt
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = velocity vx of particle n
! part(4,n) = velocity vy of particle n
! part(5,n) = velocity vz of particle n
! fxy(1,j,k) = x component of force/charge at grid (j,k)
! fxy(2,j,k) = y component of force/charge at grid (j,k)
! fxy(3,j,k) = z component of force/charge at grid (j,k)
! that is, convolution of electric field over particle shape
! bxy(1,j,k) = x component of magnetic field at grid (j,k)
! bxy(2,j,k) = y component of magnetic field at grid (j,k)
! bxy(3,j,k) = z component of magnetic field at grid (j,k)
! that is, the convolution of magnetic field over particle shape
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! kinetic energy/mass at time t is also calculated, using
! ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
!      (vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)
! idimp = size of phase space = 5
! nop = number of particles
! nx/ny = system length in x/y direction
! nxv = second dimension of field arrays, must be >= nx+1
! nyv = third dimension of field arrays, must be >= ny+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ek
      real part, fxy, bxy
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
! local data
      integer j, nn, mm, np, mp
      real qtmh, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, omxt, omyt, omzt, omt
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = 0.5*qbm*dt
      sum1 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
! find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1.0 - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
! find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    &
     &   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    &
     &   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp))                    &
     &   + amy*(dxp*fxy(3,np,mm) + amx*fxy(3,nn,mm))
! find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    &
     &   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    &
     &   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    &
     &   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
! time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
! calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
! new position
      dx = part(1,j) + dx*dtc
      dy = part(2,j) + dy*dtc
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
! normalize kinetic energy
      ek = ek + 0.5*sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine GRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx, &
     &ny,nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Boris Mover.
! scalar version using guard cells
! 131 flops/particle, 4 divides, 2 sqrts, 25 loads, 5 stores
! input: all, output: part, ek
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fx(x(t),y(t))*dt)
! py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fy(x(t),y(t))*dt)
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt) +
!    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt) +
!    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt) +
!    .5*(q/m)*fz(x(t),y(t))*dt)
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
!    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
!    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
!    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
! and om**2 = omx**2 + omy**2 + omz**2
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami,
! where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
! position equations used are:
! x(t+dt) = x(t) + px(t+dt/2)*dtg
! y(t+dt) = y(t) + py(t+dt/2)*dtg
! where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
! pz(t+dt/2)*pz(t+dt/2))*ci*ci)
! fx(x(t),y(t)), fy(x(t),y(t)), and fz(x(t),y(t))
! bx(x(t),y(t)), by(x(t),y(t)), and bz(x(t),y(t))
! are approximated by interpolation from the nearest grid points:
! fx(x,y) = (1-dy)*((1-dx)*fx(n,m)+dx*fx(n+1,m)) + dy*((1-dx)*fx(n,m+1)
!    + dx*fx(n+1,m+1))
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! similarly for fy(x,y), fz(x,y), bx(x,y), by(x,y), bz(x,y)
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = momentum px of particle n
! part(4,n) = momentum py of particle n
! part(5,n) = momentum pz of particle n
! fxy(1,j,k) = x component of force/charge at grid (j,k)
! fxy(2,j,k) = y component of force/charge at grid (j,k)
! fxy(3,j,k) = z component of force/charge at grid (j,k)
! that is, convolution of electric field over particle shape
! bxy(1,j,k) = x component of magnetic field at grid (j,k)
! bxy(2,j,k) = y component of magnetic field at grid (j,k)
! bxy(3,j,k) = z component of magnetic field at grid (j,k)
! that is, the convolution of magnetic field over particle shape
! qbm = particle charge/mass ratio
! dt = time interval between successive calculations
! dtc = time interval between successive co-ordinate calculations
! ci = reciprocal of velocity of light
! kinetic energy/mass at time t is also calculated, using
! ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
!      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
!      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
! idimp = size of phase space = 5
! nop = number of particles
! nx/ny = system length in x/y direction
! nxv = second dimension of field arrays, must be >= nx+1
! nyv = third dimension of field arrays, must be >= ny+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer idimp, nop, nx, ny, nxv, nyv, ipbc
      real qbm, dt, dtc, ci, ek
      real part, fxy, bxy
      dimension part(idimp,nop)
      dimension fxy(3,nxv,nyv), bxy(3,nxv,nyv)
! local data
      integer j, nn, mm, np, mp
      real qtmh, ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, dz, ox, oy, oz, acx, acy, acz, p2, gami, qtmg, dtg
      real omxt, omyt, omzt, omt, anorm
      real rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      double precision sum1
      qtmh = 0.5*qbm*dt
      ci2 = ci*ci
      sum1 = 0.0d0
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
! find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = part(1,j) - real(nn)
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = 1.0 - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
! find electric field
      dx = dyp*(dxp*fxy(1,np,mp) + amx*fxy(1,nn,mp))                    &
     &   + amy*(dxp*fxy(1,np,mm) + amx*fxy(1,nn,mm))
      dy = dyp*(dxp*fxy(2,np,mp) + amx*fxy(2,nn,mp))                    &
     &   + amy*(dxp*fxy(2,np,mm) + amx*fxy(2,nn,mm))
      dz = dyp*(dxp*fxy(3,np,mp) + amx*fxy(3,nn,mp))                    &
     &   + amy*(dxp*fxy(3,np,mm) + amx*fxy(3,nn,mm))
! calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
! half acceleration
      acx = part(3,j) + dx
      acy = part(4,j) + dy
      acz = part(5,j) + dz
! find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! find magnetic field
      ox = dyp*(dxp*bxy(1,np,mp) + amx*bxy(1,nn,mp))                    &
     &   + amy*(dxp*bxy(1,np,mm) + amx*bxy(1,nn,mm))
      oy = dyp*(dxp*bxy(2,np,mp) + amx*bxy(2,nn,mp))                    &
     &   + amy*(dxp*bxy(2,np,mm) + amx*bxy(2,nn,mm))
      oz = dyp*(dxp*bxy(3,np,mp) + amx*bxy(3,nn,mp))                    &
     &   + amy*(dxp*bxy(3,np,mm) + amx*bxy(3,nn,mm))
! renormalize magnetic field
      qtmg = qtmh*gami
! time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
! calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
! calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2.0/(1.0 + omt)
      omt = 0.5*(1.0 - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
! new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = dx
      part(4,j) = dy
      part(5,j) = dz
! update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = part(1,j) + dx*dtg
      dy = part(2,j) + dy*dtg
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
! normalize kinetic energy
      ek = ek + sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine GPOST2L(part,q,qm,nop,idimp,nxv,nyv)
! for 2d code, this subroutine calculates particle charge density
! using first-order linear interpolation, periodic boundaries
! scalar version using guard cells
! 17 flops/particle, 6 loads, 4 stores
! input: all, output: q
! charge density is approximated by values at the nearest grid points
! q(n,m)=qm*(1.-dx)*(1.-dy)
! q(n+1,m)=qm*dx*(1.-dy)
! q(n,m+1)=qm*(1.-dx)*dy
! q(n+1,m+1)=qm*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! q(j,k) = charge density at grid point j,k
! qm = charge on particle, in units of e
! nop = number of particles
! idimp = size of phase space = 4
! nxv = first dimension of charge array, must be >= nx+1
! nyv = second dimension of charge array, must be >= ny+1
      implicit none
      integer nop, idimp, nxv, nyv
      real qm
      real part, q
      dimension part(idimp,nop), q(nxv,nyv)
! local data
      integer j, nn, mm, np, mp
      real dxp, dyp, amx, amy
! find interpolation weights
      do 10 j = 1, nop
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1. - dyp
      np = nn + 1
! deposit charge
      q(np,mp) = q(np,mp) + dxp*dyp
      q(nn,mp) = q(nn,mp) + amx*dyp
      q(np,mm) = q(np,mm) + dxp*amy
      q(nn,mm) = q(nn,mm) + amx*amy
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine GJPOST2L(part,cu,qm,dt,nop,idimp,nx,ny,nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation
! in addition, particle positions are advanced a half time-step
! scalar version using guard cells
! 41 flops/particle, 17 loads, 14 stores
! input: all, output: part, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*vi, where i = x,y,z
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = x velocity of particle n
! part(4,n) = y velocity of particle n
! part(5,n) = z velocity of particle n
! cu(i,j,k) = ith component of current density at grid point j,k
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! nop = number of particles
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! nxv = first dimension of current array, must be >= nx+1
! nyv = second dimension of current array, must be >= ny+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt
      real part, cu
      dimension part(idimp,nop), cu(3,nxv,nyv)
! local data
      integer j, nn, mm, np, mp
      real edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
! find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
! deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      cu(1,np,mp) = cu(1,np,mp) + vx*dx
      cu(2,np,mp) = cu(2,np,mp) + vy*dx
      cu(3,np,mp) = cu(3,np,mp) + vz*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + vx*dx
      cu(2,np,mm) = cu(2,np,mm) + vy*dx
      cu(3,np,mm) = cu(3,np,mm) + vz*dx
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
! advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine GRJPOST2L(part,cu,qm,dt,ci,nop,idimp,nx,ny,nxv,nyv,ipbc&
     &)
! for 2-1/2d code, this subroutine calculates particle current density
! using first-order linear interpolation for relativistic particles
! in addition, particle positions are advanced a half time-step
! scalar version using guard cells
! 47 flops/particle, 1 divide, 1 sqrt, 17 loads, 14 stores
! input: all, output: part, cu
! current density is approximated by values at the nearest grid points
! cu(i,n,m)=qci*(1.-dx)*(1.-dy)
! cu(i,n+1,m)=qci*dx*(1.-dy)
! cu(i,n,m+1)=qci*(1.-dx)*dy
! cu(i,n+1,m+1)=qci*dx*dy
! where n,m = leftmost grid points and dx = x-n, dy = y-m
! and qci = qm*pi*gami, where i = x,y,z
! where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
! part(1,n) = position x of particle n
! part(2,n) = position y of particle n
! part(3,n) = x momentum of particle n
! part(4,n) = y momentum of particle n
! part(5,n) = z momentum of particle n
! cu(i,j,k) = ith component of current density at grid point j,k
! qm = charge on particle, in units of e
! dt = time interval between successive calculations
! ci = reciprocal of velocity of light
! nop = number of particles
! idimp = size of phase space = 5
! nx/ny = system length in x/y direction
! nxv = first dimension of current array, must be >= nx+1
! nyv = second dimension of current array, must be >= ny+1
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      integer nop, idimp, nx, ny, nxv, nyv, ipbc
      real qm, dt, ci
      real part, cu
      dimension part(idimp,nop), cu(3,nxv,nyv)
! local data
      integer j, nn, mm, np, mp
      real ci2, edgelx, edgely, edgerx, edgery, dxp, dyp, amx, amy
      real dx, dy, vx, vy, vz, p2, gami
      ci2 = ci*ci
! set boundary values
      edgelx = 0.0
      edgely = 0.0
      edgerx = real(nx)
      edgery = real(ny)
      if (ipbc.eq.2) then
         edgelx = 1.0
         edgely = 1.0
         edgerx = real(nx-1)
         edgery = real(ny-1)
      else if (ipbc.eq.3) then
         edgelx = 1.0
         edgerx = real(nx-1)
      endif
      do 10 j = 1, nop
! find interpolation weights
      nn = part(1,j)
      mm = part(2,j)
      dxp = qm*(part(1,j) - real(nn))
      dyp = part(2,j) - real(mm)
! find inverse gamma
      vx = part(3,j)
      vy = part(4,j)
      vz = part(5,j)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
! calculate weights
      nn = nn + 1
      mm = mm + 1
      amx = qm - dxp
      mp = mm + 1
      amy = 1.0 - dyp
      np = nn + 1
! deposit current
      dx = dxp*dyp
      dy = amx*dyp
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
      cu(1,np,mp) = cu(1,np,mp) + vx*dx
      cu(2,np,mp) = cu(2,np,mp) + vy*dx
      cu(3,np,mp) = cu(3,np,mp) + vz*dx
      dx = dxp*amy
      cu(1,nn,mp) = cu(1,nn,mp) + vx*dy
      cu(2,nn,mp) = cu(2,nn,mp) + vy*dy
      cu(3,nn,mp) = cu(3,nn,mp) + vz*dy
      dy = amx*amy
      cu(1,np,mm) = cu(1,np,mm) + vx*dx
      cu(2,np,mm) = cu(2,np,mm) + vy*dx
      cu(3,np,mm) = cu(3,np,mm) + vz*dx
      cu(1,nn,mm) = cu(1,nn,mm) + vx*dy
      cu(2,nn,mm) = cu(2,nn,mm) + vy*dy
      cu(3,nn,mm) = cu(3,nn,mm) + vz*dy
! advance position half a time-step
      dx = part(1,j) + vx*dt
      dy = part(2,j) + vy*dt
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j)
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j)
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   10 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine DSORTP2YL(parta,partb,npic,idimp,nop,ny1)
! this subroutine sorts particles by y grid
! linear interpolation
! parta/partb = input/output particle arrays
! parta(2,n) = position y of particle n
! npic = address offset for reordering particles
! idimp = size of phase space = 4
! nop = number of particles
! ny1 = system length in y direction + 1
      implicit none
      integer npic, idimp, nop, ny1
      real parta, partb
      dimension parta(idimp,nop), partb(idimp,nop), npic(ny1)
! local data
      integer i, j, k, m, isum, ist, ip
! clear counter array
      do 10 k = 1, ny1
      npic(k) = 0
   10 continue
! find how many particles in each grid
      do 20 j = 1, nop
      m = parta(2,j)
      m = m + 1
      npic(m) = npic(m) + 1
   20 continue
! find address offset
      isum = 0
      do 30 k = 1, ny1
      ist = npic(k)
      npic(k) = isum
      isum = isum + ist
   30 continue
! find addresses of particles at each grid and reorder particles
      do 50 j = 1, nop
      m = parta(2,j)
      m = m + 1
      ip = npic(m) + 1
      do 40 i = 1, idimp
      partb(i,ip) = parta(i,j)
   40 continue
      npic(m) = ip
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine BGUARD2L(bxy,nx,ny,nxe,nye)
! replicate extended periodic vector field bxy
! linear interpolation
! nx/ny = system length in x/y direction
! nxe = first dimension of field arrays, must be >= nx+1
! nye = second dimension of field arrays, must be >= ny+1
      implicit none
      real bxy
      integer nx, ny, nxe, nye
      dimension bxy(3,nxe,nye)
! local data
      integer j, k
! copy edges of extended field
      do 10 k = 1, ny
      bxy(1,nx+1,k) = bxy(1,1,k)
      bxy(2,nx+1,k) = bxy(2,1,k)
      bxy(3,nx+1,k) = bxy(3,1,k)
   10 continue
      do 20 j = 1, nx
      bxy(1,j,ny+1) = bxy(1,j,1)
      bxy(2,j,ny+1) = bxy(2,j,1)
      bxy(3,j,ny+1) = bxy(3,j,1)
   20 continue
      bxy(1,nx+1,ny+1) = bxy(1,1,1)
      bxy(2,nx+1,ny+1) = bxy(2,1,1)
      bxy(3,nx+1,ny+1) = bxy(3,1,1)
      return
      end
!-----------------------------------------------------------------------
      subroutine ACGUARD2L(cu,nx,ny,nxe,nye)
! accumulate extended periodic vector field cu
! linear interpolation
! nx/ny = system length in x/y direction
! nxe = first dimension of field arrays, must be >= nx+1
! nye = second dimension of field arrays, must be >= ny+1
      implicit none
      real cu
      integer nx, ny, nxe, nye
      dimension cu(3,nxe,nye)
! local data
      integer i, j, k
! accumulate edges of extended field
      do 20 k = 1, ny
      do 10 i = 1, 3
      cu(i,1,k) = cu(i,1,k) + cu(i,nx+1,k)
      cu(i,nx+1,k) = 0.0
   10 continue
   20 continue
      do 40 j = 1, nx
      do 30 i = 1, 3
      cu(i,j,1) = cu(i,j,1) + cu(i,j,ny+1)
      cu(i,j,ny+1) = 0.0
   30 continue
   40 continue
      do 50 i = 1, 3
      cu(i,1,1) = cu(i,1,1) + cu(i,nx+1,ny+1)
      cu(i,nx+1,ny+1) = 0.0
   50 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine AGUARD2L(q,nx,ny,nxe,nye)
! accumulate extended periodic scalar field q
! linear interpolation
! nx/ny = system length in x/y direction
! nxe = first dimension of field arrays, must be >= nx+1
! nye = second dimension of field arrays, must be >= ny+1
      implicit none
      real q
      integer nx, ny, nxe, nye
      dimension q(nxe,nye)
! local data
      integer j, k
! accumulate edges of extended field
      do 10 k = 1, ny
      q(1,k) = q(1,k) + q(nx+1,k)
      q(nx+1,k) = 0.0
   10 continue
      do 20 j = 1, nx
      q(j,1) = q(j,1) + q(j,ny+1)
      q(j,ny+1) = 0.0
   20 continue
      q(1,1) = q(1,1) + q(nx+1,ny+1)
      q(nx+1,ny+1) = 0.0
      return
      end
!-----------------------------------------------------------------------
      subroutine POIS23(q,fxy,isign,ffc,ax,ay,affp,we,nx,ny,nxvh,nyv,   &
     &nxhd,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! force/charge (or convolution of electric field over particle shape)
! with periodic boundary conditions.  Zeros out z component.
! for isign = 0, input: isign,ax,ay,affp,nx,ny,nxvh,nxhd,nyhd,
! output: ffc
! for isign /= 0, input: q,ffc,isign,nx,ny,nxvh,nxhd,nyhd,
! output: fxy,we
! approximate flop count is: 26*nxc*nyc + 12*(nxc + nyc)
! where nxc = nx/2 - 1, nyc = ny/2 - 1
! equation used is:
! fx(kx,ky) = -sqrt(-1)*kx*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fy(kx,ky) = -sqrt(-1)*ky*g(kx,ky)*s(kx,ky)*q(kx,ky),
! fz(kx,ky) = zero,
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! fx(kx=pi) = fy(kx=pi) = fx(ky=pi) = fy(ky=pi) = 0, and
! fx(kx=0,ky=0) = fy(kx=0,ky=0) = 0.
! q(j,k) = complex charge density for fourier mode (j-1,k-1)
! fxy(1,j,k) = x component of complex force/charge,
! fxy(2,j,k) = y component of complex force/charge,
! fxy(3,j,k) = zero,
! all for fourier mode (j-1,k-1)
! if isign = 0, form factor array is prepared
! if isign is not equal to 0, force/charge is calculated
! aimag(ffc(j,k)) = finite-size particle shape factor s
! for fourier mode (j-1,k-1)
! real(ffc(j,k)) = potential green's function g
! for fourier mode (j-1,k-1)
! ax/ay = half-width of particle in x/y direction
! affp = normalization constant = nx*ny/np, where np=number of particles
! electric field energy is also calculated, using
! we = nx*ny*sum((affp/(kx**2+ky**2))*|q(kx,ky)*s(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nxvh = first dimension of field arrays, must be >= nxh
! nyv = second dimension of field arrays, must be >= ny
! nxhd = first dimension of form factor array, must be >= nxh
! nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      real ax, ay, affp, we
      complex q, fxy, ffc
      dimension q(nxvh,nyv), fxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
! local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, at1, at2, at3, at4
      complex zero, zt1, zt2
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      if (isign.ne.0) go to 30
! prepare form factor array
      do 20 k = 1, nyh
      dky = dny*real(k - 1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 10 j = 1, nxh
      dkx = dnx*real(j - 1)
      at3 = dkx*dkx + at1
      at4 = exp(-.5*((dkx*ax)**2 + at2))
      if (at3.eq.0.0) then
         ffc(j,k) = cmplx(affp,1.0)
      else
         ffc(j,k) = cmplx(affp*at4/at3,at4)
      endif
   10 continue
   20 continue
      return
! calculate force/charge and sum field energy
   30 wp = 0.0d0
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 50 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 40 j = 2, nxh
      at1 = real(ffc(j,k))*aimag(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      zt1 = cmplx(aimag(q(j,k)),-real(q(j,k)))
      zt2 = cmplx(aimag(q(j,k1)),-real(q(j,k1)))
      fxy(1,j,k) = at2*zt1
      fxy(2,j,k) = at3*zt1
      fxy(3,j,k) = zero
      fxy(1,j,k1) = at2*zt2
      fxy(2,j,k1) = -at3*zt2
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,k)*conjg(q(j,k)) + q(j,k1)*conjg(q(j,k1)))
   40 continue
   50 continue
! mode numbers kx = 0, nx/2
      do 60 k = 2, nyh
      k1 = ny2 - k
      at1 = real(ffc(1,k))*aimag(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      zt1 = cmplx(aimag(q(1,k)),-real(q(1,k)))
      fxy(1,1,k) = zero
      fxy(2,1,k) = at3*zt1
      fxy(3,1,k) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      wp = wp + at1*(q(1,k)*conjg(q(1,k)))
   60 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 70 j = 2, nxh
      at1 = real(ffc(j,1))*aimag(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      zt1 = cmplx(aimag(q(j,1)),-real(q(j,1)))
      fxy(1,j,1) = at2*zt1
      fxy(2,j,1) = zero
      fxy(3,j,1) = zero
      fxy(1,j,k1) = zero
      fxy(2,j,k1) = zero
      fxy(3,j,k1) = zero
      wp = wp + at1*(q(j,1)*conjg(q(j,1)))
   70 continue
      fxy(1,1,1) = zero
      fxy(2,1,1) = zero
      fxy(3,1,1) = zero
      fxy(1,1,k1) = zero
      fxy(2,1,k1) = zero
      fxy(3,1,k1) = zero
      we = real(nx*ny)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine CUPERP2(cu,nx,ny,nxvh,nyv)
! this subroutine calculates the transverse current in fourier space
! input: all, output: cu
! approximate flop count is: 36*nxc*nyc
! and nxc*nyc divides
! where nxc = nx/2 - 1, nyc = ny/2 - 1
! the transverse current is calculated using the equation:
! cux(kx,ky) = cux(kx,ky)-kx*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! cuy(kx,ky) = cuy(kx,ky)-ky*(kx*cux(kx,ky)+ky*cuy(kx,ky))/(kx*kx+ky*ky)
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! except for cux(kx=pi) = cuy(kx=pi) = 0, cux(ky=pi) = cuy(ky=pi) = 0,
! and cux(kx=0,ky=0) = cuy(kx=0,ky=0) = 0.
! cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
! nx/ny = system length in x/y direction
! nxvh = first dimension of current array, must be >= nxh
! nyv = second dimension of current array, must be >= ny
      implicit none
      integer nx, ny, nxvh, nyv
      complex cu
      dimension cu(3,nxvh,nyv)
! local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dkx, dky, dky2, at1
      complex zero, zt1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
! calculate transverse part of current
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      dky2 = dky*dky
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      at1 = 1./(dkx*dkx + dky2)
      zt1 = at1*(dkx*cu(1,j,k) + dky*cu(2,j,k))
      cu(1,j,k) = cu(1,j,k) - dkx*zt1
      cu(2,j,k) = cu(2,j,k) - dky*zt1
      zt1 = at1*(dkx*cu(1,j,k1) - dky*cu(2,j,k1))
      cu(1,j,k1) = cu(1,j,k1) - dkx*zt1
      cu(2,j,k1) = cu(2,j,k1) + dky*zt1
   10 continue
   20 continue
! mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      cu(2,1,k) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
   30 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      cu(1,j,1) = zero
      cu(1,j,k1) = zero
      cu(2,j,k1) = zero
   40 continue
      cu(1,1,1) = zero
      cu(2,1,1) = zero
      cu(1,1,k1) = zero
      cu(2,1,k1) = zero
      return
      end
!-----------------------------------------------------------------------
      subroutine IBPOIS23(cu,bxy,ffc,ci,wm,nx,ny,nxvh,nyv,nxhd,nyhd)
! this subroutine solves 2-1/2d poisson's equation in fourier space for
! magnetic field, with periodic boundary conditions.
! input: cu,ffc,ci,nx,ny,nxvh,nxhd,nyhd, output: bxy,wm
! approximate flop count is: 90*nxc*nyc + 40*(nxc + nyc)
! where nxc = nx/2 - 1, nyc = ny/2 - 1
! the magnetic field is calculated using the equations:
! bx(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*ky*cuz(kx,ky),
! by(kx,ky) = -ci*ci*sqrt(-1)*g(kx,ky)*kx*cuz(kx,ky),
! bz(kx,ky) = ci*ci*sqrt(-1)*g(kx,ky)*(kx*cuy(kx,ky)-ky*cux(kx,ky)),
! where kx = 2pi*j/nx, ky = 2pi*k/ny, and j,k = fourier mode numbers,
! g(kx,ky) = (affp/(kx**2+ky**2))*s(kx,ky),
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2), except for
! bx(kx=pi) = by(kx=pi) = bz(kx=pi) = bx(ky=pi) = by(ky=pi) = bz(ky=pi) 
! = 0, and bx(kx=0,ky=0) = by(kx=0,ky=0) = bz(kx=0,ky=0) = 0.
! cu(i,j,k) = complex current density for fourier mode (j-1,k-1)
! bxy(i,j,k) = i component of complex magnetic field
! all for fourier mode (j-1,k-1)
! aimag(ffc(j,k)) = finite-size particle shape factor s
! for fourier mode (j-1,k-1)
! real(ffc(j,k)) = potential green's function g
! for fourier mode (j-1,k-1)
! ci = reciprocal of velocity of light
! magnetic field energy is also calculated, using
! wm = nx*ny*sum((affp/(kx**2+ky**2))*ci*ci*
!    |cu(kx,ky)*s(kx,ky)|**2), where
! affp = normalization constant = nx*ny/np, where np=number of particles
! this expression is valid only if the current is divergence-free
! nx/ny = system length in x/y direction
! nxvh = first dimension of field arrays, must be >= nxh
! nyv = second dimension of field arrays, must be >= ny
! nxhd = first dimension of form factor array, must be >= nxh
! nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, wm
      complex cu, bxy, ffc
      dimension cu(3,nxvh,nyv), bxy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
! local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dky, ci2, at1, at2, at3
      complex zero, zt1, zt2, zt3
      double precision wp
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      zero = cmplx(0.0,0.0)
      ci2 = ci*ci
! calculate magnetic field and sum field energy
      wp = 0.0d0
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      at1 = ci2*real(ffc(j,k))
      at2 = dnx*real(j - 1)*at1
      at3 = dky*at1
      at1 = at1*aimag(ffc(j,k))
      zt1 = cmplx(-aimag(cu(3,j,k)),real(cu(3,j,k)))
      zt2 = cmplx(-aimag(cu(2,j,k)),real(cu(2,j,k)))
      zt3 = cmplx(-aimag(cu(1,j,k)),real(cu(1,j,k)))
      bxy(1,j,k) = at3*zt1
      bxy(2,j,k) = -at2*zt1
      bxy(3,j,k) = at2*zt2 - at3*zt3
      zt1 = cmplx(-aimag(cu(3,j,k1)),real(cu(3,j,k1)))
      zt2 = cmplx(-aimag(cu(2,j,k1)),real(cu(2,j,k1)))
      zt3 = cmplx(-aimag(cu(1,j,k1)),real(cu(1,j,k1)))
      bxy(1,j,k1) = -at3*zt1
      bxy(2,j,k1) = -at2*zt1
      bxy(3,j,k1) = at2*zt2 + at3*zt3
      wp = wp + at1*(cu(1,j,k)*conjg(cu(1,j,k))                         &
     &   + cu(2,j,k)*conjg(cu(2,j,k)) + cu(3,j,k)*conjg(cu(3,j,k))      &
     &   + cu(1,j,k1)*conjg(cu(1,j,k1)) + cu(2,j,k1)*conjg(cu(2,j,k1))  &
     &   + cu(3,j,k1)*conjg(cu(3,j,k1)))
   10 continue
   20 continue
! mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      at1 = ci2*real(ffc(1,k))
      at3 = dny*real(k - 1)*at1
      at1 = at1*aimag(ffc(1,k))
      zt1 = cmplx(-aimag(cu(3,1,k)),real(cu(3,1,k)))
      zt3 = cmplx(-aimag(cu(1,1,k)),real(cu(1,1,k)))
      bxy(1,1,k) = at3*zt1
      bxy(2,1,k) = zero
      bxy(3,1,k) = -at3*zt3
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wp = wp + at1*(cu(1,1,k)*conjg(cu(1,1,k))                         &
     &   + cu(2,1,k)*conjg(cu(2,1,k)) + cu(3,1,k)*conjg(cu(3,1,k)))
   30 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      at1 = ci2*real(ffc(j,1))
      at2 = dnx*real(j - 1)*at1
      at1 = at1*aimag(ffc(j,1))
      zt1 = cmplx(-aimag(cu(3,j,1)),real(cu(3,j,1)))
      zt2 = cmplx(-aimag(cu(2,j,1)),real(cu(2,j,1)))
      bxy(1,j,1) = zero
      bxy(2,j,1) = -at2*zt1
      bxy(3,j,1) = at2*zt2
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      wp = wp + at1*(cu(1,j,1)*conjg(cu(1,j,1))                         &
     &   + cu(2,j,1)*conjg(cu(2,j,1)) + cu(3,j,1)*conjg(cu(3,j,1)))
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      wm = real(nx*ny)*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine MAXWEL2(exy,bxy,cu,ffc,ci,dt,wf,wm,nx,ny,nxvh,nyv,nxhd,&
     &nyhd)
! this subroutine solves 2-1/2d maxwell's equation in fourier space for
! transverse electric and magnetic fields with periodic boundary
! conditions.
! input: all, output: wf, wm, exy, bxy
! approximate flop count is: 286*nxc*nyc + 84*(nxc + nyc)
! where nxc = nx/2 - 1, nyc = ny/2 - 1
! the magnetic field is first updated half a step using the equations:
! bx(kx,ky) = bx(kx,ky) - .5*dt*sqrt(-1)*ky*ez(kx,ky)
! by(kx,ky) = by(kx,ky) + .5*dt*sqrt(-1)*kx*ez(kx,ky)
! bz(kx,ky) = bz(kx,ky) - .5*dt*sqrt(-1)*(kx*ey(kx,ky)-ky*ex(kx,ky))
! the electric field is then updated a whole step using the equations:
! ex(kx,ky) = ex(kx,ky) + c2*dt*sqrt(-1)*ky*bz(kx,ky)
!                       - affp*dt*cux(kx,ky)*s(kx,ky)
! ey(kx,ky) = ey(kx,ky) - c2*dt*sqrt(-1)*kx*bz(kx,ky)
!                       - affp*dt*cuy(kx,ky)*s(kx,ky)
! ez(kx,ky) = ez(kx,ky) + c2*dt*sqrt(-1)*(kx*by(kx,ky)-ky*bx(kx,ky))
!                       - affp*dt*cuz(kx,ky)*s(kx,ky)
! the magnetic field is finally updated the remaining half step with
! the new electric field and the previous magnetic field equations.
! where kx = 2pi*j/nx, ky = 2pi*k/ny, c2 = 1./(ci*ci)
! and s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)
! j,k = fourier mode numbers, except for
! ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
! ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
! ex(kx=0,ky=0) = ey(kx=0,ky=0) = ez(kx=0,ky=0) = 0.
! and similarly for bx, by, bz.
! cu(i,j,k) = complex current density
! exy(i,j,k) = complex transverse electric field
! bxy(i,j,k) = complex magnetic field
! for component i, all for fourier mode (j-1,k-1)
! real(ffc(1,1)) = affp = normalization constant = nx*ny/np,
! where np=number of particles
! aimag(ffc(j,k)) = finite-size particle shape factor s,
! s(kx,ky) = exp(-((kx*ax)**2+(ky*ay)**2)/2)
! for fourier mode (j-1,k-1)
! ci = reciprocal of velocity of light
! dt = time interval between successive calculations
! transverse electric field energy is also calculated, using
! wf = nx*ny**sum((1/affp)*|exy(kx,ky)|**2)
! magnetic field energy is also calculated, using
! wm = nx*ny**sum((c2/affp)*|bxy(kx,ky)|**2)
! nx/ny = system length in x/y direction
! nxvh = first dimension of field arrays, must be >= nxh
! nyv = second dimension of field arrays, must be >= ny
! nxhd = first dimension of form factor array, must be >= nxh
! nyhd = second dimension of form factor array, must be >= nyh
      implicit none
      integer nx, ny, nxvh, nyv, nxhd, nyhd
      real ci, dt, wf, wm
      complex exy, bxy, cu, ffc
      dimension exy(3,nxvh,nyv), bxy(3,nxvh,nyv), cu(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
! local data
      integer nxh, nyh, ny2, j, k, k1
      real dnx, dny, dth, c2, cdt, affp, anorm, dkx, dky, afdt, adt
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      double precision wp, ws
      if (ci.le.0.0) return
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
      dnx = 6.28318530717959/real(nx)
      dny = 6.28318530717959/real(ny)
      dth = 0.5*dt
      c2 = 1.0/(ci*ci)
      cdt = c2*dt
      affp = real(ffc(1,1))
      adt = affp*dt
      zero = cmplx(0.0,0.0)
      anorm = 1.0/affp
! update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
! calculate the electromagnetic fields
! mode numbers 0 < kx < nx/2 and 0 < ky < ny/2
      do 20 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      do 10 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,k))
! update magnetic field half time step, ky > 0
      zt1 = cmplx(-aimag(exy(3,j,k)),real(exy(3,j,k)))
      zt2 = cmplx(-aimag(exy(2,j,k)),real(exy(2,j,k)))
      zt3 = cmplx(-aimag(exy(1,j,k)),real(exy(1,j,k)))
      zt4 = bxy(1,j,k) - dth*(dky*zt1)
      zt5 = bxy(2,j,k) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k) - dth*(dkx*zt2 - dky*zt3)
! update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k) + cdt*(dky*zt1) - afdt*cu(1,j,k)
      zt8 = exy(2,j,k) - cdt*(dkx*zt1) - afdt*cu(2,j,k)
      zt9 = exy(3,j,k) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,j,k)
! update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k) = zt7
      exy(2,j,k) = zt8
      exy(3,j,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
      bxy(1,j,k) = zt4
      bxy(2,j,k) = zt5
      bxy(3,j,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
! update magnetic field half time step, ky < 0
      zt1 = cmplx(-aimag(exy(3,j,k1)),real(exy(3,j,k1)))
      zt2 = cmplx(-aimag(exy(2,j,k1)),real(exy(2,j,k1)))
      zt3 = cmplx(-aimag(exy(1,j,k1)),real(exy(1,j,k1)))
      zt4 = bxy(1,j,k1) + dth*(dky*zt1)
      zt5 = bxy(2,j,k1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,k1) - dth*(dkx*zt2 + dky*zt3)
! update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,j,k1) - cdt*(dky*zt1) - afdt*cu(1,j,k1)
      zt8 = exy(2,j,k1) - cdt*(dkx*zt1) - afdt*cu(2,j,k1)
      zt9 = exy(3,j,k1) + cdt*(dkx*zt2 + dky*zt3) - afdt*cu(3,j,k1)
! update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,j,k1) = zt7
      exy(2,j,k1) = zt8
      exy(3,j,k1) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg(zt9))
      zt4 = zt4 + dth*(dky*zt1)
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2 + dky*zt3)
      bxy(1,j,k1) = zt4
      bxy(2,j,k1) = zt5
      bxy(3,j,k1) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg(zt6))
   10 continue
   20 continue
! mode numbers kx = 0, nx/2
      do 30 k = 2, nyh
      k1 = ny2 - k
      dky = dny*real(k - 1)
      afdt = adt*aimag(ffc(1,k))
! update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,1,k)),real(exy(3,1,k)))
      zt3 = cmplx(-aimag(exy(1,1,k)),real(exy(1,1,k)))
      zt4 = bxy(1,1,k) - dth*(dky*zt1)
      zt6 = bxy(3,1,k) + dth*(dky*zt3)
! update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt3 = cmplx(-aimag(zt4),real(zt4))
      zt7 = exy(1,1,k) + cdt*(dky*zt1) - afdt*cu(1,1,k)
      zt9 = exy(3,1,k) - cdt*(dky*zt3) - afdt*cu(3,1,k)
! update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt3 = cmplx(-aimag(zt7),real(zt7))
      exy(1,1,k) = zt7
      exy(2,1,k) = zero
      exy(3,1,k) = zt9
      ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
      zt4 = zt4 - dth*(dky*zt1)
      zt6 = zt6 + dth*(dky*zt3)
      bxy(1,1,k) = zt4
      bxy(2,1,k) = zero
      bxy(3,1,k) = zt6
      wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
   30 continue
! mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 40 j = 2, nxh
      dkx = dnx*real(j - 1)
      afdt = adt*aimag(ffc(j,1))
! update magnetic field half time step
      zt1 = cmplx(-aimag(exy(3,j,1)),real(exy(3,j,1)))
      zt2 = cmplx(-aimag(exy(2,j,1)),real(exy(2,j,1)))
      zt5 = bxy(2,j,1) + dth*(dkx*zt1)
      zt6 = bxy(3,j,1) - dth*(dkx*zt2)
! update electric field whole time step
      zt1 = cmplx(-aimag(zt6),real(zt6))
      zt2 = cmplx(-aimag(zt5),real(zt5))
      zt8 = exy(2,j,1) - cdt*(dkx*zt1) - afdt*cu(2,j,1)
      zt9 = exy(3,j,1) + cdt*(dkx*zt2) - afdt*cu(3,j,1)
! update magnetic field half time step and store electric field
      zt1 = cmplx(-aimag(zt9),real(zt9))
      zt2 = cmplx(-aimag(zt8),real(zt8))
      exy(1,j,1) = zero
      exy(2,j,1) = zt8
      exy(3,j,1) = zt9
      ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
      zt5 = zt5 + dth*(dkx*zt1)
      zt6 = zt6 - dth*(dkx*zt2)
      bxy(1,j,1) = zero
      bxy(2,j,1) = zt5
      bxy(3,j,1) = zt6
      wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
      bxy(1,j,k1) = zero
      bxy(2,j,k1) = zero
      bxy(3,j,k1) = zero
      exy(1,j,k1) = zero
      exy(2,j,k1) = zero
      exy(3,j,k1) = zero
   40 continue
      bxy(1,1,1) = zero
      bxy(2,1,1) = zero
      bxy(3,1,1) = zero
      exy(1,1,1) = zero
      exy(2,1,1) = zero
      exy(3,1,1) = zero
      bxy(1,1,k1) = zero
      bxy(2,1,k1) = zero
      bxy(3,1,k1) = zero
      exy(1,1,k1) = zero
      exy(2,1,k1) = zero
      exy(3,1,k1) = zero
      wf = real(nx*ny)*ws
      wm = real(nx*ny)*c2*wp
      return
      end
!-----------------------------------------------------------------------
      subroutine EMFIELD2(fxy,exy,ffc,isign,nx,ny,nxvh,nyv,nxhd,nyhd)
! this subroutine either adds complex vector fields if isign > 0
! or copies complex vector fields if isign < 0
! includes additional smoothing
      implicit none
      integer isign, nx, ny, nxvh, nyv, nxhd, nyhd
      complex fxy, exy, ffc
      dimension fxy(3,nxvh,nyv), exy(3,nxvh,nyv)
      dimension ffc(nxhd,nyhd)
! local data
      integer i, j, k, nxh, nyh, ny2, k1
      real at1
      nxh = nx/2
      nyh = max(1,ny/2)
      ny2 = ny + 2
! add the fields
      if (isign.gt.0) then
         do 30 k = 2, nyh
         k1 = ny2 - k
         do 20 j = 1, nxh
         at1 = aimag(ffc(j,k))
         do 10 i = 1, 3
         fxy(i,j,k) = fxy(i,j,k) + exy(i,j,k)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
   10    continue
   20    continue
   30    continue
         k1 = nyh + 1
         do 50 j = 1, nxh
         at1 = aimag(ffc(j,1))
         do 40 i = 1, 3
         fxy(i,j,1) = fxy(i,j,1) + exy(i,j,1)*at1
         fxy(i,j,k1) = fxy(i,j,k1) + exy(i,j,k1)*at1
   40    continue
   50    continue
! copy the fields
      else if (isign.lt.0) then
         do 80 k = 2, nyh
         k1 = ny2 - k
         do 70 j = 1, nxh
         at1 = aimag(ffc(j,k))
         do 60 i = 1, 3
         fxy(i,j,k) = exy(i,j,k)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   60    continue
   70    continue
   80    continue
         k1 = nyh + 1
         do 100 j = 1, nxh
         at1 = aimag(ffc(j,1))
         do 90 i = 1, 3
         fxy(i,j,1) = exy(i,j,1)*at1
         fxy(i,j,k1) = exy(i,j,k1)*at1
   90    continue
  100    continue
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RINIT(mixup,sct,indx,indy,nxhyd,nxyhd)
! this subroutine calculates tables needed by a two dimensional
! real to complex fast fourier transform and its inverse.
! input: indx, indy, nxhyd, nxyhd
! output: mixup, sct
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! nxhyd = maximum of (nx/2,ny)
! nxyhd = one half of maximum of (nx,ny)
! written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, nxhyd, nxyhd
      integer mixup
      complex sct
      dimension mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nxy, nxhy, nxyh
      integer j, k, lb, ll, jb, it
      real dnxy, arg
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
! bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhy
      lb = j - 1
      ll = 0
      do 10 k = 1, indx1y
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
! sine/cosine table for the angles 2*n*pi/nxy
      nxyh = nxy/2
      dnxy = 6.28318530717959/real(nxy)
      do 30 j = 1, nxyh
      arg = dnxy*real(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2RX(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,    &
     &nxyhd)
! wrapper function for real to complex fft, with packed data
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,&
     &nxyhd)
! perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd&
     &,nxyhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform y fft
         call FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd&
     &,nxyhd)
! perform x fft
         call FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,&
     &nxyhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine WFFT2R3(f,isign,mixup,sct,indx,indy,nxhd,nyd,nxhyd,    &
     &nxyhd)
! wrapper function for 3 2d real to complex ffts
      implicit none
      complex f, sct
      integer mixup
      integer isign, indx, indy, nxhd, nyd, nxhyd, nxyhd
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer nxh, ny, nxi, nyi
      data nxi, nyi /1,1/
! calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
! inverse fourier transform
      if (isign.lt.0) then
! perform x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,&
     &nxyhd)
! perform y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd&
     &,nxyhd)
! forward fourier transform
      else if (isign.gt.0) then
! perform y fft
         call FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxh,nxhd,nyd,nxhyd&
     &,nxyhd)
! perform x fft
         call FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,ny,nxhd,nyd,nxhyd,&
     &nxyhd)
      endif
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RXX(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,  &
     &nxhyd,nxyhd)
! this subroutine performs the x part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of y,
! using complex arithmetic.
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(n,m) = (1/nx*ny)*sum(f(j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = first dimension of f >= nx/2
! nyd = second dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
! except for f(1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! aimag(f(1,1)) = real part of mode nx/2,0 and
! aimag(f(1,ny/2+1)) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr
      real ani
      complex t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 100
! inverse fourier transform
! bit-reverse array elements in x
      nrx = nxhy/nxh
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
! first transform in x
      nrx = nxy/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 80 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 70 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = ani*(t1 + t2)
      f(nxh2-j,k) = ani*conjg(t1 - t2)
   70 continue
   80 continue
      ani = 2.0*ani
      do 90 k = nyi, nyt
      f(nxhh+1,k) = ani*conjg(f(nxhh+1,k))
      f(1,k) = ani*cmplx(real(f(1,k)) + aimag(f(1,k)),                  &
     &                   real(f(1,k)) - aimag(f(1,k)))
   90 continue
      return
! forward fourier transform
! scramble coefficients
  100 kmr = nxy/nx
      do 120 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 110 k = nyi, nyt
      t2 = conjg(f(nxh2-j,k))
      t1 = f(j,k) + t2
      t2 = (f(j,k) - t2)*t3
      f(j,k) = t1 + t2
      f(nxh2-j,k) = conjg(t1 - t2)
  110 continue
  120 continue
      do 130 k = nyi, nyt
      f(nxhh+1,k) = 2.0*conjg(f(nxhh+1,k))
      f(1,k) = cmplx(real(f(1,k)) + aimag(f(1,k)),                      &
     &               real(f(1,k)) - aimag(f(1,k)))
  130 continue
! bit-reverse array elements in x
      nrx = nxhy/nxh
      do 150 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 150
      do 140 k = nyi, nyt
      t1 = f(j1,k)
      f(j1,k) = f(j,k)
      f(j,k) = t1
  140 continue
  150 continue
! then transform in x
      nrx = nxy/nxh
      do 190 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 170 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 160 i = nyi, nyt
      t2 = t1*f(j2,i)
      f(j2,i) = f(j1,i) - t2
      f(j1,i) = f(j1,i) + t2
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2RXY(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  &
     &nxhyd,nxyhd)
! this subroutine performs the y part of a two dimensional real to
! complex fast fourier transform and its inverse, for a subset of x,
! using complex arithmetic
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, an inverse fourier transform is performed
! f(n,m) = (1/nx*ny)*sum(f(j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, a forward fourier transform is performed
! f(j,k) = sum(f(n,m)*exp(sqrt(-1)*2pi*n*j/nx)*exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxi = initial x index used
! nxp = number of x indices used
! nxhd = first dimension of f >= nx/2
! nyd = second dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
! except for f(1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! aimag(f(1,1)) = real part of mode nx/2,0 and
! aimag(f(1,ny/2+1)) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, j1, j2, k1, k2, ns, ns2, km, kmr
      complex t1, t2
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 80
! inverse fourier transform
      nry = nxhy/ny
! bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
   10 continue
   20 continue
! then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble modes kx = 0, nx/2
      do 70 k = 2, nyh
      if (nxi.eq.1) then
         t1 = f(1,ny2-k)
         f(1,ny2-k) = 0.5*cmplx(aimag(f(1,k) + t1),real(f(1,k) - t1))
         f(1,k) = 0.5*cmplx(real(f(1,k) + t1),aimag(f(1,k) - t1))
      endif
   70 continue
      return
! forward fourier transform
! scramble modes kx = 0, nx/2
   80 do 90 k = 2, nyh
      if (nxi.eq.1) then
         t1 = cmplx(aimag(f(1,ny2-k)),real(f(1,ny2-k)))
         f(1,ny2-k) = conjg(f(1,k) - t1)
         f(1,k) = f(1,k) + t1
      endif
   90 continue
! bit-reverse array elements in y
      nry = nxhy/ny
      do 110 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 110
      do 100 j = nxi, nxt
      t1 = f(j,k1)
      f(j,k1) = f(j,k)
      f(j,k) = t1
  100 continue
  110 continue
! first transform in y
      nry = nxy/ny
      do 150 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 140 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 130 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 120 i = nxi, nxt
      t2 = t1*f(i,j2)
      f(i,j2) = f(i,j1) - t2
      f(i,j1) = f(i,j1) + t2
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2R3X(f,isign,mixup,sct,indx,indy,nyi,nyp,nxhd,nyd,  &
     &nxhyd,nxyhd)
! this subroutine performs the x part of 3 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! y, using complex arithmetic
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, three inverse fourier transforms are performed
! f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, two forward fourier transforms are performed
! f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nyi = initial y index used
! nyp = number of y indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:3,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
! except for f(1:3,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! aimag(f(1:3,1,1)) = real part of mode nx/2,0 and
! aimag(f(1:3,1,ny/2+1)) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nyi, nyp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, nxh, nxhh, nxh2, ny, nxy, nxhy, nyt
      integer nrx, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      real at1, at2, ani
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nyt = nyi + nyp - 1
      if (isign.gt.0) go to 140
! inverse fourier transform
! swap complex components
      do 20 i = nyi, nyt
      do 10 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(real(f(2,j,i)),aimag(f(3,j,i)))
      at2 = aimag(f(2,j,i))
      f(2,j,i) = cmplx(aimag(f(1,j,i)),at1)
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
   10 continue
   20 continue
! bit-reverse array elements in x
      nrx = nxhy/nxh
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   30 continue
   40 continue
! first transform in x
      nrx = nxy/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 60 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 50 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
   50 continue
   60 continue
   70 continue
   80 continue
! unscramble coefficients and normalize
      kmr = nxy/nx
      ani = 1.0/real(2*nx*ny)
      do 110 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 100 k = nyi, nyt
      do 90 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = ani*(t1 + t2)
      f(jj,nxh2-j,k) = ani*conjg(t1 - t2)
   90 continue
  100 continue
  110 continue
      ani = 2.*ani
      do 130 k = nyi, nyt
      do 120 jj = 1, 3
      f(jj,nxhh+1,k) = ani*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = ani*cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),         &
     &                      real(f(jj,1,k)) - aimag(f(jj,1,k)))
  120 continue
  130 continue
      return
! forward fourier transform
! scramble coefficients
  140 kmr = nxy/nx
      do 170 j = 2, nxhh
      t3 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 160 k = nyi, nyt
      do 150 jj = 1, 3
      t2 = conjg(f(jj,nxh2-j,k))
      t1 = f(jj,j,k) + t2
      t2 = (f(jj,j,k) - t2)*t3
      f(jj,j,k) = t1 + t2
      f(jj,nxh2-j,k) = conjg(t1 - t2)
  150 continue
  160 continue
  170 continue
      do 190 k = nyi, nyt
      do 180 jj = 1, 3
      f(jj,nxhh+1,k) = 2.0*conjg(f(jj,nxhh+1,k))
      f(jj,1,k) = cmplx(real(f(jj,1,k)) + aimag(f(jj,1,k)),             &
     &                  real(f(jj,1,k)) - aimag(f(jj,1,k)))
  180 continue
  190 continue
! bit-reverse array elements in x
      nrx = nxhy/nxh
      do 210 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 210
      do 200 k = nyi, nyt
      t1 = f(1,j1,k)
      t2 = f(2,j1,k)
      t3 = f(3,j1,k)
      f(1,j1,k) = f(1,j,k)
      f(2,j1,k) = f(2,j,k)
      f(3,j1,k) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  200 continue
  210 continue
! then transform in x
      nrx = nxy/nxh
      do 250 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 240 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 230 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 220 i = nyi, nyt
      t2 = t1*f(1,j2,i)
      t3 = t1*f(2,j2,i)
      t4 = t1*f(3,j2,i)
      f(1,j2,i) = f(1,j1,i) - t2
      f(2,j2,i) = f(2,j1,i) - t3
      f(3,j2,i) = f(3,j1,i) - t4
      f(1,j1,i) = f(1,j1,i) + t2
      f(2,j1,i) = f(2,j1,i) + t3
      f(3,j1,i) = f(3,j1,i) + t4
  220 continue
  230 continue
  240 continue
  250 continue
! swap complex components
      do 270 i = nyi, nyt
      do 260 j = 1, nxh
      at1 = real(f(3,j,i))
      f(3,j,i) = cmplx(aimag(f(2,j,i)),aimag(f(3,j,i)))
      at2 = real(f(2,j,i))
      f(2,j,i) = cmplx(at1,aimag(f(1,j,i)))
      f(1,j,i) = cmplx(real(f(1,j,i)),at2)
  260 continue
  270 continue
      return
      end
!-----------------------------------------------------------------------
      subroutine FFT2R3Y(f,isign,mixup,sct,indx,indy,nxi,nxp,nxhd,nyd,  &
     &nxhyd,nxyhd)
! this subroutine performs the y part of 3 two dimensional real to
! complex fast fourier transforms, and their inverses, for a subset of
! x, using complex arithmetic
! for isign = (-1,1), input: all, output: f
! for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)
! for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)
! where N = (nx/2)*ny
! indx/indy = exponent which determines length in x/y direction,
! where nx=2**indx, ny=2**indy
! if isign = -1, three inverse fourier transforms are performed
! f(1:3,n,m) = (1/nx*ny)*sum(f(1:3,j,k)*
!       exp(-sqrt(-1)*2pi*n*j/nx)*exp(-sqrt(-1)*2pi*m*k/ny))
! if isign = 1, two forward fourier transforms are performed
! f(1:3,j,k) = sum(f(1:3,n,m)*exp(sqrt(-1)*2pi*n*j/nx)*
!       exp(sqrt(-1)*2pi*m*k/ny))
! mixup = array of bit reversed addresses
! sct = sine/cosine table
! nxi = initial x index used
! nxp = number of x indices used
! nxhd = second dimension of f >= nx/2
! nyd = third dimension of f >= ny
! nxhyd = maximum of (nx/2,ny)
! nxyhd = maximum of (nx,ny)/2
! fourier coefficients are stored as follows:
! f(1:3,j,k) = mode j-1,k-1, where 1 <= j <= nx/2 and 1 <= k <= ny,
! except for f(1:3,1,k) =  mode nx/2,k-1, where ny/2+2 <= k <= ny, and
! aimag(f(1:3,1,1)) = real part of mode nx/2,0 and
! aimag(f(1:3,1,ny/2+1)) = real part of mode nx/2,ny/2
! written by viktor k. decyk, ucla
      implicit none
      integer isign, indx, indy, nxi, nxp, nxhd, nyd, nxhyd, nxyhd
      complex f, sct
      integer mixup
      dimension f(3,nxhd,nyd), mixup(nxhyd), sct(nxyhd)
! local data
      integer indx1, indx1y, nx, ny, nyh, ny2, nxy, nxhy, nxt
      integer nry, i, j, k, l, jj, j1, j2, k1, k2, ns, ns2, km, kmr
      complex t1, t2, t3, t4
      if (isign.eq.0) return
      indx1 = indx - 1
      indx1y = max0(indx1,indy)
      nx = 2**indx
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nxy = max0(nx,ny)
      nxhy = 2**indx1y
      nxt = nxi + nxp - 1
      if (isign.gt.0) go to 90
! inverse fourier transform
      nry = nxhy/ny
! bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
   10 continue
   20 continue
! then transform in y
      nry = nxy/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 40 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = sct(1+kmr*(j-1))
      do 30 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
   30 continue
   40 continue
   50 continue
   60 continue
! unscramble modes kx = 0, nx/2
      do 80 k = 2, nyh
      if (nxi.eq.1) then
         do 70 jj = 1, 3
         t1 = f(jj,1,ny2-k)
         f(jj,1,ny2-k) = 0.5*cmplx(aimag(f(jj,1,k) + t1),               &
     &                             real(f(jj,1,k) - t1))
         f(jj,1,k) = 0.5*cmplx(real(f(jj,1,k) + t1),                    &
     &                         aimag(f(jj,1,k) - t1))
   70    continue
      endif
   80 continue
      return
! forward fourier transform
! scramble modes kx = 0, nx/2
   90 do 110 k = 2, nyh
      if (nxi.eq.1) then
         do 100 jj = 1, 3
         t1 = cmplx(aimag(f(jj,1,ny2-k)),real(f(jj,1,ny2-k)))
         f(jj,1,ny2-k) = conjg(f(jj,1,k) - t1)
         f(jj,1,k) = f(jj,1,k) + t1
  100    continue
      endif
  110 continue
! bit-reverse array elements in y
      nry = nxhy/ny
      do 130 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 130
      do 120 j = nxi, nxt
      t1 = f(1,j,k1)
      t2 = f(2,j,k1)
      t3 = f(3,j,k1)
      f(1,j,k1) = f(1,j,k)
      f(2,j,k1) = f(2,j,k)
      f(3,j,k1) = f(3,j,k)
      f(1,j,k) = t1
      f(2,j,k) = t2
      f(3,j,k) = t3
  120 continue
  130 continue
! first transform in y
      nry = nxy/ny
      do 170 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      t1 = conjg(sct(1+kmr*(j-1)))
      do 140 i = nxi, nxt
      t2 = t1*f(1,i,j2)
      t3 = t1*f(2,i,j2)
      t4 = t1*f(3,i,j2)
      f(1,i,j2) = f(1,i,j1) - t2
      f(2,i,j2) = f(2,i,j1) - t3
      f(3,i,j2) = f(3,i,j1) - t4
      f(1,i,j1) = f(1,i,j1) + t2
      f(2,i,j1) = f(2,i,j1) + t3
      f(3,i,j1) = f(3,i,j1) + t4
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
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
!-----------------------------------------------------------------------
      function randum()
! this is a version of the random number generator dprandom due to
! c. bingham and the yale computer center, producing numbers
! in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      implicit none
      integer isc,i1,r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
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
      randum = (dble(r1) + dble(r2)*asc)*asc
      return
      end
