! Extra Fortran Library for Skeleton 2-1/2D Electromagnetic PIC Code
! contains various analytic Boris pushers
! written by Viktor K. Decyk, UCLA
!-----------------------------------------------------------------------
      subroutine AGBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ek,idimp,nop,nx,ny, &
     &nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, with magnetic field. 
! Using the non-relativistic Analytic Boris Mover,
! assumes constant E, B fields during a time step
! scalar version using guard cells
! 167 flops/particle, 29 loads, 5 stores, 2 divides, 1 sqrt, 1 tangent
! input: all, output: part, ek
! velocity equations used are:
! vx(t+dt/2) = rot(1)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t))) +
!    rot(2)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))) +
!    rot(3)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gx(x(t),y(t))))
! vy(t+dt/2) = rot(4)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t)))) +
!    rot(5)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))) +
!    rot(6)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gy(x(t),y(t))))
! vz(t+dt/2) = rot(7)*(vx(t-dt/2) + (q/m)*gx(x(t),y(t)))) +
!    rot(8)*(vy(t-dt/2) + (q/m)*gy(x(t),y(t)))t) +
!    rot(9)*(vz(t-dt/2) + (q/m)*gz(x(t),y(t))))) + (q/m)*gz(x(t),y(t))))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x,y) = 0.5*flx(x,y)*dt + (fx(x,y)-flx(x,y))*tan(0.5*om*dt)/om
! gy(x,y) = 0.5*fly(x,y)*dt + (fy(x,y)-fly(x,y))*tan(0.5*om*dt)/om
! gz(x,y) = 0.5*flz(x,y)*dt + (fz(x,y)-flz(x,y))*tan(0.5*om*dt)/om
! where flx(x,y) = fpl(x,y)*omx/om**2, fly(x,y) = fpl(x,y)*omy/om**2,
! flz(x,y) = fpl(x,y)*omz/om**2,
! and fpl(x,y) = fx(x,y)*omx+fy(x,y)*omy+fz(x,y)*omz
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
      real dth, omti, epl, x, y
      double precision sum1
      dth = 0.5*dt
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
      x = part(1,j)
      y = part(2,j)
! find interpolation weights
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
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
! normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
! half acceleration
      acx = part(3,j)
      acy = part(4,j)
      acz = part(5,j)
      omxt = acx + dx*dth
      omyt = acy + dy*dth
      omzt = acz + dz*dth
! time-centered kinetic energy
      sum1 = sum1 + (omxt*omxt + omyt*omyt + omzt*omzt)
! normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
      qtmh = dth
! correct the half-acceleration by decomposing E field into components
! parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      epl = dx*ox + dy*oy + dz*oz
      if (omt.gt.0.0) then
         omti = 1.0/omt
         qtmh = omti*tan(dth*omt)
      endif
      epl = epl*omti*omti
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      dx = omxt*dth + (dx - omxt)*qtmh
      dy = omyt*dth + (dy - omyt)*qtmh
      dz = omzt*dth + (dz - omzt)*qtmh
      acx = acx + dx
      acy = acy + dy
      acz = acz + dz
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
      dx = x + dx*dtc
      dy = y + dy*dtc
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
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
      subroutine AGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,nx,&
     &ny,nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Analytic Boris Mover,
! assumes constant E, B fields, and gamma during a time step
! scalar version using guard cells
! 255 flops/particle, 25 loads, 5 stores, 6 divides, 3 sqrts, 1 tangent
! input: all, output: part, ek
! momentum equations used are:
! px(t+dt/2) = rot(1)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
!    rot(2)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
!    rot(3)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gx(x(t),y(t)))
! py(t+dt/2) = rot(4)*(px(t-dt/2) + (q/m)*gx(x(t)),y(t)) +
!    rot(5)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))) +
!    rot(6)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gy(x(t),y(t)))
! pz(t+dt/2) = rot(7)*(px(t-dt/2) + (q/m)*gx(x(t),y(t))) +
!    rot(8)*(py(t-dt/2) + (q/m)*gy(x(t),y(t))t) +
!    rot(9)*(pz(t-dt/2) + (q/m)*gz(x(t),y(t)))) + (q/m)*gz(x(t),y(t)))
! where q/m is charge/mass, and the rotation matrix is given by:
!    rot(1) = (1 - tan(om*dt/2)**2 + 2*((omx/om)*tan(om*dt/2))**2)/norm
!    rot(2) = 2*((omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(3) = 2*(-(omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(4) = 2*(-(omz/om)*tan(om*dt/2) 
!              + (omx*omy/om**2)*tan(om*dt/2)**2)/norm
!    rot(5) = (1 - tan(om*dt/2)**2 + 2*((omy/om)*tan(om*dt/2))**2)/norm
!    rot(6) = 2*((omx/om)*tan(om*dt/2) 
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(7) = 2*((omy/om)*tan(om*dt/2)
!              + (omx*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(8) = 2*(-(omx/om)*tan(om*dt/2)
!              + (omy*omz/om**2)*tan(om*dt/2)**2)/norm
!    rot(9) = (1 - tan(om*dt/2)**2 + 2*((omz/om)*tan(om*dt/2))**2)/norm
! norm = 1 + tan(om*dt/2))**2 and om = sqrt(omx**2 + omy**2 + omz**2)
! gx(x,y) = 0.5*flx(x,y)*dt + (fx(x,y)-flx(x,y))*tan(0.5*om*dt)/om
! gy(x,y) = 0.5*fly(x,y)*dt + (fy(x,y)-fly(x,y))*tan(0.5*om*dt)/om
! gz(x,y) = 0.5*flz(x,y)*dt + (fz(x,y)-flz(x,y))*tan(0.5*om*dt)/om
! where flx(x,y) = fpl(x,y)*omx/om**2, fly(x,y) = fpl(x,y)*omy/om**2,
! flz(x,y) = fpl(x,y)*omz/om**2,
! and fpl(x,y) = fx(x,y)*omx+fy(x,y)*omy+fz(x,y)*omz
! the rotation matrix is determined by:
! omx = (q/m)*bx(x(t),y(t))*gami, omy = (q/m)*by(x(t),y(t))*gami, and
! omz = (q/m)*bz(x(t),y(t))*gami, where
! we approximate analytic gamma function with 4th order taylor series:
! gam = gam(0) + dg0*tau+d2g0*tau**2/2+d3g0*tau**3/6+d4g0*tau**3/24
! where gam(0) = sqrt(1.0 + p2(0)*ci*ci)
! p2(0) = px(t-dt/2)**2+py(t-dt/2)**2+pz(t-dt/2)**2
! dg0 = ug0/(ci*ci) = dgamma/dtau
! d2g0 = u2g0/(ci*ci) = d2gamma/dtau2
! d3g0 = u3g0/(ci*ci) = -omt*omt*(dgamma_perp/dtau)
! d4g0 = u4g0/(ci*ci) = -omt*omt*(d2gamma_perp/dtau2)
! then using the result t = integral of gamma(tau)*dtau, we can
! approximate tau(t) using a one pass Newton method and set gami = tau/t
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
! ek = sum1 + wk/(1.0 + sqrt(1.0 + wk*ci*ci)), where
! wk = px(t-dt/2)**2+py(t-dt/2)**2+pz(t-dt/2)**2+ug0*dt+u2g0*dt*dt/4
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
! weps = taylor series expansion criteria
      integer weps
      parameter (weps=1.0e-3)
      integer j, ii, nn, mm, np, mp
      real sixth, twnth, edgelx, edgely, edgerx, edgery, at1, dth, ci2
      real x, y, dxp, dyp, amx, amy, px, py, pz, dx, dy, dz, ox, oy, oz
      real p2, gam0, wk, acx, acy, acz, omt, omti, epl, omxt, omyt, omzt
      real ug0, dgp0, d2gp0, dgt0, d2gt0, gami, tau, f, fp, qtmg, dtg
      real omti2, ep2, et2, gam, gt0, wt, tn, cs, sn, csc, snc, csd, snd
      real anorm, rot1, rot2, rot3, rot4, rot5, rot6, rot7, rot8, rot9
      real t2, w2, wt2
      double precision sum1
      sixth = 1.0/6.0
      twnth = 1.0/20.0
      dth = 0.5*dt
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
      do 20 j = 1, nop
      x = part(1,j)
      y = part(2,j)
! find interpolation weights
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
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
! normalize electric field
      dx = qbm*dx
      dy = qbm*dy
      dz = qbm*dz
! normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
! read momentum
      px = part(3,j)
      py = part(4,j)
      pz = part(5,j)
! find initial gamma
      p2 = px*px + py*py + pz*pz
      gam0 = sqrt(1.0 + p2*ci2)
      gami = 1.0/gam0
! correct the half-acceleration by decomposing E field into components
! parallel and perpendicular to the B field
      omt = sqrt(ox*ox + oy*oy + oz*oz)
      omti = 0.0
      if (omt.gt.0.0) omti = 1.0/omt
      epl = dx*ox + dy*oy + dz*oz
      omti2 = omti*omti
      epl = epl*omti2
! E parallel
      omxt = epl*ox
      omyt = epl*oy
      omzt = epl*oz
      ep2 = omxt*omxt + omyt*omyt + omzt*omzt
! E perp
      dx = dx - omxt
      dy = dy - omyt
      dz = dz - omzt
      et2 = dx*dx + dy*dy + dz*dz
! dgamma/dtau
      dgp0 = omxt*px + omyt*py + omzt*pz
      dgt0 = dx*px + dy*py + dz*pz
      ug0 = dgp0 + dgt0
! dgamma_paral/dtau
      dgp0 = dgp0*ci2
! dgamma_perp/dtau
      dgt0 = dgt0*ci2
! E dot p cross omega
      acx = py*oz - pz*oy
      acy = pz*ox - px*oz
      acz = px*oy - py*ox
      gt0 = dx*acx + dy*acy + dz*acz
! time-centered kinetic energy 
      wk = p2 + ug0*dt + (ep2 + et2 + gami*gt0)*dth*dth
      sum1 = sum1 + wk/(1.0 + sqrt(1.0 + wk*ci2))
! find tau, gam with Newton method
      gam = gam0
      tau = dt*gami
! start iteration
      do 10 ii = 1, 1
      wt = omt*tau
! d2gamma_paral/dtau2
      d2gp0 = gam*ep2*ci2
! d2gamma_perp/dtau2
      d2gt0 = (gam*et2 + gt0)*ci2
! calculate trigonometric functions
      if (abs(wt).gt.weps) then
         tn = tan(0.5*wt)
         cs = tn*tn
         sn = 1.0/(1.0 + cs)
         cs = (1.0 - cs)*sn
         sn = (tn + tn)*sn
         csc = cs - 1.0
         snc = sn*omti
         csd = csc*omti2
         snd = (snc - tau)*omti2
! second order taylor series approximation for wt
      else
         wt2 = wt*wt
         t2 = tau*tau
         w2 = omt*omt
         csd = -0.5*(1.0 - 0.5*sixth*wt2)*t2
         snd = -sixth*(1.0 - twnth*wt2)*t2*tau
         snc = w2*snd + tau
!        csc = w2*csd
!        cs = csc + 1.0
!        sn = omt*snc
      endif
! calculate time integral
      f = (gam0 + (0.5*dgp0 + sixth*d2gp0*tau)*tau)*tau +               &
     &  - (dgt0*csd + d2gt0*snd) - dt
! calculate gamma
      fp = gam0 + (dgp0 + 0.5*d2gp0*tau)*tau + (dgt0*snc - d2gt0*csd)
! update tau
      tau = tau - f/fp
! update gam
      gam = dt/tau
   10 continue
      gami = 1.0/gam
! set proper time
      qtmg = dth*gami
      if (omt.gt.0.0) qtmg = omti*tan(qtmg*omt)
      at1 = qtmg/gami
! modified half acceleration
      dx = omxt*dth + dx*at1
      dy = omyt*dth + dy*at1
      dz = omzt*dth + dz*at1
      acx = px + dx
      acy = py + dy
      acz = pz + dz
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
! new momentum
      px = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      py = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      pz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
! update inverse gamma
      p2 = px*px + py*py + pz*pz
      dtg = dtc/sqrt(1.0 + p2*ci2)
! new position
      dx = x + px*dtg
      dy = y + py*dtg
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   20 continue
! normalize kinetic energy
      ek = ek + sum1
      return
      end
!-----------------------------------------------------------------------
      subroutine EAGRBPUSH23L(part,fxy,bxy,qbm,dt,dtc,ci,ek,idimp,nop,  &
     &nx,ny,nxv,nyv,ipbc)
! for 2-1/2d code, this subroutine updates particle co-ordinates and
! velocities using leap-frog scheme in time and first-order linear
! interpolation in space, for relativistic particles with magnetic field
! Using the Exact Analytic mover, a variant of the algorithm developed
! by Fei Li et al., where E and B fields are constant and gamma varies
! during a time step.
! assumes constant E, B fields, and gamma during a time step
! scalar version using guard cells
! 278 FLOPs, 11 divides, 6 sqrts per particle
! plus 48 FLOPs, 3 divides, 1 tan, 1 tanh per Newton iteration/particle
! input: all, output: part, ek
! The particle momenta are advanced in natural coordinates, defined as:
! up is in the direction ep = E parallel to the B field,
! ut is in the direction of et = E perpendicular to the B field,
! ud is in the ExB direction.  Momenta are transformed from cartesian
! coordinates to natural coordinates before the advance, then
! transformed back to cartesian coordinates.
! momenta equations used are:
! up = up + ep*dt
! ut = ut + etp*dt + (om/w)*(F*(cos(w*tau)-1.0) + G*sin(w*tau))
! ud = ved*gam(tau) + F*sin(w*tau) - H*cos(w*tau), where
! F = (w/omega)*ut(0) - om2t*(omega*w)*ci2*ep*up(0)
! G = gam0*(et - al2*om2t) - ud(0)*omega
! H = om2t*omega*gam0 - ud(0)
! om = sqrt(qbm*omx)**2 + (qbm*by)**2 + (qbm*bz)**2)
! etp = om2t*al*al, ved = om*om2t, om2t = et/(om*om + al*al
! the gyration frequency w and al are given by:
! w = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))+om*om-e*e))
! al = sqrt(0.5*(sqrt(om*om-e*e)**2 + 4.0*(ci*ep*om)**2))-om*om+e*e))
! where e = sqrt((qbm*ex)**2+(qbm*ey)**2+(qbm*ez)**2)
! gam(tau) is the relativistic factor given by:
! gam = A*sin(w*tau) - B*cos(w*tau) + C*sinh(al*tau) + D*cos(al*tau),
! where the constants A, B, C, D are defined in terms of initial
! derivatives of the gamma function, as follows:
!  w*A = (w*w*dg0 - om*om*dgp0)/(w*w + al*al)
!  w*B = (w*w*d2g0 - om*om*d2gp0)/(w*w + al*al)
!  al*C = (al*al*dg0 + om*om*dgp0)/(w*w + al*al)
!  al*C = (al*al*d2g0 + om*om*d2gp0)/(w*w + al*al)
! and where the initial derivatives of the gamma function are:
!  dgp0 = ci*ci*ep*up(0)
!  d2gp0 = ci*ci*ep*ep*gam(0)
!  dg0 = dgp0 + ci*ci*et*ut(0)
!  d2g0 = d2gp0 + ci*ci*et*(et*gam(0) - ud(0)*om)
! the proper time tau can be determined from the relativistic factor
! using the equation t = integral gam(tau)*dtau.  This involves finding
! zeros of the function f via an iteration scheme, where
! f = gam(0)*tau - dt - (A*(cos(w*tau)-1.0) + B*(sin(w*tau)-w*tau))/w
!                + (C*(cosh(al*tau)-1.0) + C*(sinh(al*tau)-w*tau))/al
! once we know tau(dt), we can evaluate the momenta up, ut, and ud
! and transform back to cartesian coordinates
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
! kinetic energy/mass at time t is also calculated, using the average
! value of energies at the beginning and end of the time step:
! ek = 0.5*((px(t-dt/2)*2+py(t-dt/2)**2+pz(t-dt/2)**2)/(1.0+gam(0))
!    +      (px(t+dt/2)*2+py(t+dt/2)**2+pz(t+dt/2)**2)/(1.0+gam(tau)))
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
      integer imax
      real half, sixth, s1, s2, s3, s4, c1, c2, c3, c4, weps, teps, deps
      real erps
! imax = maximum iteration count
      parameter (imax=25)
      parameter (half=0.5,sixth=1.0/6.0)
      parameter (s1=1.0/6.0,s2=1.0/20.0,s3=1.0/42.0,s4=1.0/72.0)
      parameter (c1=0.5,c2=1.0/12.0,c3=1.0/30.0,c4=1.0/56.0)
! weps = taylor series expansion criteria
      parameter (weps=1.0e-1)
! deps = taylor series expansion criteria for differences
      parameter (deps=1.0e-3)
      integer i, j, nn, mm, np, mp, ii
      real edgelx, edgely, edgerx, edgery
      real dxp, dyp, amx, amy, ci2, small, prec, vresult, erm
      real dx, dy, dz, ox, oy, oz, tx, ty, tz, om2, e2, omega, di, ep
      real et, px, py, pz, up, ut, ud, w2, wl2, al2, w, al, wi, ali, wi2
      real ali2, om2t, ws, p2, gam0, dgp0, dgt0, d2gp0, d2gt0, dg0, d2g0
      real wt2, ea, eb, ec, ed, fa, fb, tau, cs, sn, csh, snh, fpp, fp
      real ft, f, t2, t3, wt, wd, tn, tn2, tn2i, csd, snd, cshd, snhd
      real csc, snc, cshc, snhc, fc, fpi, gam, gami, wk, dtg, x, y
      double precision sum1
      double precision dt1, err
      data small /1.0e-12/
      prec = 1.0 + small
! teps = tau precision; erps = orthogonality precision
! detect autodouble precision
      if (vresult(prec).gt.1.0) then
         teps = 1.0e-14
         erps = 1.0e-10
! default single precision
      else
         teps = 4.0e-8
         erps = 4.0e-3
      endif
      ci2 = ci*ci
      erm = 0.0
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
      do 70 j = 1, nop
      x = part(1,j)
      y = part(2,j)
! find interpolation weights
      nn = x
      mm = y
      dxp = x - real(nn)
      dyp = y - real(mm)
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
! normalize magnetic field
      ox = qbm*ox
      oy = qbm*oy
      oz = qbm*oz
! normalize electric field
      tx = qbm*dx
      ty = qbm*dy
      tz = qbm*dz
! create direction cosines to translate from/to cartesian coordinates
! first find the direction along the magnetic field B
      om2 = ox*ox + oy*oy + oz*oz
      omega = sqrt(om2)
      if (om2.gt.0.0) then
         di = 1.0/omega
         ox = ox*di
         oy = oy*di
         oz = oz*di
! if omega = 0, then use the y direction
      else
         ox = 0.0
         oy = 1.0
         oz = 0.0
      endif
! then find the direction along the electric field Et perpendicular to B
      dt1 = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
      tx = tx - dt1*ox
      ty = ty - dt1*oy
      tz = tz - dt1*oz
      ep = dt1
      e2 = tx*tx + ty*ty + tz*tz
      et = sqrt(e2)
      if (et > 0.0) then
         di = 1.0/et
         tx = tx*di
         ty = ty*di
         tz = tz*di
! then find the direction along Et x B
         dx = ty*oz - tz*oy
         dy = tz*ox - tx*oz
         dz = tx*oy - ty*ox
! check for roundoff error
         err = dble(tx*ox) + dble(ty*oy) + dble(tz*oz)
         if (err > erps) then
            write (*,*) 'Error: Et not normal to omega=',et,err
            et = 0.0
         else
            err = dble(tx*dx) + dble(ty*dy) + dble(tz*dz)
            if (err > erps) then
               write (*,*) 'Error: Et d=',et,err
            endif
         endif
      endif
! special case Et = 0, or round off error detected
      if (et==0.0) then      
! first find direction with smallest component of B
         ii = 1
         dx = abs(ox)
         dy = abs(oy)
         dz = abs(oz)
         di = dx
         if (dy <= dx) then
            ii = 2
            di = dy
            if ((dx.gt.dy).and.(dy.eq.dz)) ii = 3
         endif
         if (dz.lt.di) ii = 3 
! then the cross product of that direction with B
         if (ii.eq.1) then
            dz = 1.0/sqrt(oy*oy + oz*oz)
            dx = 0.0
            dy = -oz*dz
            dz = oy*dz
         else if (ii.eq.2) then
            dz = 1.0/sqrt(ox*ox + oz*oz)
            dx = oz*dz
            dy = 0.0
            dz = -ox*dz
         else if (ii.eq.3) then
            dz = 1.0/sqrt(ox*ox + oy*oy)
            dx = -oy*dz
            dy = ox*dz
            dz = 0.0
         endif
! then find the direction along minus d x B
         tx = dz*oy - dy*oz
         ty = dx*oz - dz*ox
         tz = dy*ox - dx*oy
      endif
! calculate frequencies
      e2 = (e2 + ep*ep)*ci2
      w2 = om2 - e2
      wl2 = sqrt(w2*w2 + 4.0*ci2*(ep*omega)**2)
      al2 = 0.5*(wl2 - w2)
      w2 = 0.5*(wl2 + w2)
      w = sqrt(w2)
      al = sqrt(al2)
      wi = 0.0
      if (w > 0.0) wi = 1.0/w
      ali = 0.0
      if (al > 0.0) ali = 1.0/al
      wi2 = wi*wi
      ali2 = ali*ali
! calculate weights
      om2t = om2 + al2
      if (om2t > 0.0) om2t = et/om2t
      ws = 0.0
      if (omega /= 0.0) ws = w/omega
! translate momenta from cartesian coordinates
      px = part(3,j)
      py = part(4,j)
      pz = part(5,j)
      up = px*ox + py*oy + pz*oz
      ut = px*tx + py*ty + pz*tz
      ud = px*dx + py*dy + pz*dz
! find initial gamma
      p2 = up*up + ut*ut + ud*ud
      gam0 = sqrt(1.0 + p2*ci2)
! calculate initial kinetic energy
      wk = p2/(1.0 + gam0)
! partial derivatives of gamma
      dgp0 = ci2*ep*up
      dgt0 = ci2*et*ut
      d2gp0 = ci2*ep*ep*gam0
      d2gt0 = ci2*et*(et*gam0 - ud*omega)
      dg0 = dgp0 + dgt0
      d2g0 = d2gp0 + d2gt0
! calculate trigonometeric and hyperbolic coefficients
      wt2 = 0.0
      if (wl2 > 0.0) wt2 = 1.0/wl2
      ea = wt2*(w2*dg0 - om2*dgp0)
      eb = wt2*(w2*d2g0 - om2*d2gp0)
      ec = wt2*(al2*dg0 + om2*dgp0)
      ed = wt2*(al2*d2g0 + om2*d2gp0)
      if (wl2==0.0) then
         ea = dgt0
         eb = d2gt0
      endif
      fa = ws*ut - om2t*omega*wi*ci2*ep*up
      fb = gam0*(et - al2*om2t) - ud*omega
      fc = om2t*omega*gam0 - ud
! zeroth order guess for tau
      tau = dt/gam0
!
! iteration loop for finding tau(t)
      cs = 1.0
      sn = 0.0
      csh = 1.0
      snh = 0.0
      fpp = 0.0
      fp = 0.0
      ft = -tau
      do 50 ii = 1, imax
      t2 = tau*tau; t3 = t2*tau
! calculate trigonometric functions
      wt = w*tau
      wd = w*ft
      if (abs(wt) > weps) then
         if (abs(wd) > deps) then
            tn = tan(0.5*wt)
            tn2 = tn*tn
            tn2i = 1.0/(1.0 + tn2)
            cs = (1.0 - tn2)*tn2i
            sn = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 - s1*wt2)
            wt2 = 1.0 - c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 - s1*wt2*(1.0 - s2*wt2))
!           wt2 = 1.0 - c1*wt2*(1.0 - c2*wt2)
            f = cs*wt2 - sn*wd
            sn = sn*wt2 + cs*wd
            cs = f
         endif
! calculate special functions
         csc = cs - 1.0
         snc = sn*wi
         csd = csc*wi2
         snd = (snc - tau)*wi2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2))*t2
         snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2))*t3
! fifth order taylor series approximation for wt
!        csd = -c1*(1.0 - c2*wt2*(1.0 - c3*wt2*(1.0 - c4*wt2)))*t2
!        snd = -s1*(1.0 - s2*wt2*(1.0 - s3*wt2*(1.0 - s4*wt2)))*t3
         csc = w2*csd
         snc = tau + w2*snd
         cs = csc + 1.0
         sn = w*snc
      endif
! calculate hyperbolic functions
      wt = al*tau
      wd = al*ft
      if (abs(wt) > weps) then
         if (abs(wd) > deps) then
            tn = tanh(0.5*wt)
            tn2 = tn*tn;
            tn2i = 1.0/(1.0 - tn2)
            csh = (1.0 + tn2)*tn2i
            snh = (tn + tn)*tn2i
! second order taylor series approximation for wd
         else
            wt2 = wd*wd
            wd = -wd*(1.0 + s1*wt2)
            wt2 = 1.0 + c1*wt2
! third order taylor series approximation for wd
!           wd = -wd*(1.0 + s1*wt2*(1.0 + s2*wt2))
!           wt2 = 1.0d+ c1*wt2*(1.0 + c2*wt2)
            f = csh*wt2 + snh*wd
            snh = snh*wt2 + csh*wd
            csh = f
         endif
! calculate special functions
         cshc = csh - 1.0
         snhc = snh*ali
         cshd = cshc*ali2
         snhd = (snhc - tau)*ali2
! fourth order taylor series approximation for wt
      else
         wt2 = wt*wt
         cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2))*t2
         snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2))*t3
! fifth order taylor series approximation for wt
!        cshd = c1*(1.0 + c2*wt2*(1.0 + c3*wt2*(1.0 + c4*wt2)))*t2
!        snhd = s1*(1.0 + s2*wt2*(1.0 + s3*wt2*(1.0 + s4*wt2)))*t3
         cshc = al2*cshd
         snhc = tau + al2*snhd
         csh = cshc + 1.0
         snh = al*snhc
      endif
! gam = gamma(tau)
      gam = gam0 + (ea*snc - eb*csd) + (ec*snhc + ed*cshd)
      fpi = 1.0/gam
! calculate time expression whose root we seek
      f = gam0*tau - (ea*csd + eb*snd) + (ec*cshd + ed*snhd) - dt
! newton's quadratic method
      ft = f*fpi
! either add Halley's optional cubic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
!     ft = ft/(1.0d0 - 0.5d0*ft*fpp*fpi)
! or add Householder's optional quartic correction
! fpp = dgamma/dtau
!     fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     wt2 = ft*fpp*fpi
!     ft = ft*(1.0d0 - 0.5d0*wt2)/(1.0d0 - wt2 + sixth*(ft*ft)*fp*fpi)
! update tau: ft = -delta tau
      tau = tau - ft
      if (abs(f) < teps) go to 60
   50 continue

!
! convergence failure
   60 continue
      if (ii.gt.imax) then
         write (*,*) i,'tau error:f,ft=',f,ft
      endif
! update gamma
! fpp = dgamma/dt
      if (fpp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! first order taylor series
      gam = gam - fpp*ft
! second order taylor series
! fpp = dgamma/dtau
!     if (fp.eq.0.0) fpp = (ea*cs + eb*snc) + (ec*csh + ed*snhc)
! fp = d2gamma/dtau2
!     if (fp.eq.0.0) fp = (eb*cs - ea*w*sn) + (ec*al*snh + ed*csh)
!     gam = gam - (fpp - half*fp*ft)*ft
!
! update sine/cosine functions
! first order taylor series
      wd = w*ft
      f = sn*wd
      csc = csc + f
      snc = snc - cs*ft
      f = cs + f
      sn = sn - cs*wd
      cs = f
! second order taylor series
!     wd = w*ft
!     wt2 = wd*wd
!     f = -ft*(1.0d0 - sixth*wt2)
!     wd = half*wt2
!     wt2 = 1.0d0 - wd
!     sc = csc*wt2 - wd
!     wd = w*f
!     sc = csc - sn*wd
!     snc = snc*wt2 + cs*f
!     f = cs*wt2 - sn*wd
!     sn = sn*wt2 + cs*wd
!     cs = f
! update momenta
      up = up + ep*dt
      ut = ut + om2t*al2*dt + (omega*wi*fa*csc + fb*snc)
      ud = om2t*omega*gam + fa*sn - fc*cs
! calculate inverse gamma and average kinetic energy
      gami = 1.0/gam
      p2 = up*up + ut*ut + ud*ud
      sum1 = sum1 + dble(0.5*(wk + gami*p2/(1.0 + gami)))
! sanity check
!     erm = max(erm,abs(gam-sqrt(1.0d0 + p2*ci2)))
! translate momenta to cartesian coordinates
      px = up*ox + ut*tx + ud*dx
      py = up*oy + ut*ty + ud*dy
      pz = up*oz + ut*tz + ud*dz
      part(3,j) = px
      part(4,j) = py
      part(5,j) = pz
! new position
      dtg = dtc/sqrt(1.0 + p2*ci2)
      dx = x + px*dtg
      dy = y + py*dtg
! periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
! reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(3,j) = -part(3,j)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = y
            part(4,j) = -part(4,j)
         endif
! mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = x
            part(3,j) = -part(3,j)
         endif
         if (dy.lt.edgely) dy = dy + edgery
         if (dy.ge.edgery) dy = dy - edgery
      endif
! set new position
      part(1,j) = dx
      part(2,j) = dy
   70 continue
! normalize kinetic energy
      ek = ek + sum1
! sanity check
!     write (*,*) 'gamma sanity check=',erm
      return
      end
!-----------------------------------------------------------------------
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
      return
      end

