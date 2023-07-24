!-----------------------------------------------------------------------
! Skeleton 2-1/2D Electromagnetic PIC code
! written by Viktor K. Decyk, UCLA
      program bpic2
      use wbpush2
      use wabpush2
      implicit none
! indx/indy = exponent which determines grid points in x/y direction:
! nx = 2**indx, ny = 2**indy.
      integer, parameter :: indx =   9, indy =   9
! npx/npy = number of electrons distributed in x/y direction.
      integer, parameter :: npx =  3072, npy =   3072
! ndim = number of velocity coordinates = 3
      integer, parameter :: ndim = 3
! tend = time at end of simulation, in units of plasma frequency.
! dt = time interval between successive calculations.
! qme = charge on electron, in units of e.
      real, parameter :: tend = 10.0, dt = 0.04, qme = -1.0
! vtx/vty = thermal velocity of electrons in x/y direction
! vx0/vy0 = drift velocity of electrons in x/y direction.
      real, parameter :: vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0
! vtx/vz0 = thermal/drift velocity of electrons in z direction
      real, parameter :: vtz = 1.0, vz0 = 0.0
! ax/ay = smoothed particle size in x/y direction
! ci = reciprocal of velocity of light.
      real :: ax = .912871, ay = .912871, ci = 0.1
! idimp = number of particle coordinates = 5
! sortime = number of time steps between standard electron sorting
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: idimp = 5, sortime = 50, relativity = 1
! kpush = push type = (1,2,3) = (Boris,Analytic Boris,Exact Analytic)
      integer :: kpush = 1
!
! wke/we = particle kinetic/electrostatic field energy
! wf/wm/wt = magnetic field/transverse electric field/total energy
      real :: wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0
! declare scalars for standard code
      integer :: n
      integer :: np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxyh, nxhy
      integer :: ny1, ntime, nloop, isign
      real :: qbme, affp, dth
!
! declare arrays for standard code:
! part, part2 = particle arrays
      real, dimension(:,:), pointer :: part, part2, tpart
! qe = electron charge density with guard cells
      real, dimension(:,:), pointer :: qe
! cue = electron current density with guard cells
! fxyze/bxyze = smoothed electric/magnetic field with guard cells
      real, dimension(:,:,:), pointer :: cue, fxyze, bxyze
! exyz/bxyz = transverse electric/magnetic field in fourier space
      complex, dimension(:,:,:), pointer :: exyz, bxyz
! ffc = form factor array for poisson solver
      complex, dimension(:,:), pointer :: ffc
! mixup = bit reverse table for FFT
      integer, dimension(:), pointer :: mixup
! sct = sine/cosine table for FFT
      complex, dimension(:), pointer :: sct
! npicy = scratch array for reordering particles
      integer, dimension(:), pointer :: npicy
!
! declare and initialize timing data
      real :: time
      integer, dimension(4) :: itime
      real :: tdpost = 0.0, tguard = 0.0, tfft = 0.0, tfield = 0.0
      real :: tdjpost = 0.0, tpush = 0.0, tsort = 0.0
      double precision :: dtime
!
! initialize scalars for standard code
! np = total number of particles in simulation
! nx/ny = number of grid points in x/y direction
      np = npx*npy; nx = 2**indx; ny = 2**indy
      nxh = nx/2; nyh = max(1,ny/2)
      nxe = nx + 2; nye = ny + 1; nxeh = nxe/2
      nxyh = max(nx,ny)/2; nxhy = max(nxh,ny); ny1 = ny + 1
! nloop = number of time steps in simulation
! ntime = current time step
      nloop = tend/dt + .0001; ntime = 0
      qbme = qme
      affp = real(nx*ny)/real(np)
      dth = 0.0
!
! allocate data for standard code
      allocate(part(idimp,np))
      if (sortime > 0) allocate(part2(idimp,np))
      allocate(qe(nxe,nye),fxyze(ndim,nxe,nye))
      allocate(cue(ndim,nxe,nye),bxyze(ndim,nxe,nye))
      allocate(exyz(ndim,nxeh,nye),bxyz(ndim,nxeh,nye))
      allocate(ffc(nxh,nyh),mixup(nxhy),sct(nxyh))
      allocate(npicy(ny1))
!
! prepare fft tables
      call fft2_init(mixup,sct,indx,indy)
! calculate form factors
      call pois2_init(ffc,ax,ay,affp,nx,ny)
! initialize electrons
      call wdistr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,nx,ny)
!
! initialize transverse electromagnetic fields
      exyz = cmplx(0.0,0.0)
      bxyz = cmplx(0.0,0.0)
!
      if (dt > 0.45*ci) then
         write (*,*) 'Warning: Courant condition may be exceeded!'
      endif
!
! * * * start main iteration loop * * *
!
      do n = 1, nloop
      ntime = n - 1
!     write (*,*) 'ntime = ', ntime
!
! deposit current with standard procedure: updates part, cue
      call dtimer(dtime,itime,-1)
      cue = 0.0
      if (relativity==1) then
         call wgrjpost2(part,cue,qme,dth,ci,nx,ny)
      else
         call wgjpost2(part,cue,qme,dth,nx,ny)
      endif
      call dtimer(dtime,itime,1)
      tdjpost = tdjpost + real(dtime)
!
! deposit charge with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      qe = 0.0
      call wgpost2(part,qe,qme)
      call dtimer(dtime,itime,1)
      tdpost = tdpost + real(dtime)
!
! add guard cells with standard procedure: updates cue, qe
      call dtimer(dtime,itime,-1)
      call wacguard2(cue,nx,ny)
      call waguard2(qe,nx,ny)
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
!
! transform charge to fourier space with standard procedure: updates qe
      call dtimer(dtime,itime,-1)
      isign = -1
      call wfft2r(qe,isign,mixup,sct,indx,indy)
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
!
! transform current to fourier space with standard procedure: update cue
      call dtimer(dtime,itime,-1)
      isign = -1
      call w3fft2r(cue,isign,mixup,sct,indx,indy)
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
!
! take transverse part of current with standard procedure: updates cue
      call dtimer(dtime,itime,-1)
      call wcuperp2(cue,nx,ny)
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
!
! calculate electromagnetic fields in fourier space with standard
! procedure: updates exyz, bxyz, wf, wm
      call dtimer(dtime,itime,-1)
      if (ntime==0) then
         call wibpois2(cue,bxyz,ffc,ci,wm,nx,ny)
         wf = 0.0
         dth = 0.5*dt
      else
         call wmaxwel2(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny)
      endif
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
!
! calculate force/charge in fourier space with standard procedure:
! updates fxyze, we
      call dtimer(dtime,itime,-1)
      isign = -1
      call wpois2(qe,fxyze,ffc,we,nx,ny)
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
!
! add longitudinal and transverse electric fields with standard
! procedure: updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call wemfield2(fxyze,exyz,ffc,isign,nx,ny)
! copy magnetic field with standard procedure: updates bxyze
      isign = -1
      call wemfield2(bxyze,bxyz,ffc,isign,nx,ny)
      call dtimer(dtime,itime,1)
      tfield = tfield + real(dtime)
!
! transform electric force to real space with standard procedure:
! updates fxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call w3fft2r(fxyze,isign,mixup,sct,indx,indy)
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
!
! transform magnetic force to real space with standard procedure:
! updates bxyze
      call dtimer(dtime,itime,-1)
      isign = 1
      call w3fft2r(bxyze,isign,mixup,sct,indx,indy)
      call dtimer(dtime,itime,1)
      tfft = tfft + real(dtime)
!
! copy guard cells with standard procedure: updates fxyze, bxyze
      call dtimer(dtime,itime,-1)
      call wcguard2(fxyze,nx,ny)
      call wcguard2(bxyze,nx,ny)
      call dtimer(dtime,itime,1)
      tguard = tguard + real(dtime)
!
! push particles with various procedures: updates part, wke
      wke = 0.0
      call dtimer(dtime,itime,-1)
! relativistic pushers
      if (relativity==1) then
! classic Boris pusher
         if (kpush==1) then
            call wgrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,ny)
         else if (kpush==2) then
! analytic Boris pusher
            call wagrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,ny)
! exact analytic pusher
         else if (kpush==3) then
            call weagrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,ny)
         endif
! non-relativistic pushers
      else
! classic Boris pusher
         if (kpush==1) then
            call wgbpush2(part,fxyze,bxyze,qbme,dt,dth,wke,nx,ny)
! analytic Boris pusher
         else if ((kpush==2).or.(kpush==3)) then
            call wagbpush2(part,fxyze,bxyze,qbme,dt,dth,wke,nx,ny)
         endif
      endif
      call dtimer(dtime,itime,1)
      tpush = tpush + real(dtime)
!
! sort particles by cell for standard procedure
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
            call dtimer(dtime,itime,-1)
            call wdsortp2yl(part,part2,npicy)
! exchange pointers
            tpart => part
            part => part2
            part2 => tpart
            call dtimer(dtime,itime,1)
            tsort = tsort + real(dtime)
         endif
      endif
!
      if (ntime==0) then
         wt = we + wf + wm
         write (*,*) 'Initial Total Field, Kinetic and Total Energies:'
         write (*,'(3e14.7)') wt, wke, wke + wt
         write (*,*) 'Initial Electrostatic, Transverse Electric and Mag&
     &netic Field Energies:'
         write (*,'(3e14.7)') we, wf, wm
      endif
!
      enddo
      ntime = ntime + 1
!
! * * * end main iteration loop * * *
!
      write (*,*) 'ntime, relativity, kpush = ', ntime,relativity,kpush
      wt = we + wf + wm
      write (*,*) 'Final Total Field, Kinetic and Total Energies:'
      write (*,'(3e14.7)') wt, wke, wke + wt
      write (*,*) 'Final Electrostatic, Transverse Electric and Magnetic&
     & Field Energies:'
      write (*,'(3e14.7)') we, wf, wm
!
      write (*,*)
      write (*,*) 'deposit time = ', tdpost
      write (*,*) 'current deposit time = ', tdjpost
      tdpost = tdpost + tdjpost
      write (*,*) 'total deposit time = ', tdpost
      write (*,*) 'guard time = ', tguard
      write (*,*) 'solver time = ', tfield
      write (*,*) 'fft time = ', tfft
      write (*,*) 'push time = ', tpush
      write (*,*) 'sort time = ', tsort
      tfield = tfield + tguard + tfft
      write (*,*) 'total solver time = ', tfield
      time = tdpost + tpush + tsort
      write (*,*) 'total particle time = ', time
      wt = time + tfield
      write (*,*) 'total time = ', wt
      write (*,*)
!
      wt = 1.0e+09/(real(nloop)*real(np))
      write (*,*) 'Push Time (nsec) = ', tpush*wt
      write (*,*) 'Deposit Time (nsec) = ', tdpost*wt
      write (*,*) 'Sort Time (nsec) = ', tsort*wt
      write (*,*) 'Total Particle Time (nsec) = ', time*wt
!
      stop
      end program
