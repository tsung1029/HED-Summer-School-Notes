#-----------------------------------------------------------------------
from __future__ import print_function
# Skeleton 2-1/2D Electromagnetic PIC code
# written by Viktor K. Decyk, UCLA
import math
import numpy

from bpic2lib import *
from dtimer import *

if (wbpush2.fprecision()==0):
   float_type = numpy.float32
   double_type = numpy.float64
   complex_type = numpy.complex64
   print("single precision floating point library required")
elif (wbpush2.fprecision()==1):
   float_type = numpy.float64
   double_type = numpy.float64
   complex_type = numpy.complex128

#-----------------------------------------------------------------------
def main():
# indx/indy = exponent which determines grid points in x/y direction:
# nx = 2**indx, ny = 2**indy.
   indx =   9; indy =   9
# npx/npy = number of electrons distributed in x/y direction.
   npx =  3072; npy =   3072
# ndim = number of velocity coordinates = 3
   ndim = 3
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
# qme = charge on electron, in units of e.
   tend = 10.0; dt = 0.04; qme = -1.0
# vtx/vty = thermal velocity of electrons in x/y direction
# vx0/vy0 = drift velocity of electrons in x/y direction.
   vtx = 1.0; vty = 1.0; vx0 = 0.0; vy0 = 0.0
# vtx/vz0 = thermal/drift velocity of electrons in z direction
   vtz = 1.0; vz0 = 0.0
# ax/ay = smoothed particle size in x/y direction
# ci = reciprocal of velocity of light.
   ax = .912871; ay = .912871; ci = 0.1
# idimp = number of particle coordinates = 5
# sortime = number of time steps between standard electron sorting
# relativity = (no,yes) = (0,1) = relativity is used
   idimp = 5; sortime = 50; relativity = 1
# kpush = push type = (1,2,3) = (Boris,Analytic Boris,Exact Analytic)
   kpush = 1

# wke/we = particle kinetic/electrostatic field energy
# wf/wm/wt = magnetic field/transverse electric field/total energy
   wke = numpy.zeros((1),float_type)
   we = numpy.zeros((1),float_type)
   wf = numpy.zeros((1),float_type)
   wm = numpy.zeros((1),float_type)
   wt = numpy.zeros((1),float_type)

# declare and initialize timing data
   tinit = 0.0; tloop = 0.0
   itime = numpy.empty((4),numpy.int32)
   dtime = numpy.empty((1),double_type)
   tdpost = numpy.zeros((1),float_type)
   tguard = numpy.zeros((1),float_type)
   tfft = numpy.zeros((1),float_type)
   tfield = numpy.zeros((1),float_type)
   tdjpost = numpy.zeros((1),float_type)
   tpush = numpy.zeros((1),float_type)
   tsort = numpy.zeros((1),float_type)
   ws = numpy.zeros((1),float_type)

# initialize scalars for standard code
# np = total number of particles in simulation
# nx/ny = number of grid points in x/y direction
   np = npx*npy
   nx = int(math.pow(2,indx)); ny = int(math.pow(2,indy))
   nxh = int(nx/2); nyh = max(1,int(ny/2))
   nxe = nx + 2; nye = ny + 1; nxeh = int(nxe/2)
   nxyh = int(max(nx,ny)/2); nxhy = max(nxh,ny); ny1 = ny + 1
# nloop = number of time steps in simulation
# ntime = current time step
   nloop = int(tend/dt + .0001); ntime = 0
   qbme = qme
   affp = float(nx*ny)/float(np)
   dth = 0.0

# allocate data for standard code
# part, part2 = particle arrays
   part = numpy.empty((idimp,np),float_type,'F')
   if (sortime > 0):
      part2 = numpy.empty((idimp,np),float_type,'F')
# qe = electron charge density with guard cells
   qe = numpy.empty((nxe,nye),float_type,'F')
# cue = electron current density with guard cells
   cue = numpy.empty((ndim,nxe,nye),float_type,'F')
# fxyze/bxyze = smoothed electric/magnetic field with guard cells
   fxyze = numpy.empty((ndim,nxe,nye),float_type,'F')
   bxyze = numpy.empty((ndim,nxe,nye),float_type,'F')
# exyz/bxyz = transverse electric/magnetic field in fourier space
   exyz = numpy.empty((ndim,nxeh,nye),complex_type,'F')
   bxyz = numpy.empty((ndim,nxeh,nye),complex_type,'F')
# ffc = form factor array for poisson solver
   ffc = numpy.empty((nxh,nyh),complex_type,'F')
# mixup = bit reverse table for FFT
   mixup = numpy.empty((nxhy),numpy.int32)
# sct = sine/cosine table for FFT
   sct = numpy.empty((nxyh),complex_type,'F')
# npicy = scratch array for reordering particles
   npicy = numpy.empty((ny1),numpy.int32)

# prepare fft tables
   wbpush2.fft2_init(mixup,sct,indx,indy)
# calculate form factors
   wbpush2.pois2_init(ffc,ax,ay,affp,nx,ny)
# initialize electrons
   wbpush2.wdistr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,nx,ny)

# initialize transverse electromagnetic fields
   exyz.fill(complex(0.0,0.0))
   bxyz.fill(complex(0.0,0.0))

   if (dt > 0.45*ci):
      print("Warning: Courant condition may be exceeded!")

# * * * start main iteration loop * * *

   for ntime in range(0,nloop):
#     print("ntime = ", ntime)

# deposit current with standard procedure: updates part, cue
      dtimer(dtime,itime,-1)
      cue.fill(0.0)
      if (relativity==1):
         wbpush2.wgrjpost2(part,cue,qme,dth,ci,nx,ny)
      else:
         wbpush2.wgjpost2(part,cue,qme,dth,nx,ny)
      dtimer(dtime,itime,1)
      tdjpost[0] += float(dtime)

# deposit charge with standard procedure: updates qe
      dtimer(dtime,itime,-1)
      qe.fill(0.0)
      wbpush2.wgpost2(part,qe,qme)
      dtimer(dtime,itime,1)
      tdpost[0] += float(dtime)

# add guard cells with standard procedure: updates cue, qe
      dtimer(dtime,itime,-1)
      wbpush2.wacguard2(cue,nx,ny)
      wbpush2.waguard2(qe,nx,ny)
      dtimer(dtime,itime,1)
      tguard[0] += float(dtime)

# transform charge to fourier space with standard procedure: updates qe
      dtimer(dtime,itime,-1)
      isign = -1
      wbpush2.wfft2r(qe,isign,mixup,sct,indx,indy)
      dtimer(dtime,itime,1)
      tfft[0] += float(dtime)

# transform current to fourier space with standard procedure: update cue
      dtimer(dtime,itime,-1)
      isign = -1
      wbpush2.w3fft2r(cue,isign,mixup,sct,indx,indy)
      dtimer(dtime,itime,1)
      tfft[0] += float(dtime)

# take transverse part of current with standard procedure: updates cue
      dtimer(dtime,itime,-1)
      wbpush2.wcuperp2(cue,nx,ny)
      dtimer(dtime,itime,1)
      tfield[0] += float(dtime)

# calculate electromagnetic fields in fourier space with standard
# procedure: updates exyz, bxyz, wf, wm
      dtimer(dtime,itime,-1)
      if (ntime==0):
         wbpush2.wibpois2(cue,bxyz,ffc,ci,wm,nx,ny)
         wf[0] = 0.0
         dth = 0.5*dt
      else:
         wbpush2.wmaxwel2(exyz,bxyz,cue,ffc,ci,dt,wf,wm,nx,ny)
      dtimer(dtime,itime,1)
      tfield[0] += float(dtime)

# calculate force/charge in fourier space with standard procedure:
# updates fxyze, we
      dtimer(dtime,itime,-1)
      isign = -1
      wbpush2.wpois2(qe,fxyze,ffc,we,nx,ny)
      dtimer(dtime,itime,1)
      tfield[0] += float(dtime)

# add longitudinal and transverse electric fields with standard
# procedure: updates fxyze
      dtimer(dtime,itime,-1)
      isign = 1
      wbpush2.wemfield2(fxyze,exyz,ffc,isign,nx,ny)
# copy magnetic field with standard procedure: updates bxyze
      isign = -1
      wbpush2.wemfield2(bxyze,bxyz,ffc,isign,nx,ny)
      dtimer(dtime,itime,1)
      tfield[0] +=float(dtime)

# transform electric force to real space with standard procedure:
# updates fxyze
      dtimer(dtime,itime,-1)
      isign = 1
      wbpush2.w3fft2r(fxyze,isign,mixup,sct,indx,indy)
      dtimer(dtime,itime,1)
      tfft[0] += float(dtime)

# transform magnetic force to real space with standard procedure:
# updates bxyze
      dtimer(dtime,itime,-1)
      isign = 1
      wbpush2.w3fft2r(bxyze,isign,mixup,sct,indx,indy)
      dtimer(dtime,itime,1)
      tfft[0] += float(dtime)

# copy guard cells with standard procedure: updates fxyze, bxyze
      dtimer(dtime,itime,-1)
      wbpush2.wcguard2(fxyze,nx,ny)
      wbpush2.wcguard2(bxyze,nx,ny)
      dtimer(dtime,itime,1)
      tguard[0] += float(dtime)

# push particles with standard procedure: updates part, wke
      wke[:] = 0.0
      dtimer(dtime,itime,-1)
# relativistic pushers
      if (relativity==1):
# analytic Boris pusher
         if (kpush==1):
            wbpush2.wgrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,ny)
# analytic Boris pusher
         elif (kpush==2):
            wabpush2.wagrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,
                               ny)
# exact analytic pusher
         elif (kpush==3):
            wabpush2.weagrbpush2(part,fxyze,bxyze,qbme,dt,dth,ci,wke,nx,
                                ny)
# non-relativistic pushers
      else:
# classic Boris pusher
         if (kpush==1):
            wbpush2.wgbpush2(part,fxyze,bxyze,qbme,dt,dth,wke,nx,ny)
# analytic Boris pusher
         elif ((kpush==2) or (kpush==3)):
            wabpush2.wagbpush2(part,fxyze,bxyze,qbme,dt,dth,wke,nx,ny)
      dtimer(dtime,itime,1)
      tpush[0] += float(dtime)

# sort particles by cell for standard procedure
      if (sortime > 0):
         if (ntime%sortime==0):
            dtimer(dtime,itime,-1)
            wbpush2.wdsortp2yl(part,part2,npicy)
# exchange pointers
            tpart = part
            part = part2
            part2 = tpart
            dtimer(dtime,itime,1)
            tsort[0] += float(dtime)

      if (ntime==0):
         wt[0] = we[0] + wf[0] + wm[0]
         print("Initial Total Field, Kinetic and Total Energies:")
         print("%14.7e %14.7e %14.7e" % (wt, wke, wke + wt))
         print("Initial Electrostatic, Transverse Electric and Magnetic "\
               "Field Energies:")
         print("%14.7e %14.7e %14.7e" % (we, wf, wm))

      ntime += 1

# * * * end main iteration loop * * *

   print("ntime, relativity, kpush = ",ntime,relativity,kpush)
   wt[0] = we[0] + wf[0] + wm[0]
   print("Final Total Field, Kinetic and Total Energies:")
   print("%14.7e %14.7e %14.7e" % (wt, wke, wke + wt))
   print("Final Electrostatic, Transverse Electric and Magnetic "\
         "Field Energies:")
   print("%14.7e %14.7e %14.7e" % (we, wf, wm))

   print("")
   print("deposit time = ",tdpost[0])
   print("current deposit time = ",tdjpost[0])
   tdpost[0] += tdjpost[0]
   print("total deposit time = ",tdpost[0])
   print("guard time = ",tguard[0])
   print("solver time = ",tfield[0])
   print("fft time = ",tfft[0])
   print("push time = ",tpush[0])
   print("sort time = ",tsort[0])
   tfield[0] += tguard[0] + tfft[0]
   print("total solver time = ",tfield[0])
   time = tdpost[0] + tpush[0] + tsort[0]
   print("total particle time = ",time)
   wt[0] = time + tfield[0]
   print("total time = ",wt[0])
   print("")

   wt[0] = 1.0e+09/(float(nloop)*float(np))
   print("Push Time (nsec) = ",tpush[0]*wt[0])
   print("Deposit Time (nsec) = ",tdpost[0]*wt[0])
   print("Sort Time (nsec) = ",tsort[0]*wt[0])
   print("Total Particle Time (nsec) = ",time*wt[0])

if (__name__=="__main__"):
   main()
