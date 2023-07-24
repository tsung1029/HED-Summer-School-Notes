#-----------------------------------------------------------------------
from __future__ import print_function
import math
import numpy
import matplotlib.pyplot as plt

from rplasmalib import *

int_type = numpy.int32
double_type = numpy.float64
float_type = numpy.float32
   
#-----------------------------------------------------------------------
def main():
# Inputs:
# indx = exponent which determines grid points in x direction:
# nx = 2**indx.
   indx =   9
# npx = number of electrons distributed in x direction.
   npx = 18432
# nloop = number of time steps in the simulation
   nloop = 100
# dt = time interval between successive calculations.
   dt = 0.1
# vtx = thermal velocity of electrons in x direction
# vx0 = drift velocity of electrons in x direction.
   vtx = 1.0; vx0 = 0.0
# ci = reciprocal of velocity of light.
   ci = 0.1
# idimp = number of particle coordinates = 2
   idimp = 2

# initialize scalars for standard code
# nx = number of grid points in x direction
# np = total number of particles in simulation
   nx = int(math.pow(2,indx)); np = npx; 

# allocate data for standard code
# rpart = helper list for relativistic particles
#  rpart = [0.0,0.0,idimp,np]
# set hidden Fortran helper object for relativistic particles
   wrpush1.set_part1d(0.0,0.0,idimp,np)
# wke = particle kinetic energy
   wke = numpy.zeros((1),float_type)
# part = particle array
   part = numpy.empty((idimp,np),float_type,'F')

# initialize uniform plasma and maxwellian velocity: updates part
   init1.distr1(part,vtx,vx0,npx,nx,0)

   print('initial:part(:,0)=',part[:,0])

# * * * start main iteration loop * * *

   for n in range(0,nloop):

# push free-streaming relativstic particles: updates part, wke
      wke[:] = 0.0
# packing arguments from helper list
#     wrpush1.w1rpush1zf(rpart[0],rpart[1],rpart[2],rpart[3],part,dt,ci,
#                        wke,nx)
# using hidden Fortran helper object
      wrpush1.w2rpush1zf(part,dt,ci,wke,nx)

# * * * end main iteration loop * * *

   print('final:part(:,1)=',part[:,0])

   plt.plot(part[0,:],part[1,:],'o',markersize=0.1)
   plt.title("X VS VX: Final Time")
   plt.show()

if (__name__=="__main__"):
   main()
