#-----------------------------------------------------------------------
"""
test code to generate physics data

written by Viktor K. Decyk, UCLA
update: July 19, 2023
"""
from __future__ import print_function
import numpy
import comms

global sfield
sfield = None

def main():
   """ generate physics data """
   global sfield

# nx = number of grid points in x direction
   nx = 128; nxe = nx + 2
# tend = time at end of simulation, in units of plasma frequency.
# dt = time interval between successive calculations.
   tend = 10.0; dt = 0.2
# nloop = number of time steps in simulation
# nstart = initial time loop index
# ntime = current time step
   nstart = 0; ntime = 0
   nloop = int(tend/dt + .0001) + 1

   mode = 3; omega = 1.0
   dnx = 2.0*numpy.pi*float(mode)/float(nx)
   dnh = 0.5*dnx
   sfield = numpy.empty((nxe),numpy.float32,'F')

   print("starting physics")
# start GUI
# send plot_loc, sim_data dictionaries and custom event to initialize
# plots
   comms.update_gui({"SIN(KX-WT) VS X": 0, "SIN(KX+WT) VS X": 1,
                     "COS(.5*KX-WT) VS X": 2},
                    {"DT": 0.1,"TEND": tend,"NX": nx})
# wait for plot intialization
   gui_err = comms.check_run_status()
   if (gui_err=='QUIT'):
      exit()

   for ntime in range(nstart,nloop):
      time = dt*float(ntime)
# send current time to GUI
      comms.update_time(time)

# calculate sin(kx-wt)
      sfield[:nx+1] = numpy.sin(dnx*numpy.array(range(nx+1))-omega*time)
      comms.update_plot('SIN(KX-WT) VS X')
      gui_err = comms.check_plot_status()
      if gui_err != comms.plot_name:
         if (gui_err=='QUIT'):
            break

# calculate sin(kx+wt)
      sfield[:nx+1] = numpy.sin(dnx*numpy.array(range(nx+1))+omega*time)
      comms.update_plot('SIN(KX+WT) VS X')
      gui_err = comms.check_plot_status()
      if gui_err != comms.plot_name:
         if (gui_err=='QUIT'):
            break

# calculate cos(.5*kx-wt)
      sfield[:nx+1] = numpy.cos(dnh*numpy.array(range(nx+1))-omega*time)
      comms.update_plot('COS(.5*KX-WT) VS X')
      gui_err = comms.check_plot_status()
      if gui_err != comms.plot_name:
         if (gui_err=='QUIT'):
            break

# wait for GUI at each loop iteration
      gui_err = comms.check_run_status()
      if (gui_err=='QUIT'):
         break

# Check if GUI needs to quit
   if (gui_err=='QUIT'):
      exit()

   print("done physics")

if (__name__=="__main__"):
   main()
