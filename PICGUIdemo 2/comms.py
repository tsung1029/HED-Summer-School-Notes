#-----------------------------------------------------------------------
"""
Library for communicating between the MVC components of the PIC Gui
Shared global data for the MVC components is located here
The Control Panel has the tkinter GUI elements and needs to run the main
thread,  The view graphics display is running in the same thread.
The model physics is running in a different thread.

root = root window for the Control Panel is defined here so that the
       model physics can safely send events to the Control Panel
runq = run queue for commmunicating time loop information between Control
       Panel and physics which run in different threads
plotq = = plot queue for communicating plot information between Control
          Panel and physics which run in different threads
rtime = rounded current time from physics model
stime = rounded current time string for display by the Control Panel
ftime = future time when graphics is reactivated
jtime = jump time, ftime = rtime + jtime
rstatus = run status (0=idle, 1=running)
plot_name = current plot being plotted to communicate to Control Panel
plot_loc = plot name and number, used to determine location to display
sim_data = simulation parameters needed by Control Panel and display
sdt = simulation time step
stend = simulation end time
snx = simulation spatial size

written by Viktor K. Decyk, UCLA, with contributions from Aileen Wang
update: march 29, 2023
"""
global root, runq, plotq, rtime, stime, ftime, jtime, rstatus
global plot_name, plot_loc, sim_data, sdt, stend, snx

root = None
runq = None
plotq = None
rtime = 0.0
stime = ''
ftime = -1.0
jtime = 0.0
rstatus = 0
plot_name = ''
plot_loc = None
sim_data = None
sdt = 0.0
stend = 0.0
snx = 0

def update_gui(plotloc,simdata):
   global root, plot_loc, sim_data, sdt, stend, snx
   """
   Send plot_loc and sim_data dictionaries and custom event from physics
   to Control Panel
   """
   plot_loc = plotloc
   sim_data = simdata
   sdt = sim_data['DT']
   stend = sim_data['TEND']
   snx = sim_data['NX']
   root.event_generate('<<StartGUI>>',when='tail')

def update_time(time):
   """
   Send current time and custom event from physics to Control Panel
   """
   global root, rtime, stime, ftime
#round off time to 4 digits
   rtime = round(time,4)
   stime = str(round(time,3))
   root.event_generate('<<TimeChange>>',when='tail')

def update_plot(pname,plotdata=None,plotdata2=None):
   """
   Send plot_name and custom event from physics to Control Panel to
   initiate plot if rtime >= ftime or rstatus = 0.
   """
   global root, rtime, ftime, rstatus, plot_name
#request plot
   plot_name = pname
   if rtime >= ftime or rstatus==0:
      root.event_generate('<<PlotChange>>',when='tail')
#suppress plot if running and rtime < ftime
   else:
      set_plot_status(pname)

def check_run_status():
   """ wait for run to be allowed to proceed """
   global runq
   err = runq.get()
   if err != 'GO':
      if err != 'QUIT':
         print('Time Step Error:',err)
         err = 'QUIT'
   return err
   
def set_run_status():
   """ send message to allow run to proceed """
   global runq
   runq.put('GO')

def check_plot_status():
   """ wait for plot to be allowed to proceed """
   global plotq
   err = plotq.get()
   if err != plot_name:
      if err != 'QUIT':
         print(plot_name,' Plot Error:',err)
         err = 'QUIT'
   return err


def set_plot_status(pname):
   """ send message to allow plot to proceed """
   global plotq
   plotq.put(pname)

def set_pause():
   """ send message to pause if running """
   root.event_generate('<<EndLoop>>',when='tail')

