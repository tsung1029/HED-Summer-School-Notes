#-----------------------------------------------------------------------
"""
Control Panel for the MVC components of the PIC Gui
The Control Panel creates a menu and a Panel which contains:
A label to display the current physics time step
A radiobutton to select the run mode, Continuous or Time Chunks.
In Continuous mode, there is a Fast Forward Time option
If Fast Forward option is selected, then display is suppressed until the
Fast Forward Time is reached.
In Time Chunk mode, there is a Jump Time option
If Jump Time option is selected, then display is suppressed until the
Fast Forward time is reached, where Fast Forward time is the current
time plus the jump time.
There is also a Step mode where display is not suppressed and the code
pauses after every time step.

The physics component communicates with the Control Panel by
generating the custom event <<StartGUI>> and communicating the plot_loc
and sim_data dictionaries.  The custom event initializes the graphics
component.  The data being communicatedd is stored in global variables
located in the module comms.py.

The Control Panel launches the graphics (View) component in the main 
thread, and the physics (Model) component in a separate thread.
Shared data structures are global constants in the comms module, which
can be written and read by all the other components.

The physics component communicates with the Control Panel by
generating the custom event <<TimeChange>> to communicate the current
physics simulation time.  The Control Panel communicates with the
physics component by writing to the queue runq to manage the physics
time loop so that the physics and graphics stay synchronized.

The physics component also communicates with the Control Panel by
generating the custom event <<PlotChange>> to request a plot.  The
Control Panel then calls a plot function in the graphics component.
The graphics component communicates with the physics component by
writing to the queue plotq to synchronize the individual plots.

At the end of the simulation loop, the physics component communicates
with the Control Panel by generating the custom event <<EndLoop>> to
pause the display.

The physics component never communicates directly with the graphics
component.

The startup proceeds as follows:
First the cpanel.py script establishes the main thread, sets up the
control panel, binds a listener for custom events, starts up the physics
module physics.py in a separate thread, and enters a mainloop.
Then the physics.py module initializes the physics data, calls the
comms.update_gui function to copy plot dictionaries to comms.py, and
generates a custom event to start the GUI.  Finally, the cpanel.py
script responds to this custom event and initializes the plot script
graphs.py.

This version of cpanel uses a dictionary instead of an if-then-else
block in callback function on_plotstart to execute plotting functions

written by Viktor K. Decyk, UCLA, with contributions from Aileen Wang
update: may 16, 2023
"""
from __future__ import print_function
import tkinter as tk
from tkinter import ttk
from threading import Thread
from queue import Queue
import cmath

import physics
import comms
import graphs

def on_quit():
#quit physics
   comms.runq.put('QUIT')
   comms.plotq.put('QUIT')
#quit graphics
   try:
      graphs.destroy_windows()
   except:
      pass
#quit control panel
   comms.root.destroy()

def main():
   global plot_func
   plot_func = None

   comms.root = tk.Tk()
   comms.root.title('Control Panel')

#Create menu
   main_menu = tk.Menu(comms.root)
   comms.root.config(menu=main_menu)
   text_menu = tk.Menu(main_menu,tearoff=False)
   text_menu.add('command',label='Quit PIC Gui',command=on_quit)
   main_menu.add_cascade(label='PIC Gui',menu=text_menu)
   
#Set default control variable objects
#Time label
   tlabel = tk.StringVar()
#Simulation end time label
   tendlabel = tk.StringVar()
#Run mode: rmode = (Continuous=1,Chunk=2)
   rmode = tk.IntVar()
#Run mode command label
   lrunc = tk.StringVar()
#Run mode option label
   lvalue = tk.StringVar()
#Run mode option: skip graphics until some future time
   tvalue = tk.DoubleVar()
#Display messages when inputting run mode option
   errmsg = tk.StringVar()

   tlabel.set('0.0')
   tendlabel.set('0.0')
   rmode.set(1)
   lrunc.set('Run Continuously')
   lvalue.set('Fast Forward to Time:')
   tvalue.set('')
   errmsg.set('time values must be positive')

#Create queues for communication between threads
   comms.runq = Queue()
   comms.plotq = Queue()

#Define call back functions
   def on_contig():
      """ callback for continuous radiobutton """
      lrunc.set('Run Continuously')
      lvalue.set('Fast Forward to Time:')
      tvalue.set('')

   def on_chunk():
      """ callback for chunks radiobutton """
      lrunc.set('Run Chunk')
      lvalue.set('Jump Time:')
      tvalue.set('1.0')
      comms.jtime = 1.0

   def on_button():
      """ callback for runc button """
      rc = rmode.get()
#Idle, set run status to running
      if comms.rstatus==0:
#Calculate fast forward time
         if (rc==2):
            tf = comms.rtime + comms.jtime
#Round to nearest integer time step
            tf = round(comms.sdt*float(int(tf/comms.sdt + 0.5)),4)
            comms.ftime = tf
#Send a message to continue time loop
         comms.set_run_status()
         lrunc.set('Pause')
         comms.rstatus = 1
#Disable step button
         step.config(state='disabled')
#Disable Display Run Mode option entry
         descript.config(foreground='grey')
         value.config(state='disabled')
#Running, set run status to idle
      elif comms.rstatus==1:
         comms.rstatus = 0
         if rc==1:
           lrunc.set('Run Continuously')
           tvalue.set('')
         elif rc==2:
           lrunc.set('Run Chunk')
#Enable Display Run Mode option entry
         value.config(state='normal')
         descript.config(foreground='black')
#Enable step button
         step.config(state='normal')

   def on_future(*args):
      """ callback for value entry widget """
#For run chunk, add current time
      rc = rmode.get()
#For run continuous, get ftime
      if (rc==1):
         tf = tvalue.get()
#Round to nearest integer time step
         tf = round(comms.sdt*float(int(tf/comms.sdt + 0.5)),4)
         comms.ftime = tf
         errmsg.set('time values must be positive')
#For run chunk, get jump time
      elif (rc==2):
         comms.jtime = tvalue.get()

   def set_plotfunc():
      """ Create plotting function dictionary plot_func"""
      global plot_func
      plot_func = {}
      for key, value in comms.plot_loc.items():
         if key=='SIN(KX-WT) VS X':
            comms.plot_name = key
            plot_func[value] = compile("graphs.dscaler1(physics.sfield,\
                                       comms.plot_name,999,0,\
                                       comms.snx+1,'X')","<string>",
                                       "exec")
         elif key=='SIN(KX+WT) VS X':
            comms.plot_name = key
            plot_func[value] = compile("graphs.dscaler1(physics.sfield,\
                                       comms.plot_name,999,0,\
                                       comms.snx+1,'X')",
                                       "<string>","exec")
         elif key=='COS(.5*KX-WT) VS X':
            comms.plot_name = key
            plot_func[value] = compile("graphs.dscaler1(physics.sfield,\
                                        comms.plot_name,999,0,\
                                        comms.snx+1,'X')",
                                        "<string>","exec")

#Callbacks for custom events
   def on_gui_start(*args):
      """ callback for StartGUI custom event """
#Set end time label from physics, this will update the display
      tendlabel.set(comms.stend)
# Execute graphics module
      graphs.main()
#Create dictionary associating a plot location with a plotting function
      set_plotfunc()
      return

   def on_timechange(*args):
      """ callback for TimeChange custom event """
#Set time label from physics, this will update the display
      tlabel.set(comms.stime)
#If future time has arrived, pause graphics
      diff = round(abs(comms.ftime-comms.rtime),4)
      if diff < 0.5*comms.sdt:
         on_button()
#If running, send a message to continue time loop,
      if comms.rstatus==1:
         comms.set_run_status()

   def on_plotstart(*args):
      """
      This callback for PlotChange custom event uses a dictionary
      plot_func which associates an executable plotting function with a
      plot location
      """
      global plot_func
#First finds the plot location for plotname
      i = comms.plot_loc[comms.plot_name]
#Then executes the associated plotting function
      exec(plot_func[i])

   def on_endloop(*args):
      """ callback for EndLoop custom event """
#If running, pause
      if comms.rstatus==1:
         on_button()

#Bind functions to custom events
   comms.root.bind('<<StartGUI>>',on_gui_start)
   comms.root.bind('<<TimeChange>>',on_timechange)
   comms.root.bind('<<PlotChange>>',on_plotstart)
   comms.root.bind('<<EndLoop>>',on_endloop)

#Create Panels

   frame = ttk.Frame(comms.root)

#Display current simulation time
   ltime = ttk.Label(frame,text='Simulation Time:')
   ltval = ttk.Label(frame,width=8,textvariable=tlabel)

#Display final simulation time
   ltend = ttk.Label(frame,text='Simulation  End:')
   ltendval = ttk.Label(frame,width=8,textvariable=tendlabel)

#Display Run Mode choices: Continuous or Time Chunk
   lrmode = ttk.Label(frame,text='Run Mode')
   continuous = ttk.Radiobutton(frame,text='Continuous',variable=rmode,
                                value=1,command=on_contig)
   chunks = ttk.Radiobutton(frame,text='Time Chunks',variable=rmode,
                            value=2,command=on_chunk)

#Display Run Mode command: Run Continuously or Run Chunk or Pause
   runc = ttk.Button(frame,textvariable=lrunc,command=on_button)
   
#Display Run Mode option: Fast Forward or Jump Time
   descript = ttk.Label(frame,width=18,textvariable=lvalue)
   value = ttk.Entry(frame,width=8,textvariable=tvalue)
   value.bind('<FocusOut>',on_future)

#Error handling for Entry widget
   value_error = ttk.Label(frame,textvariable=errmsg)
   def only_positive(proposed):
      """ Verify keys entered represent a positive number """
      if (proposed==''):
         errmsg.set('')
         return True
      try:
         num = float(proposed)
         if (num > 0):
            errmsg.set('')
            return True
         else:
            return False
      except ValueError:
         errmsg.set('Value must be a positive number!')
         return False
   validate_ref = value.register(only_positive)
   value.configure(validate='key',validatecommand=(validate_ref,'%P'))

#Run one time step
   step = ttk.Button(frame,text='Step',command=comms.set_run_status)

#Quit
   quit = ttk.Button(frame,text='QUIT',command=on_quit)

#Layout Panel
   frame.grid(row=0,column=0)
   ltime.grid(row=0,column=0)
   ltval.grid(row=0,column=1)
   ltend.grid(row=1,column=0,pady=8)
   ltendval.grid(row=1,column=1,pady=8)
   lrmode.grid(row=2,columnspan=2)
   continuous.grid(row=3,column=0)
   chunks.grid(row=3,column=1,padx=16)
   runc.grid(row=4,columnspan=2,pady=16)
   descript.grid(row=5,column=0)
   value.grid(row=5,column=1)
   value_error.grid(row=6,columnspan=2)
   step.grid(row=7,columnspan=2,pady=16)
   quit.grid(row=8,column=1,pady=8)

   frame.update_idletasks()

# Launch physics module as a separate thread
   Thread(target=physics.main).start()

   comms.root.update_idletasks()
   comms.root.mainloop()

if (__name__=="__main__"):
   main()
