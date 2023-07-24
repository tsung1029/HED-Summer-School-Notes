#-----------------------------------------------------------------------
"""
Display program for PIC Gui
This is the graphics (View) component of the PIC GUI.
It displays multiple plots requested by the physics (model) component.
There can be up to 4 plots per window and mutiple windows.  The
dictionary comms.plot_loc contains the plot names and identifiers.
The dictionary comms.sim_data contains additional informationb that
may be needed by the plots.

The physics component communicates with the Control Panel by generating
a custom event <<PlotChange>> to request a plot.  The Control Panel then
calls a plot function in the graphics component.  The graphics component 
communicates with the physics component by writing to the queue plotq to
synchronize the plots.

Main graphic procedure is:
dscaler1 displays 1d scalar field in real space.

written by Viktor K. Decyk, UCLA, with contributions from Aileen Wang
update: may 20, 2023
"""
from __future__ import print_function
import sys
import tkinter as tk
from tkinter import ttk
import math
import cmath
import numpy
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
if (sys.version_info.major==2):
   from matplotlib.backends.backend_tkagg import \
      NavigationToolbar2TkAgg as NavigationToolbar2Tk
else:
   from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
import comms

global win, num_wins, frame, fig, canvas, toolbar
num_wins = 0
win = []; frame = []; fig = []; canvas = []; toolbar = []

def main():
   global win, num_wins, frame, fig, canvas, toolbar

#Find location of control panel
   x = comms.root.winfo_x() + comms.root.winfo_width()
   y = comms.root.winfo_y()
#Find number of plots required from plot_loc dictionary
   num_plots = len(comms.plot_loc)
   num_wins = int((num_plots - 1)/4) + 1
#Default vertical and horizontal Figure size, in inches
   vfsz = 6; hfsz = 4
#Create windows
   for j in range(0,num_wins):
      nfs = min(4,max(0,num_plots-4*j))
      win.append(tk.Tk())
      win[j].title('Plot Window %d'%(j))
#Move window to the right of control panel
      win[j].geometry('+%d+%d'%(x+10+30*j,y+30*j))
      win[j].columnconfigure(0,weight=1)
      win[j].rowconfigure(0,weight=1)
#Find number of plots in window
      nfs = min(4,max(0,num_plots-4*j))
#Create frames in window
      for i in range(0,nfs):
         ii = i + 4*j
         frame.append(ttk.Frame(win[j]))
         nc = int(i/2)
         nr = i - 2*nc
         win[j].rowconfigure(nr, weight=1)
         win[j].columnconfigure(nc, weight=1)
         fig.append(Figure(figsize=(vfsz,hfsz),dpi=100))
         canvas.append(FigureCanvasTkAgg(fig[ii],master=frame[ii]))
         toolbar.append(NavigationToolbar2Tk(canvas[ii],frame[ii]))
         frame[ii].grid(row=nr,column=nc,columnspan=1,sticky='news')
         canvas[ii].get_tk_widget().pack(expand=True,fill=tk.BOTH)
         toolbar[ii].pack(expand=True,fill=tk.BOTH)
      win[j].update_idletasks()

#Send a message to start physics time loop
   comms.set_run_status()

#-----------------------------------------------------------------------
def dscaler1(f,label,isc,ist,nx,chr):
   """ displays 1d scalar field in real space """
   global frame, fig, canvas, toolbar
# Find window and frame number
   i = comms.plot_loc[label]
   j = int(i/4)
# Clear figure
   fig[i].clf()

   fglabel = label + ', TIME = ' + comms.stime

   axes = fig[i].add_subplot(111)
   axes.set_title(fglabel)
   axes.set_xlabel(chr)
   axes.set_xlim(xmin=0.0,xmax=nx)
   axes.plot(f[:nx])
   canvas[i].draw()
   canvas[i].get_tk_widget().update_idletasks()
   toolbar[i].update_idletasks()

# set status to done
   comms.set_plot_status(label)

#-----------------------------------------------------------------------
def destroy_windows():
   """ callback for Quit """
   global win, num_wins
#Destroy windows
   for j in range(0,num_wins):
      win[j].destroy()

if (__name__=="__main__"):
   main()
