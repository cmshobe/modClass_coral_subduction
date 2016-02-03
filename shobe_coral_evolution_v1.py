# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:12:12 2016

@author: Charlie

Modeling class, week 3.

Model for coral growth and drowning on a subducting
or uplifting plate (after Galewsky 1998).

Adams-Bashforth 3 solver
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

###############Initialize
y_min = -1000 #for plot
y_max = 0 #for plot
t_plot = 10000#years
t_min = 0 #years
t_max = 10000000 #years
dt = 10 #yrs
times = np.arange(t_min, t_max + dt, dt) #time vector
#parameters:
#mean sea level (DATUM)
#msl = 0 #m relative to msl
#coral growth
max_up_growth_rate = .012 #m/yr
extinct_coeff = .16 #/m
surf_light_intens = 2250 #uE m-2 s-1
saturat_light_intens = 450 #uE m-2 s-1

#plate and flexure
initial_pos_line_load = 250000 #m, end of x domain
convergence_rate = .1 #m/yr
max_deflection = 10000 #m
g = 9.81 #m/s^2
flex_rigid = 1e22 #Nm
dens_mantle = 3300 #kg/m3
dens_ocean = 1030 #kg/m3
flex_param = np.power(4 * flex_rigid / ((dens_mantle - dens_ocean) * g), 1/4)

#time domain
#spatial domain, sloping plate with line load
#ALL DEPTHS RELATIVE TO MEAN SEA LEVEL, WHICH IS AN UNCHANGING DATUM
subducting_plate_x = np.arange(0, 200000, 100) #meters
bedrock_height = np.arange(0, -2000, -1) #meters. down-to-right tilt
coral_height = np.zeros((len(bedrock_height)))

sea_level_condition = 60 * np.sin(times / (5000 * np.pi))
up_growth_rate = np.zeros((len(bedrock_height)))
###############Run
pos_line_load = initial_pos_line_load
##within loop:
schematic_fig = plt.figure(figsize=(4,8)) #instantiate figure
schematic = plt.subplot()
#plt.gca().invert_xaxis()
schematic.set_ylim(600, 800)
schematic.set_xlim(0, 40000)
plt.ion()
plt.show()

x = pos_line_load - subducting_plate_x #distance from load to point
deflection = max_deflection * np.exp(-(x / flex_param)) * np.cos(x / flex_param)
bedrock_height = -deflection
msl = max(bedrock_height-100)
run_count = 0
for t in range(len(times)-1):
    run_count +=1 #keep track of iterations
    current_time = times[t]
    #print run_count
    rel_sl = msl + sea_level_condition[run_count - 1]
    #calculate the deflection of the plate (analytical flexure equation)
    x = abs(pos_line_load - subducting_plate_x) #distance from load to point
    deflection = max_deflection * np.exp(-(x / flex_param)) * np.cos(x / flex_param)
    bedrock_height = -deflection
    #grow coral (will need editing if using fancy solver)
    depth = rel_sl - (bedrock_height + coral_height)
    #np.where(depth > 0, up_growth_rate = 0, up_growth_rate = 99)
    up_growth_rate[depth > 50] = 0
    up_growth_rate[(depth >= 0) & (depth <= 50)] = max_up_growth_rate * np.tanh(surf_light_intens * np.exp(-extinct_coeff * depth[(depth >= 0) & (depth <= 50)]) / saturat_light_intens)
    up_growth_rate[depth < 0] = 0
    #print max(up_growth_rate)
    coral_height += (up_growth_rate * dt) #THIS IS THE SOLVER STEP!!!
    
    #advance position of the line load towards the foreland by [conv rate * dt]
    pos_line_load -= (convergence_rate * dt)
    plate_for_plot = subducting_plate_x[subducting_plate_x <= pos_line_load]
    #change RELATIVE SEA LEVEL (tectonic subsidence + eustatic sea level change)
    #plot schematic of model space
    if current_time % t_plot == 0:
        print current_time
        schematic.clear()
        schematic.plot(plate_for_plot, bedrock_height[0:len(plate_for_plot)], color='k', label='Bedrock')
        schematic.plot(plate_for_plot, bedrock_height[0:len(plate_for_plot)]+coral_height[0:len(plate_for_plot)], color='c', label='Coral')
        schematic.plot((min(plate_for_plot), max(plate_for_plot)), (rel_sl, rel_sl), color='b', label='Water')
        schematic.plot((pos_line_load, pos_line_load), (y_min, y_max), label='Line load')
        schematic.set_ylim(-4000, 800)
        schematic.set_xlim(0, 200000)
        #schematic.clear()    
        #plt.draw()
        plt.pause(0.2)
    else:
        pass
################Finalize