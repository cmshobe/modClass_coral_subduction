# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:12:12 2016

@author: Charlie

Modeling class, week 3.

Model for coral growth and drowning on a subducting
or uplifting plate (after Galewsky 1998).
"""
from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 20})
###############Initialize
y_min = -1000 #for plot
y_max = 0 #for plot
t_plot = 10000#years
t_min = 0 #years
t_max = 600000 #years
dt = 10 #yrs
times = np.arange(t_min, t_max + dt, dt) #time vector
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
x_min = 0 #m
x_max = 200000 #m
dx = 100 #m
subducting_plate_x = np.arange(x_min, x_max, dx) #meters
bedrock_height = np.zeros((len(subducting_plate_x))) #meters. down-to-right tilt
coral_height = np.zeros((len(bedrock_height)))

sea_level_condition = 60 * np.sin(times / (5000 * np.pi))
up_growth_rate = np.zeros((len(bedrock_height)))
###############Run
pos_line_load = initial_pos_line_load
##within loop:
schematic_fig = plt.figure(figsize=(10,4)) #instantiate figure
schematic = plt.subplot()
plt.gcf().subplots_adjust(bottom=0.20)
plt.xlabel('Distance from Subduction Zone [km]')
plt.ylabel('Elevation above datum [m]')
schematic.set_ylim(600, 800)
schematic.set_xlim(0, 40000)
plt.ion()
plt.show()

x = pos_line_load - subducting_plate_x #distance from load to point
x_for_plotting = np.arange(x_min/1000, x_max/1000, dx/1000)
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
    depth = rel_sl - (bedrock_height + coral_height)
    up_growth_rate[depth > 50] = 0
    up_growth_rate[(depth >= 0) & (depth <= 50)] = max_up_growth_rate * np.tanh(surf_light_intens * np.exp(-extinct_coeff * depth[(depth >= 0) & (depth <= 50)]) / saturat_light_intens)
    up_growth_rate[depth < 0] = 0
    coral_height += (up_growth_rate * dt) #THIS IS THE SOLVER STEP!!!
    pos_line_load -= (convergence_rate * dt)
    plate_for_plot = subducting_plate_x[subducting_plate_x <= pos_line_load]
    x_for_plotting = x_for_plotting[subducting_plate_x/1000 <= pos_line_load/1000]
    if current_time % t_plot == 0:
        #print current_time
        schematic.clear()
        schematic.plot(x_for_plotting, bedrock_height[0:len(plate_for_plot)], color='.5', label='Bedrock')
        schematic.plot(x_for_plotting, bedrock_height[0:len(plate_for_plot)]+coral_height[0:len(plate_for_plot)], color='c', label='Coral')
        schematic.plot((min(x_for_plotting), max(subducting_plate_x)), (rel_sl, rel_sl), color='b', label='Water')
        #schematic.plot((pos_line_load/1000, pos_line_load/1000), (y_min, y_max), label='Line load')
        schematic.fill_between(x_for_plotting, -4000, bedrock_height[0:len(plate_for_plot)], facecolor='.5', interpolate=True)
        schematic.fill_between(subducting_plate_x/1000, bedrock_height+coral_height, rel_sl, where=rel_sl >bedrock_height+coral_height, facecolor='b', interpolate=True)
        schematic.fill_between(x_for_plotting, bedrock_height[0:len(plate_for_plot)], bedrock_height[0:len(plate_for_plot)]+coral_height[0:len(plate_for_plot)], where=bedrock_height[0:len(plate_for_plot)]+coral_height[0:len(plate_for_plot)] >bedrock_height[0:len(plate_for_plot)], facecolor='c', interpolate=True)
        schematic.set_ylim(-4000, 800)
        schematic.set_xlim(0, 200)
        plt.xlabel('Distance from Subduction Zone [km]')
        plt.ylabel('Elevation above datum [m]')
        plt.text(5, -3500, 'Time [yrs]: %.1f' % current_time)
        plt.pause(0.2)
    else:
        pass
################Finalize