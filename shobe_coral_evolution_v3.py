# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 14:12:12 2016

@author: Charlie

Modeling class, week 3.

Model for coral growth and drowning on a subducting
or uplifting plate (after Galewsky 1998).

Driven by crudely "scaled" dO18 record.
Not certain whether I'm running the record backwards or forwards.

USER GUIDE:
-all user defined variables are set in __name__ block at bottom
-this will export a figure called 'shobe_coral.png' to your working directory.
"""
from __future__ import division #true division
import numpy as np #numeric python
import matplotlib #plotting
import matplotlib.pyplot as plt #plotting
matplotlib.rcParams.update({'font.size': 20})  #plot font size       
        
class Coral: #container to hold all the model functions
    def flex(self, max_deflec, x, flex_param): #does flexure calcs
        deflection = max_deflec * np.exp(-(x / flex_param)) * np.cos(x / flex_param)
        return deflection
        
    def grow_coral(self, depth, up_growth_rate, max_up_growth_rate, surf_light_intens, extinct_coeff, saturat_light_intens):
        self.up_growth_rate[depth > 50] = 0
        self.up_growth_rate[(depth >= 0) & (depth <= 50)] = max_up_growth_rate * np.tanh(surf_light_intens * np.exp(-extinct_coeff * depth[(depth >= 0) & (depth <= 50)]) / saturat_light_intens)
        self.up_growth_rate[depth < 0] = 0
        return self.up_growth_rate
        
    def initialize(self, x_min, x_max, dx, t_min, t_max, dt, init_pos_line_load, g, flex_rigid, dens_mantle, dens_ocean, file_name):
        self.times = np.arange(t_min, t_max + dt, dt) #time vector
        self.subducting_plate_x = np.arange(x_min, x_max, dx) #space vector
        self.flex_param = np.power(4 * flex_rigid / ((dens_mantle - dens_ocean) * g), 1/4)
        self.bedrock_height = np.zeros((len(self.subducting_plate_x)))
        self.coral_height = np.zeros((len(self.subducting_plate_x)))
        self.up_growth_rate = np.zeros((len(self.subducting_plate_x)))
        self.pos_line_load = init_pos_line_load   
        if isotope_or_sinusoid == 0: #isotope record drives sea level
            isotope_record = np.genfromtxt(file_name, dtype=float)
            isotope_record_times = isotope_record[0:t_max/1000+1, 0] * 1000
            isotope_record_values = isotope_record[0:t_max/1000+1, 1]
            isotope_scaled_to_sea_level = ((5.09 - isotope_record_values) * 60) - 60           
            interpolated_sea_levels = np.interp(self.times, isotope_record_times, isotope_scaled_to_sea_level)
            sea_level_condition = interpolated_sea_levels
        else: #sinusoid drives sea level
            sea_level_condition = 60 * np.sin(self.times / (5000 * np.pi))
        self.schematic_fig = plt.figure(figsize=(14,6)) #instantiate figure
        self.schematic = plt.subplot()
        plt.gcf().subplots_adjust(bottom=0.20)
        plt.xlabel('Distance from Subduction Zone [km]')
        plt.ylabel('Elevation above datum [m]')
        plt.ion()
        plt.show()
        
        return sea_level_condition
            
    def run(self, sea_level_condition, max_deflection, conv_rate, t_plot, dt, max_up_growth_rate, surf_light_intens, extinct_coeff, saturat_light_intens):
        x = self.pos_line_load - self.subducting_plate_x #distance from load to point
        x_for_plotting = np.arange(x_min/1000, x_max/1000, dx/1000)
        deflection = coral_model.flex(max_deflection, x, self.flex_param)
        bedrock_height = -deflection
        msl = max(bedrock_height-100)
        run_count = 0
        for t in range(len(self.times)-1):
            run_count +=1 #keep track of iterations
            current_time = self.times[t]
            rel_sl = msl + sea_level_condition[run_count - 1]
            x = abs(self.pos_line_load - self.subducting_plate_x) #distance from load to point
            deflection = coral_model.flex(max_deflection, x, self.flex_param)
            self.bedrock_height = -deflection #adjust bedrock height to account for deflection
            depth = rel_sl - (self.bedrock_height + self.coral_height)
            self.up_growth_rate = coral_model.grow_coral(depth, self.up_growth_rate, max_up_growth_rate, surf_light_intens, extinct_coeff, saturat_light_intens)
            self.coral_height += (self.up_growth_rate * dt) #THIS IS THE SOLVER STEP!!!
            self.pos_line_load -= (conv_rate * dt) #advance line load to the "left"
            plate_for_plot = self.subducting_plate_x[self.subducting_plate_x <= self.pos_line_load]
            x_for_plotting = x_for_plotting[self.subducting_plate_x/1000 <= self.pos_line_load/1000]
            if current_time % t_plot == 0: #do plotting and filling 
                self.schematic.clear()
                self.schematic.plot(x_for_plotting, self.bedrock_height[0:len(plate_for_plot)], color='.5', label='Bedrock')
                self.schematic.plot(x_for_plotting, self.bedrock_height[0:len(plate_for_plot)]+self.coral_height[0:len(plate_for_plot)], color='c', label='Coral')
                self.schematic.plot((min(x_for_plotting), max(self.subducting_plate_x)), (rel_sl, rel_sl), color='b', label='Water')
                self.schematic.fill_between(x_for_plotting, -4000, self.bedrock_height[0:len(plate_for_plot)], facecolor='.5', interpolate=True)
                self.schematic.fill_between(self.subducting_plate_x/1000, self.bedrock_height+self.coral_height, rel_sl, where=rel_sl >self.bedrock_height+self.coral_height, facecolor='b', interpolate=True)
                self.schematic.fill_between(x_for_plotting, self.bedrock_height[0:len(plate_for_plot)], self.bedrock_height[0:len(plate_for_plot)]+self.coral_height[0:len(plate_for_plot)], where=self.bedrock_height[0:len(plate_for_plot)]+self.coral_height[0:len(plate_for_plot)] >self.bedrock_height[0:len(plate_for_plot)], facecolor='c', interpolate=True)
                self.schematic.set_ylim(-4000, 800)
                self.schematic.set_xlim(0, 200)
                plt.title('The Stegosaurus Coral Model')
                plt.xlabel('Distance from Subduction Zone [km]')
                plt.ylabel('Elevation above datum [m]')
                plt.text(5, -3500, 'Time [yrs]: %.1f' % current_time)
                plt.pause(0.2)
            else:
                pass
            
    def finalize(self):
        self.schematic_fig.savefig('shobe_corals.png', bbox_inches='tight')
        
if __name__ == "__main__": 
    isotope_or_sinusoid = 0 #0 for isotope record, 1 for sinusoid
    x_min = 0 #m
    x_max = 200000 #m
    dx = 100 #m
    t_min = 0 #years
    t_max = 600000 #years
    dt = 10 #yrs
    t_plot = 10000#years
    init_pos_line_load = 250000 #m, end of x domain
    conv_rate = 0.1 #convergence rate, m/yr
    max_deflec = 10000 #max deflection, m
    g = 9.81 #ms-2
    flex_rigid = 1e22 #flexural rigidity, Nm
    dens_mantle = 3300 #kg/m3
    dens_ocean = 1030 #kg/m3
    max_up_growth_rate = .012 #m/yr
    extinct_coeff = .16 #/m
    surf_light_intens = 2250 #uE m-2 s-1
    saturat_light_intens = 450 #uE m-2 s-1
    sl_file_name = 'Pleist_oxy.txt'
    #instantiate classes
    coral_model = Coral()
    #INITIALIZE
    sl_condition = coral_model.initialize(x_min, x_max, dx, t_min, t_max, dt, init_pos_line_load, g, flex_rigid, dens_mantle, dens_ocean, sl_file_name)
    #RUN
    coral_model.run(sl_condition, max_deflec, conv_rate, t_plot, dt, max_up_growth_rate, surf_light_intens, extinct_coeff, saturat_light_intens)
    #FINALIZE
    coral_model.finalize()