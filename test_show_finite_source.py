# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:39:42 2016

@author: cjburke
Test of MEPL to illustrate images of a finite source
"""
import microlens as mepl
import matplotlib.pyplot as plt
import numpy as np


# Specify lens properties
lens_separation = 0.75
lens_mass1 = 1.0
lens_mass2 = 0.01

# Specify source center and radius
source_x = -0.02
source_y = -0.01
source_rad = 0.01


# Instantiate an empty lens object
lens = mepl.lens_properties()
# Calculate lens properties
# lens 1 at origin and lens 2 on x-axis
lens.set_lens_properties(0.0, 0.0, lens_mass1, \
                         lens_separation, 0.0, lens_mass2)
                         
# Instantiate an empty time stamp object
time_step = mepl.single_time_step()
# Fill in time_step object
time_step.centerx = source_x
time_step.centery = source_y
time_step.radius = source_rad
time_step.srcn = 13
time_step.lens = lens
                     
# Calculate finite image shapes and magnification
# debugLevel ==4 shows image shapes with interactive graphics
curmag, magtyp = mepl.get_time_step_mag(time_step, debugLevel=4)
print "Total Mag: {0:7.3f}".format(curmag)                    

