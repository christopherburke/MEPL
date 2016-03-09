# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:39:42 2016

@author: cjburke
Test of MEPL to illustrate making a figure of caustics
"""
import microlens as mepl

# Specify lens properties
lens_separation = 0.75
lens_mass1 = 1.0
lens_mass2 = 0.01
# Instantiate an empty lens object
lens = mepl.lens_properties()
# Calculate lens properties
# lens 1 at origin and lens 2 on x-axis
# set debugLevel == 4 in order to show figure of caustic
#  with interactive graphics
lens.set_lens_properties(0.0, 0.0, lens_mass1, \
                         lens_separation, 0.0, lens_mass2, debugLevel=4)
                    


