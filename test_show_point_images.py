# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:39:42 2016

@author: cjburke
Test of MEPL to illustrate images of a point source
"""
import microlens as mepl
import matplotlib.pyplot as plt
import numpy as np


# Specify lens properties
lens_separation = 0.75
lens_mass1 = 1.0
lens_mass2 = 0.01

# Specify source center
source_x = 0.01
source_y = -0.01

# Instantiate an empty lens object
lens = mepl.lens_properties()
# Calculate lens properties
# lens 1 at origin and lens 2 on x-axis
lens.set_lens_properties(0.0, 0.0, lens_mass1, \
                         lens_separation, 0.0, lens_mass2)
# Instantiate empty source vertex
src_obj = mepl.source_plane_vertex()
# set source position
#  for point source the radius = 0.0 and angle = 0.0
src_obj.set_source_position(source_x, source_y, 0.0, 0.0)
# Solve for image position
imageLocs, jacobians, paritys, kappas, kappaDerivs = \
                    mepl.solve_binary_roots(src_obj, lens)
# Calculate magnification of the point source approximation
magnifs = 1.0/np.abs(jacobians)
magnifTotal = np.sum(magnifs)
tmpx = np.real(imageLocs)
tmpy = np.imag(imageLocs)
for i in range(len(imageLocs)):
    print "Image {0:d} x: {1:6.4f} y: {2:6.4f} Mag: {3:7.3f}".format( \
                 i, tmpx[i], tmpy[i], magnifs[i])
print "Total Mag: {0:7.3f}".format(magnifTotal)                    

myorange=tuple(np.array([230.0,159.0,0.0])/255.0)
myskyblue=tuple(np.array([86.0,180.0,233.0])/255.0)
# plot lens1 and lens 2 as black circle
plt.plot([0.0, lens_separation],[0.0, 0.0], 'ok')
# plot source as orange circle
plt.plot([source_x], [source_y], 'o', color=myorange)
# plot images as blue circle
plt.plot(tmpx, tmpy, 'o', color=myskyblue)
plt.show()


