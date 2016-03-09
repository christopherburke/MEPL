# MEPL
Micro ExoPlanet Lens (MEPL - pronouced meeple) code
An open source python toolbox to explore the microlensing signal from an extrasolar planet.  Currently it supports binary lens scenarios with a finite source of uniform brightness.

Getting Started
Test scripts are available to ease one into the code capabilities.  I recommend one reads and runs the following in this order

1) ipython test_show_caustics.py - Start with binary lens properties and show a figure with the source plane caustics

2) ipython test_show_point_images.py - For a given binary lens and point source position, solve for the image locations, and show       a figure of the lens positions (black), source position (orange), and 5 image positions (blue).  There are 5 images because       the source is inside the caustic curve.  Using test_show_caustics.py, specify a source position outside of the caustics
      such that you only get 3 images.
      
3) ipython test_show_finite_source.py - For a given binary lens and finite point source, show the images formed.

4) ipython microlens.py  - Full featured test of a randomly selected binary lens, source trajectory that passes through 
                            caustics, and uses multiprocessing to parallelize the calculation


