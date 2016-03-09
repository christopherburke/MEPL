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

Motivation

With the upcoming K2 mission and campaign 9 along with WFIRST moving forward, I was intrigued to learn more about the microlensing extrasolar planet detection technique.  There appeared to be a lack of freely available code to play around with making pretty pictures of the caustics, finite source effects, and magnification light curves.  After reading the literature a lot I thought I had a pretty good idea of how to program up the binary lens scenario in the point source approximation.  I successfully did that without too much effort, and thought why not try to model finite source effects.  That ended up being a totally different ballgame, but I got far enough along to stick it through and decided to release MEPL for others to get started and hopefully improve and make more comprehensive.  I built the code from scratch with the main goals being an open source modern python code with an emphasis on making the code amenable to parallelization.

References

The programming itself is original, but I rely very heavily on descriptions and equations from many sources from the astronomy literature.  Here is a list of some of the more influential articles in building the code and helping me understand gravitational lensing. 
In thinking about microlensing I had sort of an aha moment from thinking about the analogy of thin lenses in optics.
The goals are similar, determine how a source looks after passing through an optical lens.   I stumbled upon this web page that opened lots of doors into my thinking. 
http://www.artemis-uk.org/Microlensing_physmath.html
By the same analogy I loved to see this lab experiment of the wine glass foot concept being used in the lab
https://vela.astro.ulg.ac.be/themes/extragal/gravlens/bibdat/engl/DE/didac.html
Apparently one can buy a plexiglass 'microlens' simulator as used in their experiments!
http://sciences.ulg.ac.be/publications/lentilles/
Anyone in in the EU want to mail me one?
The equation that describes going from the source plane to the image plane is typically described as solving a 5th order polynomial in the case of a binary lens.  When using complex numbers the lens equation becomes a simple form that I believe is attributed to Witt (1990) http://adsabs.harvard.edu/abs/1990A%26A...236..311W
The following series of papers were helpful in actually seeing how the lens equation is transformed into a 5th order polynomial in practice.

Rhie, S (1997) http://adsabs.harvard.edu/abs/1997ApJ...484...63R

Gaudi & Gould (1997) http://adsabs.harvard.edu/abs/1997ApJ...484...63R

Rhie, S (2002) http://adsabs.harvard.edu/abs/2002astro.ph..2294R

Bennett & Rhie (1996) http://adsabs.harvard.edu/abs/2002astro.ph..2294R

