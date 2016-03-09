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

MOTIVATION

With the upcoming K2 mission and campaign 9 along with WFIRST moving forward, I was intrigued to learn more about the microlensing extrasolar planet detection technique.  There appeared to be a lack of freely available code to play around with making pretty pictures of the caustics, finite source effects, and magnification light curves.  After reading the literature a lot I thought I had a pretty good idea of how to program up the binary lens scenario in the point source approximation.  I successfully did that without too much effort, and thought why not try to model finite source effects.  That ended up being a totally different ballgame, but I got far enough along to stick it through and decided to release MEPL for others to get started and hopefully improve and make more comprehensive.  I built the code from scratch with the main goals being an open source modern python code with an emphasis on making the code amenable to parallelization.  The current code does not do any data fitting!  It currently provides magnification light curves.  See next steps below.

REFERENCES

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

Michael Miller (2013) http://researcharchive.vuw.ac.nz/xmlui/bitstream/handle/10063/2775/thesis.pdf?sequence=2
Has a modern and extremely helpful explicit writing out the polynomial coefficients.

Miller's work is based upon earlier works that I also found helpful

Rhie, S (1997) http://adsabs.harvard.edu/abs/1997ApJ...484...63R

Gaudi & Gould (1997) http://adsabs.harvard.edu/abs/1997ApJ...484...63R

Rhie, S (2002) http://adsabs.harvard.edu/abs/2002astro.ph..2294R

Bennett & Rhie (1996) http://adsabs.harvard.edu/abs/2002astro.ph..2294R

Also helpful is Paul Chote's practical descriptions of the problem
http://researcharchive.vuw.ac.nz/bitstream/handle/10063/1890/thesis.pdf?sequence=1

Reading the above works allowed me to get the point source approximation with a binary lens model coded up.
For treating finite source effects, I adopted they polygon or outline of the source approach.  The  uniform brightness finite image areas are calculated using line/contour integrals around the image boundaries.  This seemed intially more amenable to parallelization than the approaches taken by Bennett et al.  I believe contour integral approach is originaly attributed to 

Gould & Gaucherel (1997) http://adsabs.harvard.edu/abs/1997ApJ...477..580G  in the case of binary lenses

Bozza, V. (2010) http://adsabs.harvard.edu/abs/2010MNRAS.408.2188B  Takes this approach to higher order approximation.  I have yet to implement the treatment from here, but this paper is also a good read for the contour integral approach for finite sources.

PERFORMANCE

There are two avenues to improve performance.  From profiling most time is spent in solve_binary_roots() in particular the call to numpy polyroots where the 5th order polynomial roots are determined.  One could try compiled root solving codes and interface them with MEPL.  There is Skowron & Gould (2012) open source code, but it seems tailored for binary lens, I dont know if it would support triple lens scenario should MEPL be extended to support triple lenses.  So maybe a more general compiled root finding code could be tried.

The other avenue is to reduce the number of source vertices needed.  The fill_image_gap() functions job is to determine when large stretching of the source image was occurring and add more source vertices where needed.  The fill_image_gap was tuned in order to get sufficient resolution such that the subsequent call to finite_source_mag() would successfully be able to trace the image loops in the complicated case of the finite source strattling caustics.  Gould & Gaucherel (1997) show that actually tracing the image loops isnt necessary.  Thus, in principal one could forgo the complicated loop matching and deal with end points following Gould & Gaucherl.  This may reduce the requisite number of source vertices needed just to fullfill the loop matching code and one can then focus on a tuning fill_image_gaps() in order to achieve a desired area precision rather the current tuning of fill_image_gaps() which is tuned to make sure finite_source_mag() runs successfully.  Although I am sure that the area integral precision and an ability to trace the loops are highly correlated.  Another way to improve area integral precision is to go to a higher order curve integration outlined by Bozza (2010).  There is some precursor code commented out and untested to calculate the requisite kappa derivatives needed for Bozza (2010).  

The parllelization is currently implented at the time stamp level.  That provides significant speed up for multi core computers.  I originally developed the code for supporting parallelization at the individual source vertex solve_binary_image() level.  Even though most time is spent here each vertex calculation is very short, its just a time sink because it is called many many times.  The result of parallizing using pytons pool.map appproach at this level was very disappointing because the parllelization overheads were comparable to the indivdual function call resulting in slowing down of the runtime.  However, allowing each core to work on a single time step / source position works really well for speedups.  On a single core machine with large 0.1 source radii in particularly complicated caustic crossing scenarios it can take up to 40secs to run the ipython microlens.py test code with the 3000 source position time stamps.  Most of the time for typical source radii the single core run time is of order 10 seconds for 3000 source position times.  You will find run times on large multi-core computer signficantly shorter.
  
WHAT TO WORK ON NEXT

-Implementing a fitter.  That is a large undertaking, but would be needed to actually make the code useful for scientific purposes.

-Runtimes can be improved.  The lengthy discussion above points out the easiest ways to improve runtime outside of porting all the code to a compiled language.

-The test_point_approx() function which determines whether to use the point source approximation or finite source is very rudimentary and sometimes results in larger than desired jumps at the transitions.  That could be updated to something smarter.

-Benchmarking against other codes.  I have only eyeballed the results seem consistent with published figures.  I have not rigorously checked the validity of my code against others.  Anyone with software to compare to I would love hear what you find.

-The code does not have any error checking or rigorous magnification precision monitoring.  In other words the code doesnt know if it is messing up.

-There is a scenario in trace_three_images() that is not implemented yet, but I have only encountered a crash in the unimplemented routine once in testing and wasn't able to recreate it to debug it.  This is a rare occurrence, but if the code does crash it will likely be here.

-Implementing limb darkened sources.  Bozza (2010) describes how to implement limb darkened sources with the current polygon approach.  I wonder if it would actually be faster to switch to a 2D integral approach ala Bennett et al. once the edges of the source are outlined in the image plane, rather then continue to the polygon approach with multiple smaller source outlines.

-Support triple lenses.  The Miller and Chote work in principle provides the nitty gritty details of implementing point source approximation for a triple lens.  I would definitely try and get the Gould & Gaucherel approach to the contour integral working that avoids the nauseating bookkeeping done in finite_source_mag() for sources strattling caustics before trying to extend the code to handle triple lens.
