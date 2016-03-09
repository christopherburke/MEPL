""" Module for binary microlens modeling
    Author: Christopher J Burke 
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import numpy.polynomial.polynomial as ply
import time
from operator import attrgetter
import multiprocessing
from meeple_shape import meeple_outline

class Timer(object):
    """Useful timer object used in solve_binary_roots
    """
    def __init__(self, verbose=False):
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        return self

    def __exit__(self, *args):
        self.end = time.time()
        self.secs = self.end - self.start
        self.msecs = self.secs * 1000  # millisecs
        if self.verbose:
            print 'elapsed time: %f ms' % self.msecs

class source_plane_vertex:
    """Storage class to keep parameters of a single source_vertex
       CONTENTS
         centerx  - [float] source center x-axis position 
         centery - [float] source center y-axis position
         radius - [float] source radius
         angle - [Radians] source vertex angle position measured counter
                           clockwise from x-axis
         x - [float] source vertex x-axis position centerx + radius*Cos(angle)
         y - [float] source vertex y-axis position centery + radius*Sin(angle)
         w - [complex] complex form of source plane vertex x and y positions x + y*1j
         stretchFactor - [int] When augmenting the source vertices
                             add this many new source vertices
                            in the gap between this and its neighbor
         nImages - [Int 3 or 5] number of images in image plane.  For binary lens
                     the number of images is 3 or 5
         imageLocs - [complex] complex image plane vertex location there are nImages
         imageParity - [1 / -1] Parity of image vertex
         jacobians - [complex] Surface brightness jacobian determinant at vertex location
         kappas - [complex] auxillary variable for vertex location
         kappaDerivs - [complex] auxillary variable vertex location
       METHODS
         set_source_position() -  From the source center position, radius, and angle
            calculate the x and y position of the circular source outline
         sort() -  Used to sort a list of source objects by their angle
            This is important because the routines that trace image
            boundaries assume they come in in order
    """
    def __init__(self):
        self.centerx = 0.0
        self.centery = 0.0
        self.radius = 1.0e-2
        self.angle = 0.0
        self.x = 0.0
        self.y = 0.0
        self.w = 0.0 + 1j
        self.stretchFactor = 1
        self.nImages = 0
        self.imageLocs = np.array([], dtype=np.complex)
        self.imageParity = np.array([])
        self.jacobians = np.array([], dtype=np.complex)
        self.kappas = np.array([], dtype=np.complex)
        self.kappaDerivs = np.array([], dtype=np.complex)

    def set_source_position(self, centerx, centery, radius, angle):
        self.centerx = centerx
        self.centery = centery
        self.radius = radius
        self.angle = angle
        # Normal circle source shape
        self.x = centerx + radius * np.cos(angle)
        self.y = centery + radius * np.sin(angle)
        # Use hidden easter egg meeple shape
        #if radius > 0.0:
        #    tmpx, tmpy = meeple_outline(angle, radius)
        #    self.x = (centerx + tmpx)[0]
        #    self.y = (centery + tmpy)[0]
        
        self.w = self.x + self.y*1j

    def sort(self, a, b):
        a.angle < b.angle
        
    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

class lens_properties:
    """Storage class for the lens properties 
       CONTENTS
       x1 - [float] x-axis lens 1 position
       y1 - [float] y-axis lens 1 position 
       m1 - [float] lens 1 mass
       x2 - [float] x-axis lens 2 position
       y2 - [float] y-axis lens 2 position
       m2 - [float] lens 2 mass
       mtot - [float] summed lens mass
       r1 - [complex] complex lens 1 position
       r2 - [complex] complex lens 2 position
       e1 - [float] m1 / mtot
       e2 - [float] m2 / mtot
       s - [float] lens image plane separation
       q - [float] lens mass ratio
       GBig, HBig, and eta lens auxillary values used in the polynomial
            function coefficients.  These dont depend on source properties
       caustLocs [complex] The source plane caustics get set by
               set_lens_properties.  These are either an analytic approximation
               for the cusp locations or full caustic curve locations found 
               numerically.
       criticLocs [complex] critical curves in the source plane
                    These only get populated when the full caustic curve is
                    found numerically.  Check non-empty before using
       centCaustSize [float] - The central caustic size
       centCaustLoc [float] - The central caustic location
       METHODS
         set_lens_properties() - Sets up all aspects of the problem that depend upon the lens
           configuration independent of the source position
           This only needs to be done once for a lens that is stationary
           Calculates auxillary arrays GBig, HBig, and eta as well as
            determining caustic and crictical curves
    """
    def __init__(self):
        self.x1 = 0.0
        self.y1 = 0.0
        self.m1 = 1.0
        self.x2 = 0.5
        self.y2 = 0.0
        self.m2 = 0.1
        self.mtot = 1.1
        self.r1 = 0.0 + 0j
        self.r2 = 0.5 + 0j
        self.e1 = 1.0 / 1.1
        self.e2 = 0.1 / 1.1
        self.s = 0.5
        self.q = 0.1
        self.GBig = np.zeros((3, 1), dtype=np.complex)
        self.HBig = np.zeros((3, 1), dtype=np.complex)
        self.eta = np.zeros((3, 5), dtype=np.complex)
        self.caustLocs = np.array([], dtype=np.complex)
        self.criticLocs = np.array([], dtype=np.complex)
        self.centCaustSize = 1.0
        self.centCaustLoc = 0.0 + 0j
        
    def set_lens_properties(self, x1, y1, m1, x2, y2, m2, debugLevel = 0):
        self.mtot = m1 + m2
        self.r1 = x1 + y1*1j
        self.r2 = x2 + y2*1j
        self.e1 = m1 / self.mtot
        self.e2 = m2 / self.mtot
        self.s = np.abs(self.r1 - self.r2)
        self.q = self.e2 / self.e1
        self.GBig = np.array([-self.e1*self.r2 - self.e2*self.r1, 1.0, 0.0])
        self.HBig = np.array([self.r1*self.r2, -self.r1 - self.r2, 1.0])
        self.eta[0,0] = self.HBig[0] * self.HBig[0]
        self.eta[0,1] = 2.0 * self.HBig[0] * self.HBig[1]
        self.eta[0,2] = self.HBig[1] * self.HBig[1] + 2.0 * self.HBig[0] * self.HBig[2]
        self.eta[0,3] = 2.0 * self.HBig[1] * self.HBig[2]
        self.eta[0,4] = self.HBig[2] * self.HBig[2]
        self.eta[1,0] = self.GBig[0] * self.HBig[0]
        self.eta[1,1] = self.GBig[0] * self.HBig[1] + self.GBig[1] * self.HBig[0]
        self.eta[1,2] = self.GBig[1] * self.HBig[1] + self.GBig[0] * self.HBig[2] + self.HBig[0] * self.GBig[2]
        self.eta[1,3] = self.GBig[1] * self.HBig[2] + self.GBig[2] *self. HBig[1]
        self.eta[1,4] = self.GBig[2] * self.HBig[2]
        self.eta[2,0] = self.GBig[0] * self.GBig[0]
        self.eta[2,1] = 2.0 * self.GBig[0] * self.GBig[1]
        self.eta[2,2] = self.GBig[1] * self.GBig[1] + 2.0 * self.GBig[0] * self.GBig[2]
        self.eta[2,3] = 2.0 * self.GBig[1] * self.GBig[2]
        self.eta[2,4] = self.GBig[2] * self.GBig[2]
        # Get the caustic curves analytically or numerically
        # Use analytic cusp locations where they are valid
        if (self.q <= 5.0e-3 and np.abs(np.log10(self.s)) > 0.1):
            self.caustLocs = get_cusp_locations(self.q, self.s)
        else:
            # Get numerical caustic points
            # Get critical curve
            theta = np.linspace(0.0, 2.0 * np.pi, 25)
            for i in range(len(theta)):    
                zlocs = solve_binary_critcurve(theta[i], self.r1, self.r2, \
                            self.e1, self.e2)
                wlocs = zlocs - (self.e1 / (np.conj(zlocs) - np.conj(self.r1))) - \
                       (self.e2 / (np.conj(zlocs) - np.conj(self.r2)))
                self.caustLocs = np.append(self.caustLocs, wlocs)
                self.criticLocs = np.append(self.criticLocs, zlocs)
        # Determine size of central caustic
        idx = np.where(np.abs(self.caustLocs) < 0.3)[0]
        tmpx = np.real(self.caustLocs)[idx]
        tmpy = np.imag(self.caustLocs)[idx]
        # Get the distance between all pairs of image positions
        dist1d = dist.pdist(np.hstack((tmpx[:,None], \
                    tmpy[:,None])), 'euclidean')
        self.centCaustSize = max(dist1d)
        self.centCaustLoc = np.mean(self.caustLocs[idx])
        if debugLevel > 1:
            fig, ax = plt.subplots()
            ax.plot(np.real(self.caustLocs), np.imag(self.caustLocs), '+r')
            #plt.xlim([-0.2, 0.2])
            #plt.ylim([-0.2, 0.2])
            if debugLevel == 2:
                plt.pause(0.5)
            if debugLevel == 3:
                ax.set_title('Noninteractive Window, click on window to close and proceed:')
                plt.waitforbuttonpress()
            if debugLevel == 4:
                fig.canvas.mpl_connect('key_press_event', on_key_event)    
                # Disable default Matplotlib shortcut keys:
                keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
                for key in keymaps:
                    plt.rcParams[key] = ''
                ax.set_title('Interactive Window, select a vowel on the keyboard to close window and proceed:')
                plt.show()
                
                
            plt.close()
            print "Plotting in lens_properties"

    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''
        
class single_time_step:
    """Storage class for a single time/observations
         that you want to calculate the microlens magnification for
       CONTENTS
         centerx - [float] source x position
         centery - [float] source y position
         radius - [float] source radius in units of einstein radius
         srcn - [int] Number of base source vertices to start with
                      The runtime is very sensitive to this parameter
                      srcn=13 is what has worked in testing as stable while
                        providing decent runtime.  If you run into troubles
                        turn on debugging to watch images and bump this value up
         lens - This needs to be a valid lens_properties class object that
                 is filled in before 
         debugLevel - [0,1,2,3,4] 0 no output, 1 printed output diagnostics
                     2 figures plotted with minimal pausing
                     3 non interactive figures plotted with requiring click on window to proceed
                     4 Some figures plotted that can be interacted/zoomed,etc. with
                        pressing a vowel key is needed to proceed
    """
    def __init__(self):
        self.centerx = 0.0
        self.centery = 0.0
        self.radius = 1.0e-2
        self.srcn = 13
        self.lens = None
        self.debugLevel = 0
       
    def __str__(self):
        for k in self.__dict__:
            print k, self.__dict__[k]
        return ''

def get_cusp_locations(q, s):
    """ Han, C (2006) http://adsabs.harvard.edu/abs/2006ApJ...638.1080H
        and references therein
        gives quick analytic equations for planetary and central caustic
        locations.  This returns the complex cusp locations
    """
    sp2 = s * s
    qphalf = np.sqrt(q)
    # Get central caustic cusp locations
    cosphi = 0.75*(s + 1.0/s)*(1.0 - np.sqrt(1.0 - 32.0/9.0/(s + 1.0/s)**2))
    sinphi = np.sqrt(1.0 - cosphi*cosphi)
    a = 2.0 * q * np.abs(sinphi)**3 / (s + 1.0/s - 2.0*cosphi)**2
    all_cusp = np.array([q/(1.0+s)/(1.0+1.0/s) + 0j, \
                    -q/(1.0-s)/(1.0-1.0/s) + 0j, \
                    a*1j, -a*1j])
    # Get planetary caustic cusp locations
    #  for s ~ 1 the planetary and central caustic have merged just skip
    #   their calculation
    if s > 1.1:
        cen_plan_pos = s - 1.0/s + 0j
        del_12 = 2.0 * qphalf / s / np.sqrt(sp2 - 1.0) + 0j
        del_34 = 2.0 * qphalf / s / np.sqrt(sp2 + 1.0)*1j
        cusp_plan_pos = np.array([cen_plan_pos + del_12, cen_plan_pos - del_12, \
                            cen_plan_pos + del_34, cen_plan_pos - del_34])
        all_cusp = np.append(all_cusp, cusp_plan_pos)

    if s < 0.9:
        cen_plan_pos = (s - 1.0/s) + 0j
        del_1 = 2.0 * qphalf / s / np.sqrt(1.0 + sp2) *1j
        #del_2 = 2.0 * qphalf * np.sqrt(1 - sp2) / s *1j
        a = 0.5 * np.arcsin(np.sqrt(3.0) * 0.5 * sp2)
        #th0p = np.pi / 2.0 + a
        th0n = np.pi / 2.0 - a
        #b0p = np.sqrt(sp2*sp2 - np.sin(2.0*th0p)**2)
        b0n = np.sqrt(sp2*sp2 - np.sin(2.0*th0n)**2)
        varc = sp2 - 1.0 / sp2
        #k0p = np.sqrt((np.cos(2.0 * th0p) + b0p) / varc)
        k0n = np.sqrt((np.cos(2.0 * th0n) - b0n) / varc)
        cusp_plan_pos = np.array([cen_plan_pos + del_1, cen_plan_pos - del_1, \
#                        cen_plan_pos + del_2, cen_plan_pos - del_2])
                            cen_plan_pos + \
                            qphalf*(k0n - 1.0/k0n + k0n/sp2)*np.cos(th0n) + \
                            qphalf*(k0n - 1.0/k0n - k0n/sp2)*np.sin(th0n)*1j, \
                            cen_plan_pos + \
                           -qphalf*(k0n - 1.0/k0n + k0n/sp2)*np.cos(th0n) + \
                            qphalf*(k0n - 1.0/k0n - k0n/sp2)*np.sin(th0n)*1j, \
                            cen_plan_pos + \
                            -qphalf*(k0n - 1.0/k0n + k0n/sp2)*np.cos(th0n) + \
                            -qphalf*(k0n - 1.0/k0n - k0n/sp2)*np.sin(th0n)*1j, \
                            cen_plan_pos + \
                            qphalf*(k0n - 1.0/k0n + k0n/sp2)*np.cos(th0n) + \
                            -qphalf*(k0n - 1.0/k0n - k0n/sp2)*np.sin(th0n)*1j])
        all_cusp = np.append(all_cusp, cusp_plan_pos)
    return all_cusp

def solve_binary_roots(obj, lens):
    """Determine the image locations for a given source vertex by solving the
       roots of the fifth order polynomial.
       This is the work horse of the code and where all the time is spent.
       In particular, it is spent in the call to polyroots.
       The polynomial coefficient expressions come from Miller M.(2013)
        http://researcharchive.vuw.ac.nz/xmlui/handle/10063/2775
       which come from Rhie, S (2002)
       http://adsabs.harvard.edu/abs/2002astro.ph..2294R
       Timing has shown a small performance increase by explicitly writing
         out some of the arrays that are used.  Thus,
         compact expressions are commented out while the expanded
         explicit expressions are used.
       INPUTS:
          obj - source_plane_vertex object 
          lens - lens_properties object
       OUTPUTS:
          zlocs - [complex] numpy array of length 3 or 5 image locations
          jacs - jacobian determinant at image locations
          paritys - [1,-1] gives parity for each image location
          zkappas - related to jacobian
          zkappaPs - kappa derivative potentially useful in higher order
                     area integral.  Currently deactivated
       
    """
    #with Timer() as t:
    # Keep indent to use the timer
    if True:
        wbar = np.conj(obj.w)
        r1bar = np.conj(lens.r1)
        r2bar = np.conj(lens.r2)
        wbar1 = r1bar - wbar
        wbar2 = r2bar - wbar
#        XBig = np.array([wbar1*wbar2, -wbar1 - wbar2, 1.0])
        XBig = np.zeros((3,), dtype=np.complex)
        XBig[0] = wbar1*wbar2
        XBig[1] = -wbar1 - wbar2
        XBig[2] = 1.0
#        VBig = np.array([-lens.e1*wbar2 - lens.e2*wbar1, 1.0, 0.0])
        VBig = np.zeros((3,), dtype=np.complex)
        VBig[0] = -lens.e1*wbar2 - lens.e2*wbar1
        VBig[1] = 1.0
#        WBig = np.array([obj.w*XBig[0] + VBig[0], obj.w*XBig[1] + VBig[1], obj.w])
        WBig = np.zeros((3,), dtype=np.complex)
        WBig[0] = obj.w*XBig[0] + VBig[0]
        WBig[1] = obj.w*XBig[1] + VBig[1]
        WBig[2] = obj.w
#        etaX = np.reshape(np.repeat(XBig, 5), (3,5))
        etaX = np.zeros((3,5), dtype=np.complex)
        etaX[0,0] = XBig[0]
        etaX[1,0] = XBig[1]
        etaX[2,0] = XBig[2]
        etaX[0,1] = XBig[0]
        etaX[1,1] = XBig[1]
        etaX[2,1] = XBig[2]
        etaX[0,2] = XBig[0]
        etaX[1,2] = XBig[1]
        etaX[2,2] = XBig[2]
        etaX[0,3] = XBig[0]
        etaX[1,3] = XBig[1]
        etaX[2,3] = XBig[2]
        etaX[0,4] = XBig[0]
        etaX[1,4] = XBig[1]
        etaX[2,4] = XBig[2]        
#        etaW = np.reshape(np.repeat(WBig, 5), (3,5))
        etaW = np.zeros((3,5), dtype=np.complex)
        etaW[0,0] = WBig[0]
        etaW[1,0] = WBig[1]
        etaW[2,0] = WBig[2]
        etaW[0,1] = WBig[0]
        etaW[1,1] = WBig[1]
        etaW[2,1] = WBig[2]
        etaW[0,2] = WBig[0]
        etaW[1,2] = WBig[1]
        etaW[2,2] = WBig[2]
        etaW[0,3] = WBig[0]
        etaW[1,3] = WBig[1]
        etaW[2,3] = WBig[2]
        etaW[0,4] = WBig[0]
        etaW[1,4] = WBig[1]
        etaW[2,4] = WBig[2]        
        etaX = lens.eta * etaX
        etaW = lens.eta * etaW
#        c0 = -np.sum(etaW[:,0])
        c0 = -(etaW[0,0] + etaW[1,0] + etaW[2,0])
#        c1 = np.sum(etaX[:,0]) - np.sum(etaW[:,1])
        c1 = etaX[0,0] + etaX[1,0] + etaX[2,0] - (etaW[0,1] + etaW[1,1] + etaW[2,1])
#        c2 = np.sum(etaX[:,1]) - np.sum(etaW[:,2])
        c2 = etaX[0,1] + etaX[1,1] + etaX[2,1] - (etaW[0,2] + etaW[1,2] + etaW[2,2])        
#        c3 = np.sum(etaX[:,2]) - np.sum(etaW[:,3])
        c3 = etaX[0,2] + etaX[1,2] + etaX[2,2] - (etaW[0,3] + etaW[1,3] + etaW[2,3])
#        c4 = np.sum(etaX[:,3]) - np.sum(etaW[:,4])
        c4 = etaX[0,3] + etaX[1,3] + etaX[2,3] - (etaW[0,4] + etaW[1,4] + etaW[2,4])
#        c5 = np.sum(etaX[:,4])
        c5 = etaX[0,4] + etaX[1,4] + etaX[2,4]
    #print "=> elasped a: %s s" % t.secs
    #with Timer() as t:    
        # This is numpy version
        zlocs = ply.polyroots((c0, c1, c2, c3, c4, c5))
    #print "=> elasped b: %s s" % t.secs
    #with Timer() as t:
        # return the image locations that are consistent with lens equation
        zwdiff = np.abs(zlocs - (lens.e1 / (np.conj(zlocs) - np.conj(lens.r1))) - \
                        (lens.e2 / (np.conj(zlocs) - np.conj(lens.r2))) - obj.w)
        idx = np.where(zwdiff < 1.0e-4)[0]
        zlocs = zlocs[idx]
        # Now get magnifications
    #    z0 = np.full_like(zlocs, 1.0)
        z2 = zlocs * zlocs
        nz = len(zlocs)
        #zall = np.reshape(np.append(z0, np.append(zlocs, z2)),(3,nz))
    #    zall = np.vstack((z0, zlocs, z2))
        GBigP = np.array([1.0, 0.0])
        HBigP = np.array([-lens.r1 -lens.r2, 2.0])
        
        if nz == 3:
            zall = np.zeros((3,3), dtype=np.complex)
            zall[0,0]=1.0
            zall[0,1]=1.0
            zall[0,2]=1.0
            zall[1,0]=zlocs[0]
            zall[1,1]=zlocs[1]
            zall[1,2]=zlocs[2]
            zall[2,0]=z2[0]
            zall[2,1]=z2[1]
            zall[2,2]=z2[2]
            HBig2d = np.zeros((3,3), dtype=np.complex)
            HBig2d[0,0] = lens.HBig[0]
            HBig2d[0,1] = lens.HBig[0]
            HBig2d[0,2] = lens.HBig[0]
            HBig2d[1,0] = lens.HBig[1]
            HBig2d[1,1] = lens.HBig[1]
            HBig2d[1,2] = lens.HBig[1]
            HBig2d[2,0] = lens.HBig[2]
            HBig2d[2,1] = lens.HBig[2]
            HBig2d[2,2] = lens.HBig[2]
            GBig2d = np.zeros((3,3), dtype=np.complex)
            GBig2d[0,0] = lens.GBig[0]
            GBig2d[0,1] = lens.GBig[0]
            GBig2d[0,2] = lens.GBig[0]
            GBig2d[1,0] = lens.GBig[1]
            GBig2d[1,1] = lens.GBig[1]
            GBig2d[1,2] = lens.GBig[1]
            GBig2d[2,0] = lens.GBig[2]
            GBig2d[2,1] = lens.GBig[2]
            GBig2d[2,2] = lens.GBig[2]
            HBigP2d = np.zeros((2,3), dtype=np.complex)
            HBigP2d[0,0] = HBigP[0]
            HBigP2d[0,1] = HBigP[0]
            HBigP2d[0,2] = HBigP[0]
            HBigP2d[1,0] = HBigP[1]
            HBigP2d[1,1] = HBigP[1]
            HBigP2d[1,2] = HBigP[1]
            GBigP2d = np.zeros((2,3), dtype=np.complex)
            GBigP2d[0,0] = GBigP[0]
            GBigP2d[0,1] = GBigP[0]
            GBigP2d[0,2] = GBigP[0]
            GBigP2d[1,0] = GBigP[1]
            GBigP2d[1,1] = GBigP[1]
            GBigP2d[1,2] = GBigP[1]
        else:
            zall = np.zeros((3,5), dtype=np.complex)
            zall[0,0]=1.0
            zall[0,1]=1.0
            zall[0,2]=1.0
            zall[0,3]=1.0
            zall[0,4]=1.0
            zall[1,0]=zlocs[0]
            zall[1,1]=zlocs[1]
            zall[1,2]=zlocs[2]
            zall[1,3]=zlocs[3]
            zall[1,4]=zlocs[4]
            zall[2,0]=z2[0]
            zall[2,1]=z2[1]
            zall[2,2]=z2[2]
            zall[2,3]=z2[3]
            zall[2,4]=z2[4]
            HBig2d = np.zeros((3,5), dtype=np.complex)
            HBig2d[0,0] = lens.HBig[0]
            HBig2d[0,1] = lens.HBig[0]
            HBig2d[0,2] = lens.HBig[0]
            HBig2d[0,3] = lens.HBig[0]
            HBig2d[0,4] = lens.HBig[0]
            HBig2d[1,0] = lens.HBig[1]
            HBig2d[1,1] = lens.HBig[1]
            HBig2d[1,2] = lens.HBig[1]
            HBig2d[1,3] = lens.HBig[1]
            HBig2d[1,4] = lens.HBig[1]
            HBig2d[2,0] = lens.HBig[2]
            HBig2d[2,1] = lens.HBig[2]
            HBig2d[2,2] = lens.HBig[2]
            HBig2d[2,3] = lens.HBig[2]
            HBig2d[2,4] = lens.HBig[2]
            GBig2d = np.zeros((3,5), dtype=np.complex)
            GBig2d[0,0] = lens.GBig[0]
            GBig2d[0,1] = lens.GBig[0]
            GBig2d[0,2] = lens.GBig[0]
            GBig2d[0,3] = lens.GBig[0]
            GBig2d[0,4] = lens.GBig[0]
            GBig2d[1,0] = lens.GBig[1]
            GBig2d[1,1] = lens.GBig[1]
            GBig2d[1,2] = lens.GBig[1]
            GBig2d[1,3] = lens.GBig[1]
            GBig2d[1,4] = lens.GBig[1]
            GBig2d[2,0] = lens.GBig[2]
            GBig2d[2,1] = lens.GBig[2]
            GBig2d[2,2] = lens.GBig[2]
            GBig2d[2,3] = lens.GBig[2]
            GBig2d[2,4] = lens.GBig[2]
            HBigP2d = np.zeros((2,5), dtype=np.complex)
            HBigP2d[0,0] = HBigP[0]
            HBigP2d[0,1] = HBigP[0]
            HBigP2d[0,2] = HBigP[0]
            HBigP2d[0,3] = HBigP[0]
            HBigP2d[0,4] = HBigP[0]
            HBigP2d[1,0] = HBigP[1]
            HBigP2d[1,1] = HBigP[1]
            HBigP2d[1,2] = HBigP[1]
            HBigP2d[1,3] = HBigP[1]
            HBigP2d[1,4] = HBigP[1]
            GBigP2d = np.zeros((2,5), dtype=np.complex)
            GBigP2d[0,0] = GBigP[0]
            GBigP2d[0,1] = GBigP[0]
            GBigP2d[0,2] = GBigP[0]
            GBigP2d[0,3] = GBigP[0]
            GBigP2d[0,4] = GBigP[0]
            GBigP2d[1,0] = GBigP[1]
            GBigP2d[1,1] = GBigP[1]
            GBigP2d[1,2] = GBigP[1]
            GBigP2d[1,3] = GBigP[1]
            GBigP2d[1,4] = GBigP[1]
#        HBig2d = np.reshape(np.repeat(lens.HBig, nz), (3, nz))
#        GBig2d = np.reshape(np.repeat(lens.GBig, nz), (3, nz))
        HBigScalar = np.sum(HBig2d * zall, axis=0)
        GBigScalar = np.sum(GBig2d * zall, axis=0)
        zallP = zall[0:2,:]
        #HBigP2d = np.reshape(np.repeat(HBigP, nz), (2, nz))
        #GBigP2d = np.reshape(np.repeat(GBigP, nz), (2, nz))
        HBigPScalar = np.sum(HBigP2d * zallP, axis=0)
        GBigPScalar = np.sum(GBigP2d * zallP, axis=0)
        #HBigPP = np.array([2.0])
        #zallPP = zall[0,:]
        #HBigPP2d = np.reshape(np.repeat(HBigPP, nz), (1, nz))
        #HBigPPScalar = np.sum(HBigPP2d * zallPP, axis=0)
        zkappas = (GBigPScalar * HBigScalar - GBigScalar * HBigPScalar) / \
                 (HBigScalar * HBigScalar)
        zkappaPs = 0.0
#        zkappaPs = (2.0 * GBigScalar * HBigPScalar * HBigPScalar - \
#                    HBigScalar * (2.0 * GBigPScalar * HBigPScalar + \
#                    GBigScalar * HBigPPScalar)) / \
#                    (HBigScalar * HBigScalar * HBigScalar)
        jacs = 1.0 - np.abs(zkappas) * np.abs(zkappas)
        paritys  = np.where(jacs >= 0.0, 1, -1)
    #print "=> elasped c: %s s" % t.secs
    
        return zlocs, jacs, paritys, zkappas, zkappaPs

def solve_binary_critcurve(theta, r1, r2, e1, e2):
    """The image plane critical curves for binary lens can be found
       by solving a 4th order polynomial.  The critical curve complex 
       positions in the image plane are parameterized by the polar angle, theta.
       This function solves for the critical points for a single theta.
       The polynomial coefficient expressions come from Miller M.(2013)
        http://researcharchive.vuw.ac.nz/xmlui/handle/10063/2775

       INPUTS:
         theta - [radians] polar angle to find critical points for
         r1 - [complex] - complex position of lens 1
         r2 - [complex] - complex position of lens 2
         e1 - [float] - mass of lens 1 relative to total mass
         e2 - [float] - mass of lens 2 relative to total mass
       OUTPUTS:
         zlocs - [complex] numpy array of critical curve points 
    """
    et = np.exp(theta*1j)
    GBig = np.array([-e1*r2 - e2*r1, 1.0, 0.0])
    HBig = np.array([r1*r2, -r1 - r2, 1.0])
    c0 = (GBig[1] * HBig[0] - GBig[0] * HBig[1]) - HBig[0] * HBig[0] * et
    c1 = 2.0 * (GBig[2] * HBig[0] - GBig[0] * HBig[2]) - \
            2.0 * HBig[0] * HBig[1] * et
    c2 = (GBig[2] * HBig[1] - GBig[1] * HBig[2]) - \
            (HBig[1] * HBig[1] + 2.0 * HBig[0] * HBig[2]) * et
    c3 = -2.0 * HBig[1] * HBig[2] * et
    c4 = -HBig[2] * HBig[2] * et
    # direct coeffs from Paul Chote
    #a0 = e1*r2*r2 + e2*r1*r1 - r1*r1*r2*r2*et
    #a1 = -2.0*(e1*r2+e2*r1) + 2.0*r1*r2*(r1+r2)*et
    #a2 = 1.0 - (r1*r1 + r2*r2 + 4.0*r1*r2)*et
    #a3 = 2.0*(r1+r2)*et
    #a4 = -et    
    zlocs = ply.polyroots((c0, c1, c2, c3, c4))
    z0 = np.full_like(zlocs, 1.0)
    z2 = zlocs * zlocs
    zall = np.reshape(np.append(z0, np.append(zlocs, z2)),(3,4))
    HBig2d = np.reshape(np.repeat(HBig, 4), (3, 4))
    GBig2d = np.reshape(np.repeat(GBig, 4), (3, 4))
    HBigScalar = np.sum(HBig2d * zall, axis=0)
    GBigScalar = np.sum(GBig2d * zall, axis=0)
    GBigP = np.array([1.0, 0.0])
    HBigP = np.array([-r1 -r2, 2.0])
    zallP = zall[0:2,:]
    HBigP2d = np.reshape(np.repeat(HBigP, 4), (2, 4))
    GBigP2d = np.reshape(np.repeat(GBigP, 4), (2, 4))
    HBigPScalar = np.sum(HBigP2d * zallP, axis=0)
    GBigPScalar = np.sum(GBigP2d * zallP, axis=0)
    zdiff = np.abs((GBigPScalar * HBigScalar - GBigScalar * HBigPScalar) / \
             (HBigScalar * HBigScalar) - et)
    # return the image locations that are consistent with lens equation
    idx = np.where(zdiff < 1.0e-6)[0]
    #print zwdiff
    return zlocs[idx]

def trace_two_images(dblzlocs):
    """ Given a flat list of two image locations 
        identify where the image locations have flipped order
        Flipping order identification is accomplished by finding
        jumps in neighboring image locations
        INPUTS:
          dblzlocs - [complex] (2*n, 1) numpy array 
                     If image i from source vertex j is defined by z_ij
                     dblzlocs is a flat/1-d array of image positions
                     [z11, z21, z12, z22, z13, z23, .., z1n, z2n]
                     where there are n source vertices
        OUTPUTS:
          imageLocs - [complex] (2, n) numpy array where 
                     the image positions that belong in the same curve
                     are separated out in separate rows of the 2d array
    """
    # get every other image location
    topzlocs = dblzlocs[0::2]
    srcn=len(topzlocs)
    imageLocs = np.zeros((2, len(topzlocs)), dtype=np.complex)
    topzlocs2d = np.reshape(np.repeat(topzlocs, 2), (srcn, 2))
    dblzlocs2d = np.reshape(dblzlocs, (srcn, 2))
    topzlocs2d = np.roll(topzlocs2d, 1, axis=0)
    dist2d = np.abs(topzlocs2d - dblzlocs2d)
    idxmap = np.where(np.argmin(dist2d, axis=1) == 1)[0]
    # No ordering changes detected
    if len(idxmap) == 0:
        imageLocs[0,:] = topzlocs
        imageLocs[1,:] = dblzlocs[1::2]
    # Handle cases where ordering has flipped
    else:
        istrt = np.append(np.array([0]), idxmap)
        iend = np.append(idxmap, np.array([srcn]))
        rowone = 0
        rowtwo = 1
        for k in range(0,len(istrt)):
            imageLocs[0, istrt[k]:iend[k]] = dblzlocs2d[istrt[k]:iend[k], rowone]
            imageLocs[1, istrt[k]:iend[k]] = dblzlocs2d[istrt[k]:iend[k], rowtwo]
            if rowone == 0:
                rowone = 1
                rowtwo = 0
            else:
                rowone = 0
                rowtwo = 1
    return imageLocs

def trace_three_images(trpzlocs):
    """ Given a flat list of three image locations 
        identify where the image locations have flipped order
        Flipping order identification is accomplished by finding
        jumps in neighboring image locations.
        Uses a divide and conquer approach to find at least
        one of the three images that hasnt flipped order, then
        it calls trace_two_images to handle to remaining two images
         that have flipped order.
        So far I have not encountered a scenario where more than one
        pair of images has flipped.  The code at the end almost certainly
        doesnt work correctly and will need to debugged if a scenario arises
        that it is needed.   So far I have not found it to occur, thus I have
        not been able to confirm the code there.
        INPUTS:
          trpzlocs - [complex] (3*n, 1) numpy array 
                     If image i from source vertex j is defined by z_ij
                     trpzlocs is a flat/1-d array of image positions
                     [z11, z21, z31, z12, z22, z32, z13, z23, z33, .., z1n, z2n, z3n]
                     where there are n source vertices
        OUTPUTS:
          imageLocs - [complex] (3, n) numpy array where 
                     the image positions that belong in the same curve
                     are separated out in separate rows of the 2d array
    """
    # get every third image location
    topzlocs = trpzlocs[0::3]
    midzlocs = trpzlocs[1::3]
    botzlocs = trpzlocs[2::3]
    imgn = len(topzlocs)
    imageLocs = np.zeros((3, imgn), dtype=np.complex)
    topzlocs2d = np.reshape(np.repeat(topzlocs, 2), (imgn, 2))
    midzlocs2d = np.reshape(np.repeat(midzlocs, 2), (imgn, 2))
    botzlocs2d = np.reshape(np.repeat(botzlocs, 2), (imgn, 2))
    topzlocsoff2d = np.roll(topzlocs2d, 1, axis=0)
    midzlocsoff2d = np.roll(midzlocs2d, 1, axis=0)
    botzlocsoff2d = np.roll(botzlocs2d, 1, axis=0)

    dist_12_2d = np.abs(topzlocsoff2d - np.transpose(np.vstack((topzlocs, midzlocs))))
    dist_13_2d = np.abs(topzlocsoff2d - np.transpose(np.vstack((topzlocs, botzlocs))))
    dist_21_2d = np.abs(midzlocsoff2d - np.transpose(np.vstack((midzlocs, topzlocs))))
    #dist_22_2d = np.abs(midzlocsoff2d - midzlocs2d)
    dist_23_2d = np.abs(midzlocsoff2d - np.transpose(np.vstack((midzlocs, botzlocs))))
    dist_31_2d = np.abs(botzlocsoff2d - np.transpose(np.vstack((botzlocs, topzlocs))))
    dist_32_2d = np.abs(botzlocsoff2d - np.transpose(np.vstack((botzlocs, midzlocs))))
    #dist_33_2d = np.abs(botzlocsoff2d - botzlocs2d)
    
    idxmap = {}
    idxmap[0] = np.where(np.argmin(dist_12_2d, axis=1) == 1)[0]
    idxmap[1] = np.where(np.argmin(dist_13_2d, axis=1) == 1)[0]
    idxmap[2] = np.where(np.argmin(dist_21_2d, axis=1) == 1)[0]
    idxmap[3] = np.where(np.argmin(dist_23_2d, axis=1) == 1)[0]
    idxmap[4] = np.where(np.argmin(dist_31_2d, axis=1) == 1)[0]
    idxmap[5] = np.where(np.argmin(dist_32_2d, axis=1) == 1)[0]
    
    foundThem = False
    # Check if image one never changed
    if len(idxmap[0]) == 0 and len(idxmap[1]) == 0:
        imageLocs[0,:] = topzlocs
        dblzlocs = np.zeros((2*imgn,), dtype=np.complex)
        dblzlocs[0::2] = midzlocs
        dblzlocs[1::2] = botzlocs
        tempImageLocs = trace_two_images(dblzlocs)
        imageLocs[1,:] = tempImageLocs[0,:]
        imageLocs[2,:] = tempImageLocs[1,:]
        foundThem = True
    # Check if image two never changed
    if len(idxmap[2]) == 0 and len(idxmap[3]) == 0 and (not foundThem):
        imageLocs[0,:] = midzlocs
        dblzlocs = np.zeros((2*imgn,), dtype=np.complex)
        dblzlocs[0::2] = topzlocs
        dblzlocs[1::2] = botzlocs
        tempImageLocs = trace_two_images(dblzlocs)
        imageLocs[1,:] = tempImageLocs[0,:]
        imageLocs[2,:] = tempImageLocs[1,:]
        foundThem = True
    # Check if image three never changed
    if len(idxmap[4]) == 0 and len(idxmap[5]) == 0 and (not foundThem):
        imageLocs[0,:] = botzlocs
        dblzlocs = np.zeros((2*imgn,), dtype=np.complex)
        dblzlocs[0::2] = topzlocs
        dblzlocs[1::2] = midzlocs
        tempImageLocs = trace_two_images(dblzlocs)
        imageLocs[1,:] = tempImageLocs[0,:]
        imageLocs[2,:] = tempImageLocs[1,:]
        foundThem = True
    # Handle case of at least two images changed order
    # This needs more testing and probably slow
    if (not foundThem):
        print "This is preliminary code not thoroughly tested"
        print "Please test it if you actually get here!"
        blcks = 1
        curloc = 0
        currow = 0
        rows = [0, 1, 2]
        nexrow2d = np.array([1, 2], [0, 2], [0, 1])
        zlocs2d = np.vstack((topzlocs, midzlocs, botzlocs))
        tmpdblzlocs = np.zeros((2,imgn), dtype=np.complex)
        while curloc < imgn:
            chk1 = currow * 2
            idx1 = np.where(idxmap[chk1] > curloc)[0]
            idx2 = np.where(idxmap[chk1+1] > curloc)[0]
            tmp = [idxmap[chk1][idx1], idxmap[chk1+1][idx2]]
            blckend = min(tmp)
            nexrow = np.argmin(tmp)
            imageLocs[0, curloc:blckend+1] = zlocs2d[currow, curloc:blckend+1]
            tmpdblzlocs[0, curloc:blckend+1] =  \
                    zlocs2d[np.delete(rows, currow)[0], curloc:blckend+1]
            tmpdblzlocs[1, curloc:blckend+1] = \
                    zlocs2d[np.delete[rows, currow][1], curloc:blckend+1]
            blcks = blcks + 1
            curloc = blckend + 1
            currow = nexrow2d[currow, nexrow]
        # Fill in last block
        imageLocs[0, curloc:] = zlocs2d[currow, curloc:]
        tmpdblzlocs[0, curloc:] = \
                    zlocs2d[np.delete(rows, currow)[0], curloc:]
        tmpdblzlocs[1, curloc:] = \
                    zlocs2d[np.delete[rows, currow][1], curloc:]
        dblzlocs = np.zeros((2*imgn,), dtype=np.complex)
        dblzlocs[0::2] = tmpdblzlocs[0,:]
        dblzlocs[1::2] = tmpdblzlocs[1,:]
        tempImageLocs = trace_two_images(dblzlocs)
        imageLocs[1,:] = tempImageLocs[0,:]
        imageLocs[2,:] = tempImageLocs[1,:]                                
    return imageLocs

def handle_three_images(nSrc, zlocs, pars):
    """ Given a flat list of image locations, zlocs,
         and parities, pars, for nSrc source vertices return an array
         where the image border tracings match.
         Uses a divide and conquer approach.  With 3 images
         two will have the same parity.  Identify image with unique parity
         and handle that trivial case.  Send the other 2 image arrays to
           trace_two_images to handle.
        INPUTS:
          nSrc - [int] number of source vertices
          zlocs - [complex] (3*nSrc, 1) numpy array 
                     If image i from source vertex j is defined by z_ij
                     trpzlocs is a flat/1-d array of image positions
                     [z11, z21, z31, z12, z22, z32, z13, z23, z33, .., z1n, z2n, z3n]
                     where there are n source vertices
          pars - [1, -1] same length and shape as zocs given parity of image
                
        OUTPUTS:
          imageLocs - [complex] (3, n) numpy array where 
                     the image positions that belong in the same curve
                     are separated out in separate rows of the 2d array
          imagePars - [1, -1] same shape as imageLocs giving parities of images
    """
    imageLocs = np.zeros([3, nSrc], dtype=np.complex)
    imagePars = np.zeros([3, nSrc])
    # Find parity with only one image
    if (pars[0] + pars[1] + pars[2]) == -1:
        singleparity = 1
        doubleparity = -1
    else:
        singleparity = -1
        doubleparity = 1
    # trivial image locations with unique parity copy it over
    idx = np.where(pars == singleparity)[0]
    imageLocs[0,:] = zlocs[idx]
    imagePars[0,:] = np.full((nSrc,), singleparity)
    # get image locations with same parity together and send to trace_two_images 
    idx = np.where(pars == doubleparity)[0]
    dblzlocs = zlocs[idx]
    tempImageLocs = trace_two_images(dblzlocs)
    imageLocs[1,:] = tempImageLocs[0,:]
    imagePars[1,:] = np.full((nSrc,), doubleparity)
    imageLocs[2,:] = tempImageLocs[1,:]            
    imagePars[2,:] = np.full((nSrc,), doubleparity)

    return imageLocs, imagePars

def handle_five_images(nSrc, zlocs, pars):
    """ Given a flat list of image locations, zlocs,
         and parities, pars, for nSrc source vertices return an array
         where the image border tracings match.
         Uses a divide and conquer approach.  With 5 images
         two will have the same parity and there are another 3 with other parity.
         Identify images with doubled parity and send to trace_two_images to handle.
         Identify images with tripled parity and send to trace_three_images to handle
        INPUTS:
          nSrc - [int] number of source vertices
          zlocs - [complex] (5*nSrc, 1) numpy array 
                     If image i from source vertex j is defined by z_ij
                     trpzlocs is a flat/1-d array of image positions
                     [z11, z21, z31, z41, z51, .., z1n, z2n, z3n, z4n, z5n]
                     where there are n source vertices
          pars - [1, -1] same length and shape as zocs given parity of image
                
        OUTPUTS:
          imageLocs - [complex] (5, n) numpy array where 
                     the image positions that belong in the same curve
                     are separated out in separate rows of the 2d array
          imagePars - [1, -1] same shape as imageLocs giving parities of images
    """
    imageLocs = np.zeros([5, nSrc], dtype=np.complex)
    imagePars = np.zeros([5, nSrc])
    # Find parity with two images
    if (pars[0] + pars[1] + pars[2] + pars[3] + pars[4]) == -1:
        doubleparity = 1
        tripleparity = -1
    else:
        doubleparity = -1
        tripleparity = 1
    # Send doubled parity to trace_two_images
    idx = np.where(pars == doubleparity)[0]
    dblzlocs = zlocs[idx]
    tempImageLocs = trace_two_images(dblzlocs)
    imageLocs[0,:] = tempImageLocs[0,:]
    imagePars[0,:] = np.full((nSrc,), doubleparity)
    imageLocs[1,:] = tempImageLocs[1,:]
    imagePars[1,:] = np.full((nSrc,), doubleparity)
    # Send tripled parity to trace_three_images
    idx = np.where(pars == tripleparity)[0]
    trpzlocs = zlocs[idx]
    tempImageLocs = trace_three_images(trpzlocs)
    for k in range(0,3):
        imageLocs[k+2,:] = tempImageLocs[k,:]
        imagePars[k+2,:] = np.full((nSrc,), tripleparity)

    return imageLocs, imagePars
    
def finite_source_mag(src_list, debugLevel=0):
    """From the list of source vertex objects perform the line integral
       around the source images to calculate their area and magnification in 
       the uniform surface brightness finite source scenartio.
       Before calculating areas one has to connect the image curves in the
       scenario of the finite source being occulted/split up by a caustic curve.
       This function could really be called handle the nauseating bookkeeping of
       sources strattling a cautstic curve because 90% of this functions code
       is actually doing that rather than calculating the magnitude.
       For the image curve matching in this function to work, 
       the source outline has to have sufficient resolution of source vertex objects.
       So errors in properly tracing the images result in this function, but 
       it may actually be sympotmatic of not having enough source vertex resolution,
        which is handled by fill_image_gaps().  So when debugging increasing
        source_vertex resolution may solve the issue you have with this function.
        In other words, the inputs have to be sane for this function to work
        properly, and this function doesnt check that the inputs are sane ahead
        of time.
       INPUTS:
        src_list - list of source_plane_vertex objects
                    It is assumed that the source vertices have been sorted
                    by their parametric polar angle.
       OUTPUTS:
         curmag - [float] magnification of image area relative to the unlensed
                     source area
    """
    # The information for each source vertex object is unrolled into the 1-d arrays
    #   that are actually used in calculation 
    aux_arrays = get_auxillary_arrays(src_list)
    zlocs = aux_arrays[0]
    stpns = np.int64(aux_arrays[1])
    jacs = aux_arrays[3]
    pars = np.int64(aux_arrays[4])
    zkaps = aux_arrays[5]
    zkapPs = aux_arrays[6]
    # Identify how many images there are for every source vertex
    stpidx = np.where(np.diff(stpns) == 1)[0]
    stpidx = np.append(np.array([-1]), stpidx)
    stpidx = np.append(stpidx, np.array([len(stpns)-1]))
    # This vector tells how many images there are at every source vertex
    stplen = np.diff(stpidx)
    maxNImages = max(stplen)
    nSrc = max(stpns)+1
    imageLocs = np.full([maxNImages, nSrc], -99.0, dtype=np.complex)
    imagePars = np.zeros([maxNImages, nSrc])
    foundThem = False
    # Handle no caustic crossing case.  All sources have 3 and only 3 images
    if maxNImages == 3:
        #print "3 Full"
        imageLocs, imagePars = handle_three_images(nSrc, zlocs, pars)
        foundThem = True
        # put image tracings into the loop dictionary format
        loopN = 3
        loopLocs = {}
        for kk in range(3):
            loopLocs[kk] = imageLocs[kk,:]
    # Handle fully in caustic case.  All sources have 5 and only 5 images
    if np.all(stplen == 5):
        #print "5 Full"
        imageLocs, imagePars = handle_five_images(nSrc, zlocs, pars)
        foundThem = True
        # put image tracings into the loop dictionary format
        loopN = 5
        loopLocs = {}
        for kk in range(5):
            loopLocs[kk] = imageLocs[kk,:]
    if (not foundThem):
        # Here is the neuseating bookkeeping scenario of source
        #  strattling a caustic and multiple images are created and destroyed
        #  along the source outline! Ick.
        # In summary it works by blocking up streams of 3 and 5 images
        #  Handles the 3 and 5 image blocks independently
        #  Then heuristically matches the image curves across the 3->5 or 5->3
        # image boundaries through a close distance and angle metric
        #print "Mix"
        # First rearrange things to start on a outside caustic 3 image source veretex
        idx = np.where(stplen == 3)[0]
        if (idx[0] != 0):
            delayidx = np.where(stpns == idx[0])[0]
            totndelay = len(zlocs) - delayidx[-1] + 2
            zlocs = np.roll(zlocs, totndelay)
            stpns = np.roll(stpns, totndelay)
            basestp = stpns[0]
            stpns = np.where(stpns >= basestp, stpns-basestp, \
                            stpns+nSrc-basestp)
            jacs = np.roll(jacs, totndelay)
            pars = np.roll(pars, totndelay)
            zkaps = np.roll(zkaps, totndelay)
            zkapPs = np.roll(zkapPs, totndelay)
            stpidx = np.where(np.diff(stpns) == 1)[0]
            stpidx = np.append(np.array([-1]), stpidx)
            stpidx = np.append(stpidx, np.array([len(stpns)-1]))
            # This vector tells how many images there are at every source vertex
            stplen = np.diff(stpidx)
        else:
            if (idx[-1] == len(stplen)-1):
                idx2 = np.where(stplen == 5)[0]      
                delayidx = np.where(stpns == idx2[-1])[0]
                totndelay = len(zlocs) - delayidx[-1] -1
                zlocs = np.roll(zlocs, totndelay)
                stpns = np.roll(stpns, totndelay)
                basestp = stpns[0]
                stpns = np.where(stpns >= basestp, stpns-basestp, \
                            stpns+nSrc-basestp)
                jacs = np.roll(jacs, totndelay)
                pars = np.roll(pars, totndelay)
                zkaps = np.roll(zkaps, totndelay)
                zkapPs = np.roll(zkapPs, totndelay)
                stpidx = np.where(np.diff(stpns) == 1)[0]
                stpidx = np.append(np.array([-1]), stpidx)
                stpidx = np.append(stpidx, np.array([len(stpns)-1]))
                # This vector tells how many images there are at every source vertex
                stplen = np.diff(stpidx)
                
        # Identify blocks of 3 or 5 images where stplen changes number                       
        idxBlck = np.where(np.abs(np.diff(stplen)) > 0)[0]
        strtidx = 0
        strtidx2 = 0
        lstLocs = np.full((5,), -99.0, dtype=np.complex)
        curimg = 3
        allstrts = np.array([])
        allends = np.array([])
        # handle each block and try to match the previous block
        # if you step through the code we always start with 3 images
        #  then alternate between 5, 3, 5, 3... until the last block
        for k in range(0,len(idxBlck)):
            idx = np.where(stpns == idxBlck[k])[0]
            endidx = idx[-1]
            tz = zlocs[strtidx:endidx+1]
            ps = pars[strtidx:endidx+1]
            curn = idxBlck[k] - strtidx2 + 1
            if curimg == 3:
                endidx2= strtidx2 + curn
                tmpLocs, tmpPars = handle_three_images(curn, tz, ps)
                rows = range(0,3)
                # On initial block just fill in values
                if k == 0:
                    for kk in range(0,3):
                        imageLocs[kk,strtidx2:endidx2] = tmpLocs[kk,:]
                        imagePars[kk,strtidx2:endidx2] = tmpPars[kk,:]
                        lstLocs[kk] = tmpLocs[kk,-1]
                # Otherwise match blocks to last positions
                else:
                    frstLocs = tmpLocs[:,0]
                    for kk in range(0,3):
                        useidx = np.argmin(np.abs(lstLocs - frstLocs[kk]))# + \
                                            #np.abs(lstPars[kk] - frstPars))
                        imageLocs[useidx, strtidx2:endidx2] = tmpLocs[kk, :]
                        imagePars[useidx, strtidx2:endidx2] = tmpPars[kk, :]
                        lstLocs[useidx] = tmpLocs[kk, -1]
                curimg = 5
            else:
                endidx2 = strtidx2 + curn
                tmpLocs, tmpPars = handle_five_images(curn, tz, ps)
                rows = range(0,5)
                # On first 5 image block match first three then fill in values
                #  for other 2 images
                if k == 1:
                    frstLocs = tmpLocs[:,0]
                    for kk in range(0,3):
                        useidx = np.argmin(np.abs(lstLocs[kk] - frstLocs))# + \
 #                                           np.abs(lstPars[kk] - frstPars))
                        imageLocs[kk,strtidx2:endidx2] = tmpLocs[useidx,:]
                        imagePars[kk,strtidx2:endidx2] = tmpPars[useidx,:]
                        lstLocs[kk] = tmpLocs[useidx, -1]
                        idx = np.where(rows == useidx)[0]
                        rows = np.delete(rows, idx)
                    for kk in range(3,5):
                        imageLocs[kk,strtidx2:endidx2] = tmpLocs[rows[kk-3], :]
                        imagePars[kk,strtidx2:endidx2] = tmpPars[rows[kk-3], :]
                        lstLocs[kk] = tmpLocs[rows[kk-3], -1]
                else:
                    frstLocs = tmpLocs[:,0]
                    lastIdx = np.where(imageLocs[:, strtidx2-1] > -10.0)[0]
                    finIdx = np.where(imageLocs[:, strtidx2-1] < -10.0)[0]
                   # first match the 3 rows that continue image
                    for kk in range(0,3):
                        useidx = np.argmin(np.abs(lstLocs[lastIdx[kk]] - frstLocs))# + \
                                           # np.abs(lstPars[kk] - frstPars))
                        imageLocs[lastIdx[kk],strtidx2:endidx2] = tmpLocs[useidx,:]
                        imagePars[lastIdx[kk],strtidx2:endidx2] = tmpPars[useidx,:]
                        lstLocs[lastIdx[kk]] = tmpLocs[useidx, -1]
                        idx = np.where(rows == useidx)[0]
                        rows = np.delete(rows, idx)
                    # For the last twwo rows put them in previously blank
                        # rows.  These extra curves will get put in the 
                        #correct place later
                    for kk in range(0,2):
                        useidx = rows[0]
                        imageLocs[finIdx[kk],strtidx2:endidx2] = tmpLocs[useidx,:]
                        imagePars[finIdx[kk],strtidx2:endidx2] = tmpPars[useidx,:]
                        lstLocs[finIdx[kk]] = tmpLocs[useidx, -1]
                        idx = np.where(rows == useidx)[0]
                        rows = np.delete(rows, idx)
                curimg = 3
            allstrts = np.append(allstrts, strtidx2)
            allends = np.append(allends, endidx2)
            strtidx = endidx+1
            strtidx2 = endidx2
        # Do last block
        tz = zlocs[strtidx:]
        ps = pars[strtidx:]
        curn = nSrc - strtidx2
        if curimg == 3:
            print "Warning untested code!"
            print "One shouldnt be here I dont think"
            endidx2 = strtidx2 + curn
            tmpLocs, tmpPars = handle_three_images(curn, tz, ps)
            rows = range(0,3)
            frstLocs = tmpLocs[:,0]
            for kk in range(0,3):
                useidx = np.argmin(np.abs(lstLocs[kk] - frstLocs))# + \
                                #np.abs(lstPars[kk] - frstPars))
                imageLocs[kk, strtidx2:endidx2] = tmpLocs[useidx, :]
                imagePars[kk, strtidx2:endidx2] = tmpPars[useidx, :]
                #idx = np.where(rows == useidx)[0]
                #rows = np.delete(rows, idx)
            #if (len(rows) != 0):
                #print "Warning Not all rows used!! 3 Last"
        else:
            endidx2 = strtidx2 + curn
            tmpLocs, tmpPars = handle_five_images(curn, tz, ps)
            rows = range(0,5)
            frstLocs = tmpLocs[:,0]
            lastIdx = np.where(imageLocs[:, strtidx2-1] > -10.0)[0]
            finIdx = np.where(imageLocs[:, strtidx2-1] < -10.0)[0]

            for kk in range(0,3):
                useidx = np.argmin(np.abs(lstLocs[lastIdx[kk]] - frstLocs))# + \
                                    #np.abs(lstPars[kk] - frstPars))
                imageLocs[lastIdx[kk],strtidx2:endidx2] = tmpLocs[useidx,:]
                imagePars[lastIdx[kk],strtidx2:endidx2] = tmpPars[useidx,:]
                idx = np.where(rows == useidx)[0]                
                rows = np.delete(rows, idx)
            for kk in range(2):
                useidx = rows[0]
                imageLocs[finIdx[kk],strtidx2:endidx2] = tmpLocs[useidx,:]
                imagePars[finIdx[kk],strtidx2:endidx2] = tmpPars[useidx,:]
                idx = np.where(rows == useidx)[0]
                rows = np.delete(rows, idx)
        allstrts = np.append(allstrts, strtidx2)
        allends = np.append(allends, endidx2)
        
        # Now we have a bunch of image pieces
        # Some rows are complete continuous images
        #  Other rows are pieces of images with gaps at 3/5 image boundaries
        #  and the ordering is flipped because of parity differences
        # but the pieces belong to the same 'loop'/outline of image
        rows = range(0, maxNImages)
        # First check which rows are complete independent images
        frstLocs = imageLocs[:,0]
        lstLocs = imageLocs[:,-1]
        indrows = np.where(np.abs(frstLocs - lstLocs) < 5.0e-3, 1, 0)
        for kk in range(5):
            # Check for row with gaps in it
            idx = np.where(imageLocs[kk,:] < -10.0)[0]
            if (len(idx) > 0): 
                indrows[kk] = 0
            # Check for row with jumps in it
            idx = np.where(np.abs(np.diff(imageLocs[kk,:])) > 0.4)[0]
            if (len(idx) > 0):
                indrows[kk] = 0
        idx = np.where(indrows == 1)[0]
        #print indrows
        # put these complete independent images into loops
        loopN = 0
        loopLocs = {}
        if len(idx) > 0:
            for kk in range(len(idx)):
                loopN = loopN + 1
                loopLocs[kk] = imageLocs[idx[kk],:]
                rows = np.delete(rows, np.where(rows == idx[kk])[0])
        # What remains are partial images
        # Break remaining rows into pieces
        crvN = 0
        crvLocs = {}
        crvPars = np.array([])        
        for kk in range(len(rows)):
            idxBlck = np.where(np.abs(np.diff(imageLocs[rows[kk],:])) > 0.4)[0]
            if (len(idxBlck) == 0):
                crvLocs[crvN] = imageLocs[rows[kk],:]
                crvPars = np.append(crvPars, imagePars[rows[kk],0])
                crvN = crvN + 1
                
            else:
                strtidx = 0
                for jj in range(len(idxBlck)):
                    endidx = idxBlck[jj]
                    if imageLocs[rows[kk],strtidx] > -10.0 and endidx != nSrc-2:
                        crvLocs[crvN] = imageLocs[rows[kk],strtidx:endidx+1]
                        crvPars = np.append(crvPars, imagePars[rows[kk],strtidx])
                        crvN = crvN + 1
                    if endidx != nSrc-2:
                        strtidx = endidx + 1
                # do last block
                if imageLocs[rows[kk],strtidx] >  -10.0:
                    crvLocs[crvN] = imageLocs[rows[kk],strtidx:nSrc-1]
                    crvPars = np.append(crvPars, imagePars[rows[kk],strtidx])
                    crvN = crvN + 1
        # We have all the pieces separated out take the first one in the
        # list and try to find the best match until a loop is 'closed'
        # by closed the end of the next piece is close to the beginning
        # repeat process until the pieces dictionary is empty
        # This is where pieces are matched by distance between end points
        #  and also angle/slope of end points agree.  The angle constraint
        #  was needed to handle very thin images stretched along einstein ring
        #  parity also has to be checked for each image piece in order to 
        #  match to the beginning or end of the piece.
        crvs = np.arange(0, len(crvLocs))
        #for kk in range(len(crvs)):
        #    print( [len(crvLocs[kk]), i])
        while len(crvs) > 0:
            lpstrtidx = crvs[0]
            tmp = crvLocs[lpstrtidx]
            lppar = crvPars[lpstrtidx]
            strtLoc = tmp[0]
            endLoc = tmp[-1]
            endAng = np.angle(tmp[-1] - tmp[-2])
            # normalize angle differences such that an angle shift
            #  of pi/2 is equivalent to the typical spacing
            #  of vertices
            angNorm = np.mean(np.abs(np.diff(tmp))) * 2.0 / (np.pi / 2.0)
            crvs = np.delete(crvs, np.where(crvs == lpstrtidx)[0])
            atStrt = False
            while not atStrt and len(crvs) > 0:
                compLocs = np.array([])
                compAngs = np.array([])
                for kk in range(len(crvs)):
                    curidx = crvs[kk]
                    ax1 = 1
                    ax2 = -2
                    if (len(crvLocs[curidx]) == 1):
                        ax1 = 0
                        ax2 = 0
                        
                    if crvPars[curidx] == lppar:
                        compLocs = np.append(compLocs, crvLocs[curidx][0])
                        compAngs = np.append(compAngs, \
                            np.angle(crvLocs[curidx][ax1] - crvLocs[curidx][0]))
                    else:
                        compLocs = np.append(compLocs, crvLocs[curidx][-1])
                        compAngs = np.append(compAngs, \
                            np.angle(crvLocs[curidx][ax2] - crvLocs[curidx][-1]))
                useidx = crvs[np.argmin(np.abs(endLoc - compLocs) + \
                                        np.abs(endAng - compAngs) * angNorm)]
                #print([lpstrtidx, useidx])
                # add curve to loop
                if crvPars[useidx] == lppar:
                    tmp = np.append(tmp, crvLocs[useidx])
                else:
                    tmp = np.append(tmp, np.flipud(crvLocs[useidx]))
                # update end location
                endLoc = tmp[-1]
                endAng = np.angle(tmp[-1] - tmp[-2])
                angNorm = np.mean(np.abs(np.diff(tmp))) * 2.0 / (np.pi / 2.0)
                # delete this curve from availability
                crvs = np.delete(crvs, np.where(crvs == useidx)[0])
                # Check if at start of loop
                if np.abs(endLoc - strtLoc) < min(np.abs(endLoc - compLocs)):
                        atStrt = True
            # Finished Loop add it to list
            loopLocs[loopN] = tmp
            loopN = loopN + 1

    if debugLevel > 1:
        # Show each individual image plane loop
        colorWheel = {0:'+r', 1:'+b', 2:'+k', 3:'+c', 4:'+m'}
        for kk in range(loopN):
            fig, ax = plt.subplots()
            ax.plot(np.real(loopLocs[kk]), np.imag(loopLocs[kk]), colorWheel[kk])
            if debugLevel == 2:
                plt.pause(0.5)
            if debugLevel == 3:
                print "Click anywhere in window to proceed"
                plt.waitforbuttonpress()
            if debugLevel == 4:
                fig.canvas.mpl_connect('key_press_event', on_key_event)    
                # Disable default Matplotlib shortcut keys:
                keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
                for key in keymaps:
                    plt.rcParams[key] = ''
                ax.set_title('Interactive Window, select a vowel on the keyboard to close window and proceed:')
                plt.show()
            plt.close()
            print "Debugging in finite_source_mag"
            
    # Whew!  Actually calculate the area of the source image loop
    indmags = np.zeros((loopN,))
    for k in range(loopN):
        x = np.real(loopLocs[k])
        y = np.imag(loopLocs[k])
        delx = np.roll(x, 1)
        dely = np.roll(y, 1)
        indmags[k] = np.abs(0.5 * np.sum(y*delx - x*dely))
    #print(indmags)
    curmag = np.sum(indmags)
    # Handle ring case
    if max(indmags) > 0.8:
        tmp = np.sort(indmags)
        tmp[-2] = -tmp[-2]
        curmag = np.sum(tmp)
        
    return curmag

def get_auxillary_arrays(src_list):
    """Go through source_plan_vertex objects and place image plane
       metrics into individual numpy arrays
    """
    # preallocate the maximum size possible
    maxN = len(src_list) * 5
    allzlocs = np.zeros((maxN,), dtype=np.complex)
    allstpns = np.zeros((maxN,), dtype=np.int64)
    allnimages = np.zeros((maxN,), dtype=np.int)
    alljacs = np.zeros((maxN,), dtype=np.complex)
    allparitys = np.zeros((maxN,), dtype=np.int)
    allkappas = np.zeros((maxN,), dtype=np.complex)
    allkappaDs = np.zeros((maxN,), dtype=np.complex)
    allstrtch = np.zeros((maxN,), dtype=np.int64)
    strt = 0
    end = 0
    for i in range(len(src_list)):
        obj = src_list[i]
        nImg = obj.nImages
        end = strt + nImg
        allzlocs[strt:end] = obj.imageLocs
        allstpns[strt:end] = i
        allnimages[strt:end] = nImg
        alljacs[strt:end] = obj.jacobians
        allparitys[strt:end] = obj.imageParity
        allkappas[strt:end] = obj.kappas
        allkappaDs[strt:end] = obj.kappaDerivs
        allstrtch[strt:end] = obj.stretchFactor
        strt = end
    # Truncate results to actual size
    # The preallocation and truncation is way faster than using np.append
    idx = range(end)
    return allzlocs[idx], allstpns[idx], allnimages[idx], alljacs[idx], \
            allparitys[idx],  allkappas[idx], allkappaDs[idx], allstrtch[idx]
                          
def fill_image_gaps(src_list, lens, debugLevel=0):
    """This is the routine to find gaps in neighboring image
       vetices and determine how many source vertices
       need to be added to fill the gap.  If the code provides
       innaccurate results this routine may need some fiddling. Because
       the image tracing that happens in finite_source_mag() relies
       on having the source outline be well filled in the image plane
       by this routine.
       INPUT:
         src_list - List of source_plane_vertex objects
       OUTPUT:
         new_src_list - new list of source_planet_vertiex objects
                        augmented to fill in image plane gaps
    """
    # unroll the source vertex object data into 1-d arrays of image locations
    aux_arrays = get_auxillary_arrays(src_list)
    allzlocs = aux_arrays[0]
    allstpns = np.int64(aux_arrays[1])
    actcnt = max(allstpns) + 1 # Number of sources
    
    srcx = src_list[0].centerx
    srcy = src_list[0].centery
    srcsz = src_list[0].radius
    # Handle a very special case that needs large boosting in resolution
    # to work correct
    #  For large source sizes relative to central caustic size
    # extreme source stretching occurs when the central caustic is crossing
    #  the large source limb.  For practical purposes source sizes are too small
    #  for this scenario to be realistic, but for testing it occurs and handled here
    sclfactor = 1.0
    if (srcsz > 3*lens.centCaustSize):
        srcCaustSep = np.abs(srcx + srcy*1j - lens.centCaustLoc)
        if (np.abs(srcCaustSep - srcsz) < 3.0*lens.centCaustSize):
            sclfactor = 6.0
            if debugLevel > 0:
                print([i, -999])
    # Here starts the main routine to identify vertices with large gaps
    # between their neighbor vertex relative to the typical spacing between
    #  vertices
    doneFilling = False
    fillCount = 0
    while (not doneFilling) and (fillCount < 2):
        aux_arrays = get_auxillary_arrays(src_list)
        allzlocs = aux_arrays[0]
        allstpns = np.int64(aux_arrays[1])
        actcnt = max(allstpns) + 1 # Number of sources
        # Clip out last source as it is identical to first\
        idx = np.where(allstpns < actcnt-1)[0]
        # Break the complex image positions into separate x and y components
        tmpx = np.real(allzlocs)[idx]
        tmpy = np.imag(allzlocs)[idx]
        tmpn  = len(tmpx)
        # Get the distance between all pairs of image positions
        dist1d = dist.pdist(np.hstack((tmpx[:,None], \
                    tmpy[:,None])), 'euclidean')
        # Get the distances into NxN form
        dist2d = dist.squareform(dist1d)
        # Now we need to find the distance between each image vertex
        # and its neighbor that comes from the next higher source angle
        # clip the distance array to the next 5 image positions in source angle
        tmpidxx = np.arange(tmpn)
        tmpidyy = np.arange(tmpn)
        subDist2d = np.zeros((7, tmpn))
        for kk in range(7):
            tmpidyy = np.roll(tmpidyy, -1)
            subDist2d[kk,:] = dist2d[tmpidxx, tmpidyy]
        neighborDist = np.min(subDist2d, axis=0)
        # neighborDist contains for every image vertex the distance
        # between its next image vertex with higher source angle
        # Next is to find the typical disances amongst image vertices
        #  nearby.  For every image vertex get the closest distances
        srtdist2dArg = np.argsort(dist2d, kind='heapsort')
        #srtdist2d = np.sort(dist2d, kind='heapsort')
        srtdist2d = dist2d[srtdist2dArg, np.tile(np.arange(tmpn), (tmpn, 1)).T]
        closestDist = srtdist2d[:,1]
        # Get the typical distance between nearest neighbors from the 
        #  number of source vertices. 
        # Use slightly fewer than the total number of source vertices
        #  so we reduce the chance of contaminating from a different
        #   image location with different scales for the typical distances
        useN = np.fix(actcnt/1.5)
        srtdist2dArgClip = srtdist2dArg[:,1:useN+1]
        usedistClip = closestDist[srtdist2dArgClip]
        typdist = np.mean(usedistClip, axis=1)
        # This is where we add some protection 
        # Dont let typical distance to get too big
        # Also apply sclfactor if you want to really boost resoltuion
        # Dont let typical distance get too small with midist
        maxdist = 0.1 / sclfactor
        mindist = 1.0e-3
        idx = np.where(typdist > maxdist)[0]
        if (len(idx) > 0):
            #print([i, len(idx)])
            typdist[idx] = maxdist
        idx = np.where(typdist < mindist)[0]
        if (len(idx) > 0):
            typdist[idx] = mindist
        # stretch factor is how many more vertices we should add between neighbors
        # in order to maintain the typical image vertice distance
        stretchfac = np.int64(np.round(neighborDist / typdist  ))
        # put the stretch factors into the source vertex objects
        # and add new sources if there are any that need augmenting
        # check if this is the last round of filling
        if max(stretchfac) < 3:
            doneFilling = True
        idx = np.where(stretchfac > 1)[0]
        if len(idx) > 0:
            srcids = allstpns[idx]
            stretchfac = stretchfac[idx]
            for kk in range(len(srcids)):
                curidx = srcids[kk]
                curstrtch = src_list[curidx].stretchFactor
                if (curstrtch < stretchfac[kk]):
                    src_list[curidx].stretchFactor = stretchfac[kk]
            # Add new sources
            new_src_list = []
            # Get base source angle separation
            srcx = src_list[0].centerx
            srcy = src_list[0].centery
            srcsz = src_list[0].radius
            srcn = len(src_list)
            for kk in range(len(src_list)):
                obj=src_list[kk]
                nexobj = src_list[np.mod(kk+1, srcn)]
                if obj.stretchFactor > 1:
                    newAngles = np.linspace(obj.angle, nexobj.angle, \
                                obj.stretchFactor + 2)
                    # trim off first and last newAngles that match previous srcs
                    newAngles = newAngles[1:-1]
                    for zz in range(len(newAngles)):
                        newobj = source_plane_vertex()
                        newobj.set_source_position(srcx, srcy, srcsz, newAngles[zz])
                        zlocs, jacs, paritys, zkappas, zkappaPs = \
                            solve_binary_roots(newobj, lens)
                        newobj.nImages = len(zlocs)
                        newobj.imageLocs = zlocs
                        newobj.imageParity = paritys
                        newobj.jacobians = jacs
                        newobj.kappas = zkappas
                        newobj.kappaDerivs = zkappaPs
                        # set stretch Factor -1 to indicate this was an added source
                        newobj.stretchFactor = -1
                        new_src_list.append(newobj)
            # Add new sources to old ones and sort them by source angle
            for kk in range(len(new_src_list)):
                src_list.append(new_src_list[kk])
            src_list = sorted(src_list, key=attrgetter('angle'))

    
        # Now augment source vertices at caustic crossings
        # Add caustboost source vertices in boundary between
        #  3 and 5 images
        caustboost = 6
        aux_arrays = get_auxillary_arrays(src_list)
        allzlocs = aux_arrays[0]
        allnImag = np.int64(aux_arrays[2])
        allstpns = np.int64(aux_arrays[1])
        idximg = np.where(np.diff(allnImag) != 0)[0]
        if len(idximg) > 0:
            srcids = allstpns[idximg]
            # Add new sources
            new_src_list = []
            srcx = src_list[0].centerx
            srcy = src_list[0].centery
            srcsz = src_list[0].radius
            for kk in range(len(srcids)):
                curid = srcids[kk]
                obj=src_list[curid]
                deltaAngle = src_list[curid+1].angle - src_list[curid].angle
                newAngles = np.linspace(obj.angle, obj.angle+deltaAngle, \
                                caustboost + 2)
                # trim off first and last newAngles that match previous srcs
                newAngles = newAngles[1:-1]
                for zz in range(len(newAngles)):
                    newobj = source_plane_vertex()
                    newobj.set_source_position(srcx, srcy, srcsz, newAngles[zz])
                    zlocs, jacs, paritys, zkappas, zkappaPs = \
                            solve_binary_roots(newobj, lens)
                    newobj.nImages = len(zlocs)
                    newobj.imageLocs = zlocs
                    newobj.imageParity = paritys
                    newobj.jacobians = jacs
                    newobj.kappas = zkappas
                    newobj.kappaDerivs = zkappaPs
                    # set stretch Factor -2 to indicate this was an added source
                    newobj.stretchFactor = -2
                    new_src_list.append(newobj)
            # Add new sources to old ones and sort them by source angle
            for kk in range(len(new_src_list)):
                src_list.append(new_src_list[kk])
            src_list = sorted(src_list, key=attrgetter('angle'))
        fillCount = fillCount + 1

    if debugLevel > 1:
        aux_arrays = get_auxillary_arrays(src_list)
        allzlocs = aux_arrays[0]
        stretchfac = np.int64(aux_arrays[7])
        fig, ax = plt.subplots()
        colorWheel = {0:'or', 1:'ob', 2:'ok', 3:'oc', 4:'om', 5:'og'}

        for kk in range(5):
            idx = np.where(stretchfac == kk)[0]
            if (len(idx) > 0):
                ax.plot(np.real(allzlocs)[idx], np.imag(allzlocs)[idx], colorWheel[kk])
        idx = np.where(stretchfac > 4)[0]
        if (len(idx) > 0):
            ax.plot(np.real(allzlocs)[idx], np.imag(allzlocs)[idx], colorWheel[5])
        idx = np.where(stretchfac == -1)[0]
        if (len(idx) > 0):
            ax.plot(np.real(allzlocs)[idx], np.imag(allzlocs)[idx], '+')
        idx = np.where(stretchfac == -2)[0]
        if (len(idx) > 0):
            ax.plot(np.real(allzlocs)[idx], np.imag(allzlocs)[idx], '+r')
        if debugLevel == 2:
            plt.pause(0.5)
        if debugLevel == 3:
            print "Click anywhere in window to proceed"
            plt.waitforbuttonpress()
        if debugLevel == 4:
            fig.canvas.mpl_connect('key_press_event', on_key_event)    
            # Disable default Matplotlib shortcut keys:
            keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
            for key in keymaps:
                plt.rcParams[key] = ''
            ax.set_title('Interactive Window, select a vowel on the keyboard to close window and proceed:')
            plt.show()
        plt.close()
        
        print "Debugging on in fill_image_gaps"
                
    return src_list

def test_point_approx(sx, sy, sr, lens):
    """ Determine whether to use point source approximation
        If source center or its edge is too close to the caustics
        dont use point source approximation
    """
    # Get caustic locations from lens object
    dx = (sx - np.real(lens.caustLocs))
    dy = (sy - np.imag(lens.caustLocs))
    mindist = np.min(np.sqrt(dx*dx + dy*dy) - sr)
    usePoint = False
    if mindist > 8.0*max([sr, 1.0e-2]):
        usePoint = True
    return usePoint

def on_key_event(event):
    '''Keyboard interaction
    http://central.scipy.org/item/84/1/simple-interactive-matplotlib-plots
    courtesy thomas haslwanter'''

    #print('you pressed %s'%event.key)        
    key = event.key

    # In Python 2.x, the key gets indicated as "alt+[key]"
    # Bypass this bug:
    if key.find('alt') == 0:
        key = key.split('+')[1]

    curAxis = plt.gca()
    if key in 'aeiou':
        curAxis.set_title('Well done!')
        plt.pause(0.01)
        plt.close()
    else:
        curAxis.set_title(key + ' is not a vowel: try again to find a vowel ....')
        plt.draw()

def get_time_step_mag(tsobj, debugLevel=0):
    """Get the magnification of a particular time step
       It proceeds by determining whether to use point source
       approximation or finite source.  If finite source
       make the base srcn source vertices.  Call fill_image_gaps()
       to make sure the resolution of the source outline is sufficient.
       Then, call finite_source_mag() to calculate finite source magnification
       INPUTS:
         tsobj - single_time_step object
       OUTPUTS:
         curmag - magnification over all images at the time step
         magtyp - [0/1] 0=point source approx used; 1=finite source used
    """
    srcx = tsobj.centerx
    srcy = tsobj.centery
    srcr = tsobj.radius
    lens = tsobj.lens
    srcn = tsobj.srcn
    
    # Decide if one should use point source or polygon approximation for 
    # calculating magnification
    use_point_approx = test_point_approx(srcx, srcy, srcr, lens)
    if (use_point_approx):
        obj = source_plane_vertex()
        obj.set_source_position(srcx, srcy, 0.0, 0.0)
        zlocs, jacs, paritys, zkappas, zkappaPs = \
                    solve_binary_roots(obj, lens)
        curmag = np.sum(1.0/np.abs(jacs))
        magtyp = 0
    else:
        #Do finite source calculation
        # loop over base source vertices
        source_locations = []
        srcAngs = np.linspace(0.0, 2.0*np.pi, srcn)
        
        for j in range(srcn):
            obj = source_plane_vertex()
            obj.set_source_position(srcx, srcy, srcr, srcAngs[j])
            zlocs, jacs, paritys, zkappas, zkappaPs = \
                        solve_binary_roots(obj, lens)
            obj.nImages = len(zlocs)
            obj.imageLocs = zlocs
            obj.imageParity = paritys
            obj.jacobians = jacs
            obj.kappas = zkappas
            obj.kappaDerivs = zkappaPs
            source_locations.append(obj)
        # Got the base source vertice images look for image locations
        #  That have large gaps between neighboring image locations
        # Then add source vertices to fill gaps
        source_locations = fill_image_gaps(source_locations, lens, debugLevel=debugLevel)
        # now calculate flux of images
        curmag = finite_source_mag(source_locations, debugLevel=debugLevel) /np.pi/srcr**2
        # For meeple source shape need to adjust area
        #  see source_plane_vertex for using meeple source shape
        #curmag = curmag / 0.0580416
        magtyp = 1
    return curmag, magtyp

    
# Run the test of the root finding
if __name__ == "__main__":
    # Assume positions and masses etc
    # Source coordinates in terms of Re
    debugLevel = 0
    strtOverride = 0
    # Turn off multiprocessin with multiProc = False
    multiProc = True
    
    NRandom=10
    for ii in range(NRandom):
        lensep = np.random.random_sample() * (1.2-0.1) + 0.1 # Binary lens separation
        lenm2 = 10.0**(np.random.random_sample() * (0.0+3.0) - 3.0) # mass of lens 2 scaled to mass of lens 1 (i.e. lenm1=1.0)
        # finite source size
        src_sz = 10.0**(np.random.random_sample() * (-0.8 + 3.0) -3.0)
        #lensep = 0.732948
        #lenm2 = 0.264296
        #src_sz = 0.147522
        str1 = 'Lens Sep: {:f}'.format(lensep)
        str2 =  'Lens2 Mass: {:f}'.format(lenm2)
        str3 = 'Source Sz: {:f}'.format(src_sz)
        print(str1, str2, str3)
        # Define the lens
        lens = lens_properties()
        # lens 1 at origin and lens 2 on x-axis
        lens.set_lens_properties(0.0, 0.0, 1.0, lensep, 0.0, lenm2)
    
        # Get a source track that goes through the caustics
        caustx = np.real(lens.caustLocs)
        causty = np.imag(lens.caustLocs)
        # Fit a line to the caustic curve points in the upper plane
        idx = np.where(causty >= 0.0)[0]
        lincoeffs = np.polyfit(caustx[idx], causty[idx], 1)
        lincoeffs[1] = lincoeffs[1] + 0.05
        srcSteps = 3000 # Number of source steps along source path
        src_x = np.linspace(-2.5, 2.5, srcSteps)
        src_y = np.polyval(lincoeffs, src_x)

        # Show source track and cuastics
        fig, ax = plt.subplots()
        fig.canvas.mpl_connect('key_press_event', on_key_event)
        # Disable default Matplotlib shortcut keys:
        keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
        for key in keymaps:
            plt.rcParams[key] = ''
        ax.plot(caustx, causty, '+r')
        ax.plot(src_x, src_y, '-k')
        ax.set_title('First, enter a vowel:')
        plt.show()
    

        # number of sides to approximate source
        # srcn is the base number of source vertices
        srcn = 13
        
        # Define the initial source vertice list        
        start = time.clock()
        all_time_steps = []
        for i in range(strtOverride, srcSteps):
            new_time = single_time_step()
            new_time.centerx = src_x[i]
            new_time.centery = src_y[i]
            new_time.radius = src_sz
            new_time.srcn = srcn
            new_time.lens = lens
            all_time_steps.append(new_time)
        # This is the main loop that moves the source relative to the binary lens
        if multiProc:
            pool = multiprocessing.Pool()
            resultslist = pool.map(get_time_step_mag, all_time_steps)
        else:
            resultslist = []
            for i in range(strtOverride, srcSteps):
                resultslist.append(get_time_step_mag(all_time_steps[i]))
        tmp = zip(*resultslist)
        allmags = np.array(tmp[0])
        magtyp = np.array(tmp[1])
        
        print (time.clock() -start)
        
        # Show magnitude curve
        fig, ax = plt.subplots()
        fig.canvas.mpl_connect('key_press_event', on_key_event)
    
        # Disable default Matplotlib shortcut keys:
        keymaps = [param for param in plt.rcParams if param.find('keymap') >= 0]
        for key in keymaps:
            plt.rcParams[key] = ''
        allidx = np.arange(len(allmags))
        ax.plot(allidx, allmags, '-k')
        idx = np.where(magtyp == 1)[0]
        if len(idx)> 0:
            ax.plot(allidx[idx], allmags[idx],'ok')
        idx = np.where(magtyp == 0)[0]
        if len(idx)> 0:
            ax.plot(allidx[idx], allmags[idx], 'or')
        ax.set_title('Interactive Window, select a vowel on the keyboard to close window and proceed:')
        plt.show()

