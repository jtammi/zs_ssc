#!/usr/bin/env python
"""zs_ssc: Zacharias & Schlickeiser (2012) SSC model. 

Purpose:      Multi-zone SED plotting with synchrotron and SSC.
Requirements: Parameters in separate file zs_params.py.
Version:      v.1.0 [20131029] 
Author:       Joni Tammi  (joni.tammi@aalto.fi)

Details: Physical model and equations described in the paper 
"A new ordering parameter of spectral energy distributions from synchrotron 
self-Compton emitting blazars", by M. Zacharias & R. Schlickeiser, in 
Monthly Notices of the Royal Astronomical Society, Vol. 420, p. 84 (2012)
http://mnras.oxfordjournals.org/content/420/1/84

N.b. the original article contains some typos in the equations; they are 
corrected in this code based on private communication between J.T. and M.Z.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys

# Input parameters
import zs_params as zs

##############################################################################

class JetComponent:
    """Class for jet components."""

    def __init__(self, id):
        self.__id = id

##############################################################################

class Calculated_SED:
    """Class for simulated SEDs."""

    def __init__(self,sed_id,params,cnu):
        self.__id = sed_id
        self.__nu = cnu
        self.__nuF = []
        self.__sync = []
        self.__ssc = []
        self.__alpha = params[0]
        self.__doppler = params[1]
        self.__b = params[2]
        self.__gamma_0 = params[3]
        self.__Norm = params[4]
        self.__ls = ' '

        if sed_id == -1:
            self.__ls = '-k'
        elif sed_id == 0:
            self.__ls = '--r'
        elif sed_id == 1:
            self.__ls = '--c'
        elif sed_id == 2:
            self.__ls = '--g'
        elif sed_id == 3:
            self.__ls = ':r'
        elif sed_id == 4:
            self.__ls = ':c'
        elif sed_id == 5:
            self.__ls = ':g'
        else:
            self.__ls = ':k'

    def calculate_SED(self):
        alpha   = self.__alpha
        doppler = self.__doppler
        b       = self.__b
        gamma_0 = self.__gamma_0
        Norm    = self.__Norm
        cnu     = self.__nu

        ## Conveniencies
        gamma_4 = gamma_0 / 1.0e4
        d4 = doppler**4.0
        K = 0.136 * b * gamma_4**3 # Klein-Nishina parameter, Eq. (44.)

        # Frequencies needed in the calculations
        nu_syn = 4.1E14 * doppler * b * gamma_4**2
        nu_T   = 1.682E23 * doppler * b * gamma_4**4
        nu_B   = 2.3E24 * doppler * b**(-1./3.)
        nu_g0  = 1.2E24 * doppler * gamma_4
        nu_c = 2.9E14 * doppler * b * gamma_4**2 / alpha**2

        if alpha < 1.0:
            f_synch  = 5.9e12 * d4 * alpha**2 / gamma_4 * \
                (cnu / nu_syn)**(1./2.)* np.exp(-cnu / nu_syn)
            if K < 1:
                f_ssc = 3.0e13 * d4 * alpha**4 / gamma_4 * \
                    (cnu / nu_T  )**(3./4.) * np.exp(-cnu / nu_T)
            else:
                f_ssc = 2.2e14 * d4 * alpha**4 / gamma_4**4 / \
                    b * (cnu / nu_B  )**(3./4.) / \
                    (1.+ cnu / nu_B)**(7./4.) * \
                    ( 1.0 - (cnu / nu_g0)**(7./12.) )
        else:
            f_synch  = 5.9e12 * d4 * alpha**2 / gamma_4 * \
                (cnu / nu_syn)**(1./2.)\
                / (1.0 + cnu/nu_c)  * np.exp(-cnu / nu_syn)
            if K < 1:
                f_ssc = 3.0e13 * d4 * alpha**4 / gamma_4 * \
                    (cnu / nu_T  )**(3./4)\
                    / (1.0 + alpha**4 * cnu / nu_T)**(1./2) * \
                    np.exp(-cnu / nu_T)
            elif K < alpha**3: # i.e., 1 < K < alpha**3
                f_ssc = 3.0e13 * d4 * alpha**4 / gamma_4   * \
                    (cnu / nu_T )**(3./4.)
                gt_indices = np.where( cnu > nu_T / alpha**4 )
                f_ssc[gt_indices] = 3.1e13 * d4 * alpha**2 / \
                    (gamma_4**2 * b**(1./3.)) * \
                    (cnu[gt_indices] / nu_B )**(1./4.) / \
                    (1.0 + cnu[gt_indices] / nu_B )**(13./4.) * \
                    (1.0-(cnu[gt_indices]/nu_g0)**(13./3.))
            else: # i.e., K > alpha**3
                f_ssc = 2.2e14 * d4 * alpha**4 / (gamma_4**4 * b)   * \
                    (cnu / nu_B )**(3./4.) / (1.+ cnu / nu_B)**(7./4.)
                gt_indices = np.where( cnu > nu_g0 / alpha ) # nu > values
                f_ssc[gt_indices] = 3.1e13 * d4 * alpha**2 / \
                    (gamma_4**2 * b**(1./3)) * \
					(cnu[gt_indices] / nu_B )**(-3.)\
                    * ( 1.0 - (cnu[gt_indices]/nu_g0)**(13./3) )

        self.__sync = f_synch * Norm
        self.__ssc = f_ssc * Norm
        self.__nuF = self.__sync + self.__ssc

    def add_to_SED(self,SED):
        #print np.shape(SED), np.shape(self.__nuF)
        #total = np.add(self.__nuF, SED)
        if np.size(self.__nuF) < 1:
            self.__nuF = SED
        else:
            self.__nuF = np.add(self.__nuF, SED)

    def get_SED(self):
        return self.__nuF

    def get_max(self):
        return np.max(self.__nuF)

    def plot_SED(self):
        #print self.__id, self.__ls
        plt.loglog( self.__nu, self.__nuF, self.__ls )

##############################################################################

class Observations:
    """Class or all data-related things.
    Use get_nu(), get_nuF() and get_err() to get the nu, nu*F, and nu_F_err."""

    def __init__(self,datafile):
        self.__file = datafile
        self.__name = ""
        self.__rawdata = []
        self.__nu = []
        self.__nuF = []
        self.__err = []
        self.__W = []
        self.__yerr = []
        self.__maxY = 0.0
        self.__minY = 0.0

    def read_data(self):
        """Reads the data from the file."""
        data = np.genfromtxt(self.__file) # Planck SED
        self.__nu = 10.0**data[:,0]
        self.__nuF = 10.0**data[:,2]
        self.__err = 10.0**data[:,3]
        #self.__W = 10.0**data[:,4]
        self.__yerr = [ self.__nuF - self.__nuF / self.__err, \
                            self.__nuF * self.__err - self.__nuF  ]
        self.__maxY = max( self.__nuF )
        self.__minY = min( self.__nuF )
    
    def get_ymax(self):
        return self.__maxY
    
    def get_ymin(self):
        return self.__minY
    
    def get_name(self):
        return self.__name

    def get_nu(self):
        return self.__nu

    def get_nuF(self):
        return self.__nuF

    def get_err(self):
        return self.__err

    def plot_observations(self):
        plt.errorbar(self.__nu, self.__nuF, yerr=self.__yerr, fmt='.b');

##############################################################################

def get_cnu(nu_min, nu_max, n_nu):
    """Create the logarithmic grid."""
    ## Frequency grids; border, difference, centre [b, d, c]
    bnu = nu_min * (nu_max/nu_min)**(np.arange(n_nu+1)/float(n_nu))
    dnu = bnu[1:] - bnu[0:n_nu]
    cnu = np.sqrt( bnu[1:] * bnu[0:n_nu] )
    return cnu

def run_fits(datafile, plotrange, SED_fits, cnu, src_name):
    """This is the actual functional main program."""

    n_SEDs = np.shape(SED_fits)[0]  # How many SED fits we are making

    # Create the observational data object and read in the data
    obs_data = Observations(datafile)
    obs_data.read_data()

    all_SEDs = []
    for i in range(0, n_SEDs):
        # Create a new SED fit, store it the list of all SEDs
        all_SEDs.append( Calculated_SED(i,SED_fits[i], cnu) )
        all_SEDs[i].calculate_SED()
    
    # Create a total SED, and sum up the fluxes of all SEDs
    sum_SED = Calculated_SED(-1, SED_fits[0], cnu)
    for i in range(0, n_SEDs):
        sum_SED.add_to_SED(all_SEDs[i].get_SED())

    ##-------------------------------------------------------------
    # Plotting
    plt.clf()
    
    # Plot the data points
    obs_data.plot_observations()

    # Plot the SEDs from each components
    for i in range(0, n_SEDs):
        all_SEDs[i].plot_SED()

    # Plot the total SED
    sum_SED.plot_SED()

    plt.title(src_name)
    plt.xlabel(r'$\nu$')
    plt.ylabel(r'$\nu F_\nu$')

    # Increase vertical plot range if data doesn't fit in the window 
    new_ymin = min( plotrange[2], obs_data.get_ymin() )
    new_ymax = max( plotrange[3], obs_data.get_ymax(), sum_SED.get_max() )

    if (new_ymin < plotrange[2]):
        print "Lowering the y_min to fit all data inside the window."
        plotrange[2] = 10.0**np.floor(np.log10(new_ymin))
    
    if (new_ymax > plotrange[3]):
        print "Increasing the y_max to fit all data inside the window."
        plotrange[3] = 10.0**np.ceil(np.log10(new_ymax))

    plt.axis(plotrange)

    # Save images in files and show the plot on the screen.
    plt.savefig('pics/'+src_name+'.png')
    plt.savefig('pics/'+src_name+'.eps')
    #    plt.savefig('pics/'+src_name+'.png')
    #    plt.savefig('pics/'+src_name+'.eps')
    print "Images (eps and png) saved in directory pics"
    plt.show() 

##############################################################################

def main():

    # Ask for input file if one is not given
    if len(sys.argv) > 1: # or sys.argv[1] == "-h":
        print ""
        print "Usage:   python zs_ssc.py"
        print "Give fit and data parameters in the file zs_params.py"
        print ""
        sys.exit()
    else:
        print "Close plotting window to exit."
        cnu = get_cnu(zs.nu_min, zs.nu_max, zs.n_nu)
        run_fits(zs.datafile, zs.plotrange, zs.SED_fit_parameters,\
                 cnu,zs.src_name )

				 
main()
