""" This file contains the input parameters needed to run multizone.py.

Below are examples for the parameters needed:

# Plot window settings:
# Format:   [nu_min, nu_max, nuF_min, nuF_max]
plotrange = [1e9,    1e26,   1e-13,   1e-9]

# Parameters for creating the energy vector. 
# nu_min (real) : minimum frequency for the plot
# nu_max (real) : maximum frequency for the plot
# n_nu   (int)  : number of frequency bins  
nu_min = 1E8  ;  nu_max = 1E28  ;  n_nu = 200

# Location of the data file
datafile = "./data/3C279.dat"

# Model fits. 
# Each on its own line, inside brackets, ending with comma.
# Parameters: [ alpha, doppler, B, gamma, Norm ]
SED_fit_parameters = [  #     You can give one or more parameter lines.
    [ 1.0, 8.0, 0.2, 1.0E4, 1.0E-28 ],   # Emission zone 1
    [ 1.0, 8.0, 0.2, 2.0E4, 1.0E-27 ],   # Emission zone 2
    [ 1.0, 8.0, 0.2, 3.0E4, 1.0E-28 ],   # Emission zone 3
    ] # This ending bracket needs to be here.

"""
##############################################################################

# Plot window settings. [nu_min, nu_max, nuF_min, nuF_max]
plotrange = [1e9, 1e24, 1e-12, 1e-9]

# Parameters for creating the energy vector
nu_min = 1E8  ;  nu_max = 1E28  ;  n_nu = 200

##############################################################################

## Test case
src_name = "Example"
datafile = "./data/example.dat"
SED_fit_parameters = [ 
#    alpha, doppler, B,  gamma,  Norm
    [ 0.6,  11.0,  0.1,  1.5E3,  1.0E-27 ], # Emission zone 1
    [ 1.0,  11.0,  0.1,  8.0E3,  8.0E-27 ]  # Emission zone 2
]
