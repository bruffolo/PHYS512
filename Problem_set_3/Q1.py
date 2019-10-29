# Numpy, scipy for data manipulation
import numpy as np
import scipy as sc

# Matplotlib for plotting
from matplotlib import pyplot as plt

import camb
import corner

import json
config = json.load(open('matplotlib_config.json'))
plt.rcParams.update(config)


def get_cmb(params,tau = 0.05,lmax=2000):
    
    if (tau_fixed):
        H0,ombh2,omch2,As,ns = params
    else:
        H0,ombh2,omch2,tau,As,ns = params

    params = camb.CAMBparams()
    params.set_cosmology(H0=H0,ombh2=ombh2,omch2=omch2,mnu=0.06,omk=0,tau=tau)
    params.InitPower.set_params(As=As,ns=ns,r=0)
    params.set_for_lmax(lmax,lens_potential_accuracy=0)
    
    results = camb.get_results(params)
    powers = results.get_cmb_power_spectra(params,CMB_unit='muK')
    cmb = powers['total']
    cmb = cmb[:,0]  
    
    # Remove the monopole and dipole
    cmb = cmb[2:]
    
    return cmb


# Import the data
a = np.loadtxt("TT_power_spectra.txt").T

# Unpack the data
col1,col2,col3,col4,col5 = a # Ignore col4 and col5

tau_fixed = False

plt.figure()

# CMB parameters
pars = np.asarray([65,0.02,0.1,0.05,2e-9,0.96])

# Calculate the CMB from the parameters above.
cmb = get_cmb(pars)

# Plot the data
plt.plot(col1,col2,'.')

# Plot the calculated model
plt.plot(cmb, label = "Model")

# Calculate Chi-Sqaure
r = (col2-cmb[:1199])/col3
chisq = (r**2).sum()

# Display the Chi-Sqaure to the user
print("chi^2 = %f"%(chisq))

# Setup the plot
plt.xlim(0)
plt.legend()
plt.show()
