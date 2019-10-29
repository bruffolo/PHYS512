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
    
    H0,ombh2,omch2,As,ns = params

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
    
    return cmb[:1199]

def calc_grad(p,l):
    y0 = get_cmb(p)

    # Small increment used to calculate the derivative of the get_cmb function with 
    # respect to each parameter
    dx = p/10000
    grad = np.zeros([y0.size,p.size])

    grad[:,0]= ( get_cmb(p+[dx[0],0,0,0,0]) - y0 )/dx[1]
    grad[:,1]= ( get_cmb(p+[0,dx[1],0,0,0]) - y0 )/dx[1]
    grad[:,2]= ( get_cmb(p+[0,0,dx[2],0,0]) - y0 )/dx[2]
    grad[:,3]= ( get_cmb(p+[0,0,0,dx[3],0]) - y0 )/dx[3]
    grad[:,4]= ( get_cmb(p+[0,0,0,0,dx[4]]) - y0 )/dx[4]
    
    return y0[:l.size],grad[:l.size]


# Import the data
a = np.loadtxt("TT_power_spectra.txt").T

# Unpack the data
x,y,error,col4,col5 = a # Ignore col4 and col5

p0 = np.asarray([65,0.02,0.1,2e-9,0.96]) #np.asarray([65,0.02,0.1,2e-9,0.56])

# Load in the initial guess at the best fit parameters
p = p0.copy()

plt.figure()
plt.plot(x,y,'o', markersize = 4)

# Calculate the chisquare and step for the initial guess.

# Calculate and plot the initial guess 
pred,grad = calc_grad(p,x)
plt.plot(x,pred,label = "Initial guess")
print()

# Get ready for Gauss-Newton loop.
d_chisq = 1
chisq = 0

cycle = 0  # Keep track of the number of cycles of the Gauss-Newton algorithim.

# Impliment the Gauss-Newton algorithm 
while(np.abs(d_chisq) > 0.001):
    cycle += 1
    pred,grad = calc_grad(p,x)
    r = y-pred
    err = (r**2).sum()
    new_chisq = (((y-pred)/error)**2).sum()
    r = np.matrix(r).transpose()
    grad = np.matrix(grad)

    lhs = grad.transpose()*grad
    rhs = grad.transpose()*r
    
    # Calculate parameter step
    dp = np.linalg.inv(lhs)*(rhs)
    # Step each parameter
    for i in range(p.size):
        p[i] = p[i]+dp[i]
    print("Cycle %d: χ^2 = %f"%(cycle,new_chisq))
    print()
    d_chisq = new_chisq - chisq
    chisq = new_chisq
    


plt.plot(x,pred, label= "Gauss-Newton")    
plt.legend()
plt.savefig("Q2.png")
plt.show()