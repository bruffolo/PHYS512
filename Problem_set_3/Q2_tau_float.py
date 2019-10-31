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

def get_cmb(params,lmax=2000):
    
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
    
    return cmb[:1199]

def calc_grad(p):
    y0 = get_cmb(p)

    # Small increment used to calculate the derivative of the get_cmb function with 
    # respect to each parameter
    dx = p0/10000#np.array([0.00002,0.00002,0.0001,2e-12,0.00096])
    grad = np.zeros([y0.size,p.size])

    grad[:,0]= ( get_cmb(p+[dx[0],0,0,0,0,0]) - y0 )/dx[0]
    grad[:,1]= ( get_cmb(p+[0,dx[1],0,0,0,0]) - y0 )/dx[1]
    grad[:,2]= ( get_cmb(p+[0,0,dx[2],0,0,0]) - y0 )/dx[2]
    grad[:,3]= ( get_cmb(p+[0,0,0,dx[3],0,0]) - y0 )/dx[3]
    grad[:,4]= ( get_cmb(p+[0,0,0,0,dx[4],0]) - y0 )/dx[4]
    grad[:,5]= ( get_cmb(p+[0,0,0,0,0,dx[5]]) - y0 )/dx[5]
    
    return y0,grad


# Import the data
a = np.loadtxt("TT_power_spectra.txt").T

# Unpack the data
x,y,error,col4,col5 = a # Ignore col4 and col5

p0 = np.asarray([65,0.02,0.1,0.05,2e-9,.96]) #np.asarray([45,0.04,0.1,2e-9,0.86])

# Load in the initial guess at the best fit parameters
p = p0.copy()

plt.figure()
plt.plot(x,y,'o', markersize = 4)

# Calculate the chisquare and step for the initial guess.

# Calculate and plot the initial guess 
#pred,grad = calc_grad(p,x)
#plt.plot(x,pred,label = "Initial guess")
print()

# Gauss-Newton/Levenberg-Marquardt loop initialization.

chisq = 0   #
λ0 = 1e-2      # Initial value for "Damping" parameter
cycle = 0   # Keep track of the number of cycles of the Gauss-Newton algorithim.
ε = 1e-1

λ = λ0
W = np.diag(1/error**2)

d_chisq = 1

# Impliment the Gauss-Newton algorithm 
while(d_chisq > 0.001):
    cycle += 1

    pred,J = calc_grad(p)
    r = y-pred
    err = (r**2).sum()
    chisq = (((y-pred)/error)**2).sum()

    r = np.matrix(r).T
    J = np.matrix(J)

    lhs = J.T*W*J + λ*np.matrix(np.diag(J.T*W*J))
    rhs = J.T*W*r

    # Calculate parameter step
    dp = np.linalg.inv(lhs)*(rhs)

    ####
    h = np.zeros(len(p))
    for i in range(p.size):
        h[i] = p[i]+dp[i]
    D = get_cmb(h)
    new_chisq = (((y-D)/error)**2).sum()
    d_chisq = chisq - new_chisq

    # Calculate metric for step evaluation
    f = dp.T*(λ*np.matrix(np.diag(J.T*W*J)*dp+J.T*W*r))
    ρ = (d_chisq)/np.abs(f)
    
    if(ρ>ε):
        for i in range(p.size):
            p[i] = p[i]+dp[i]
        λ = max(λ/9,10**(-7))

    else:
         λ = min(λ*11,10**(7))
        
    if(ρ>ε):
        print("Cycle %d: χ^2 = %f"%(cycle,new_chisq))
    else:
        print("Cycle %d: χ^2 = %f"%(cycle,chisq))

print("\nFinal parameter set: ")
print("p = ",p)
sigma_p = np.sqrt(np.diag(np.linalg.inv(J.T*W*J)))
print("σ_p = ",sigma_p)

plt.plot(x,pred, label= "Levenberg–Marquardt")    
plt.legend()
plt.savefig("Q2_tau_float.png")
plt.show()