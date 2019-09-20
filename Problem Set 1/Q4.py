#########################################
##### Assignment #1: Problem #4    ######
#########################################

import numpy as np
from Q3 import optimized_integrater
from scipy.integrate import quad
from matplotlib import pyplot as plt

cos = lambda x:np.cos(x)
sin = lambda x:np.sin(x)

# Function to be integrated 
E = lambda theta,z_R: (0.5)*(sin(theta)*(z_R - cos(theta)))/(1+(z_R)**2 - 2*(z_R)*cos(theta))**(3/2.)

# Analytic solution to the integral
inverse_square =  lambda x:1/x**2 if x>1 else (0.5 if x==1 else 0)

E_z = lambda theta:E(theta,z_R=1)
#result,errors,neval=simple_integrate(E_z,0,np.pi,1e-7)
#print(quad(E_z, 0, np.pi))
#print(optimized_integrater(E_z,0,np.pi,1e-7))

plt.figure()

z = np.linspace(0,20,1000)


analytic = np.zeros([len(z)])
for i in range(len(z)):
   analytic[i] = inverse_square(z[i])
plt.plot(z,analytic,'o', markersize = 1, color = 'deepskyblue', label = "Analytic solution")
results = np.zeros(len(z))
errors = np.zeros(len(z))

for i in range(len(z)):
    E_z = lambda theta:E(theta,z_R=z[i])
    results[i],errors[i],neval=optimized_integrater(E_z,0,np.pi,1e-7)

plt.plot(z,results,"-", color = "orange", label = "Variable step size")


results = np.zeros(len(z))
errors = np.zeros(len(z))

for i in range(len(z)):
    E_z = lambda theta:E(theta,z_R=z[i])
    results[i],errors[i] = quad(E_z, 0, np.pi)
plt.plot(z,results,'--', color = "lawngreen", label = "Quad")

# Make the plot pretty
plt.legend()
plt.ylabel("E (Normalized)")
plt.xlabel("z/R")
plt.xlim(0,z[len(z)-1])

# Show the plot to user
plt.show()