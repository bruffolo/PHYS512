#########################################
##### Assignment #1: Problem #1    ######
#########################################

import numpy as np
from matplotlib import pyplot as plt

'''
The four point derivative 


'''

# 3 numerical derivative schemes 
def center_derivative(f,x,h):
    return (f(x+h)-f(x-h))/(2*h)

def forward_derivative(f,x,h):
    return (f(x+h)-f(x))/(h)

def four_point_derivative(f,x,h):
    return ( 8*(f(x+h)-f(x-h))-f(x+2*h)+f(x-2*h) )/(12*h)

# Computes relative error between computed derivative and the analytic sol
def error(computed,analytic):
    return np.abs((computed-analytic)/analytic)

# Test functions
func1 = lambda x:np.exp(x)
func2 = lambda x:np.exp(0.01*x)

# Analytic derivative values at x = 0 for func1 and func2
ans1 = 1
ans2 = 0.01

# Setting up the values of h we want tested.
num_orders = 18
init_order = -17

# numpy array to be filled with h values.
h = np.array([])

# Fill the h array with values.
# Fill the h array with values.
for i in range(num_orders):
    filler = np.linspace(10**(init_order+i),10**(init_order+i+1),100)
    h = np.concatenate((h,filler), axis = None)

plt.figure()
plt.yscale('log')
plt.xscale('log')

# Compute four point derivative with h array and plot
der = four_point_derivative(func1,0,h)
rel_err = error(der,ans1)
plt.plot(h,rel_err,color = "Orange", label = "f(x) = $e^{x}$")


# Setting up the values of h we want tested.
num_orders = 20
init_order = -17

# numpy array to be filled with h values.
h = np.array([])

k = 0 #index variable

# Fill the h array with values.
for i in range(num_orders):
    filler = np.linspace(10**(init_order+i),10**(init_order+i+1),100)
    h = np.concatenate((h,filler), axis = None)

# Compute four point derivative with h array and plot
der = four_point_derivative(func2,0,h)
rel_err = error(der,ans2)
plt.plot(h,rel_err,color = "lawngreen", label = "f(x) = $e^{0.01x}$")


# Add axis labels, legend and grid to the plot
plt.ylabel("Relative error")
plt.xlabel("Step size [h]")
plt.grid()
plt.ylim(1e-15)
plt.legend()

# Display the plot to user
plt.show()
