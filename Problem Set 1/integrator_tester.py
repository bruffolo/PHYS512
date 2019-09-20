import numpy as np
from scipy.integrate import quad
from Q3 import optimized_integrater,simple_integrater


sig=0.1

# Test functions
function_1 = lambda x:1.0/(1.0+x**2)
function_2 = lambda x:1.0+np.exp(-0.5*x**2/(sig**2))
function_3 = lambda x:np.cos(x)

print("f(x) = 1/(1+x^2)")
f,err_estimate,neval=optimized_integrater(function_1,-10,10,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the OPTIMIZED integrator")
f,err_estimate,neval=simple_integrater(function_1,-10,10,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the SIMPLE integrator")
pred=quad(function_1, -10,10)
print("Integral computed by Scipy quad = ",pred[0], "with estimated error ", pred[1])

print()

print("f(x) = 1+e^{-0.5*(x/0.1)^2}")
f,err_estimate,neval=optimized_integrater(function_2,-1,1,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the OPTIMIZED integrator")
f,err_estimate,neval=simple_integrater(function_2,-1,1,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the SIMPLE integrator")
pred=quad(function_2, -1,1)
print("Integral computed by Scipy quad = ",pred[0], "with estimated error ", pred[1])

print()

print("f(x) = cos(x)}")
f,err_estimate,neval=optimized_integrater(function_3,-1,1,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the OPTIMIZED integrator")
f,err_estimate,neval=simple_integrater(function_3,-1,1,1e-7)
print("Computed Integral = ",f," with estimated error of ",err_estimate," and ",neval," function evaluations on the SIMPLE integrator")
pred=quad(function_3, -1,1)
print("Integral computed by Scipy quad = ",pred[0], "with estimated error ", pred[1])
