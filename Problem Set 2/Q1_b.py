import numpy as np
from matplotlib import pyplot as plt
from numpy.polynomial import chebyshev

plot = True

c = True
while (c): 
    try:
        print("Enter a number t greater than 0: ")
        t = input("t = ")
        t = float(t)
        if (t< 0):
                raise ValueError
        c = False
    except ValueError:
        print("Try again.")

########################
#######  Setup   #######
########################

# Setup the test points
x = np.linspace(0.5,1.0,100)

# Data points
y = np.log2(x)

# Order of the polynomial
npoly = 8

################################
#######  Chebyshev Fit   #######
################################

# Number of data points
ndata=x.size

# Model matrix A
A=np.zeros([ndata,npoly])

# Fill A with the Chebyshev polynomials
A[:,0]=1.0
A[:,1]=x
for i in range(1,npoly-1):
    A[:,i+1]=2*x*A[:,i]-A[:,i-1]

A = np.matrix(A) # Convert A to numpy matrix
d = np.matrix(y).transpose()

# SVD decomposition
u,s,vt = np.linalg.svd(A,False)
sinv = np.matrix(np.diag(1.0/s)) #s comes back as a 1-d array, turn it into a 2-d matrix

cheb_fit_params = vt.transpose()*sinv*(u.transpose()*d) # Best fit Chebyshev coeffs.
cheb_fit = A*cheb_fit_params # Best fit truncated Chebyshev poly.

##########################################
#######  Standard polynomial Fit   #######
##########################################

# Number of data points
ndata=x.size

# Model matrix A
A=np.zeros([ndata,npoly])

# Fill A with the standard polynomials
A[:,0]=1.0
for i in range(0,npoly-1):
    A[:,i+1] = A[:,i]*x

A = np.matrix(A) # Convert A to numpy matrix
d = np.matrix(y).transpose()

# SVD decomposition
u,s,vt = np.linalg.svd(A,False)
sinv = np.matrix(np.diag(1.0/s)) #s comes back as a 1-d array, turn it into a 2-d matrix

fit_params = vt.transpose()*sinv*(u.transpose()*d) # Best fit Chebyshev coeffs.
poly_fit = A*fit_params # Best fit truncated Chebyshev poly.

##############################
#######   PLOTTING     #######
##############################

####################################
#######   Return log2(t)     #######
####################################

# Break down the input number into x = y1*2**y2, where y1 is in [0,1].
y1, y2 = np.frexp(t)

# Take the log2(x)
answer = np.polynomial.chebyshev.chebval(y1,cheb_fit_params)[0] + y2

# Print the result to the user
print("log2(t) = ", answer )