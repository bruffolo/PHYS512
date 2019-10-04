import numpy as np
from matplotlib import pyplot as plt
from numpy.polynomial import chebyshev

plot = True

c = True
while (c): 
    try:
        print("Enter a number t between 0.5 and 1.0: ")
        t = input("t = ")
        t = float(t)
        if (t< 0.5 or t>1.0 ):
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

if(plot):
    # Main figure
    fig1 = plt.figure(1) # Add main frame
    frame1=fig1.add_axes((.13,.3,.83,.6))

    plt.plot(x,y,'o',markersize = 3) # plot data
    plt.plot(x,poly_fit, color = "Orange", label = " $%d^{th}$ order polynomial fit"%npoly) # plot polynomial fit 
    plt.plot(x,cheb_fit,'--',color = "Red", label = "Truncated Chebyshev fit") # plot Chebyshev fit 

    frame1.set_xticklabels([]) # remove x-tic labels for the first frame
    plt.ylabel("$log_{2}(x)$")
    plt.grid()
    plt.legend()

    # Residual figure
    frame2=fig1.add_axes((.13,.05,.83,.2))  # Add frame for Residuals
    y = np.matrix(y).T
    diff_chev = cheb_fit - y # Compute Chebyshev residuals
    diff_poly = poly_fit - y # Compute standard polynomial residuals
    plt.plot(x,diff_poly*10**6,'or', color = "Orange",markersize = 4) # plot the residuals
    plt.plot(x,diff_chev*10**6,'or', color = "red", markersize = 2, label = "Chebyshev Residuals") # plot the residuals
    plt.ylabel("Residuals (x$10^{-6}$)")

    plt.grid()
    diff_chev = np.abs(diff_chev)
    diff_poly = np.abs(diff_poly)
    #print("Maximum error of the Truncated Chebychev polynomial fit on the interval [0.5,1] is: ",diff_chev.max())
    #print("Maximum error of the standard polynomial fit on the interval [0.5,1] is: ",diff_poly.max())

    # Save the plot
    plt.savefig("log2(x)")

####################################
#######   Return log2(t)     #######
####################################

print("log2(t) = ", np.polynomial.chebyshev.chebval(t,cheb_fit_params)[0])