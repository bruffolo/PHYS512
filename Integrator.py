import numpy as np

# Calculate polynomial matrix
def poly_matrix(npt):
    x = np.linspace(-1,1,npt)
    mat = np.zeros([npt,npt])
    mat[:,0] = 1.0
    mat[:,1] = x

    if (npt>2):
        for i in range(1,npt-1):
            mat[:,i+1] = x*mat[:,i]
    return np.matrix(mat)

# Calculate a_k (as defined in problem 3 as integral of x^k from -1 to 1)
def area_coeffs(npt):
    vec = np.zeros([npt])
    for k in range(npt):
        vec[k] = (1-(-1)**(k+1))/(k+1)
    return np.matrix(vec)

def print_weights(npoints):
    print("---------------------------------------------\n")
    for npt in npoints: 
        a = area_coeffs(npt)

        # Calculate inverse of polynomial matrix
        A_inv = np.linalg.inv(poly_matrix(npt))
        
        # Get weights
        w =  a*A_inv

        print("Weights for %d point fit: "%(npt))
        print(w)
        print("\nNormalized Weights for %d point fit: "%(npt))
        print(w*(npt-1)/2) # Normalized weights

        x = np.ones([npt])
        x = np.matrix(x)

        # Calculate the area for f(x) = 1 (should = 2)!
        print("\nArea for unity input (validity check): %.1f"%(w*x.T)[0,0])
        print("---------------------------------------------\n")

def integrate(f,npoints,order):

    if((npoints-1)%order != 0):
        print("Number of points submitted are incorrect for this order method.")
        print("Examples of correct number of points for method of order %d are: %d, %d, %d, ect.."%(order,order+1,2*order+1,3*order+1))
        return

    x = np.linspace(-1,1,npoints)
    x = f(x)
    dx = 2/(npoints-1)
    coeffs = area_coeffs(order+1)
    A_inv = np.linalg.inv(poly_matrix(order+1))
    
    # Calculate weights
    w = coeffs*A_inv 
    # Normalize the weights
    w = w*(order)/2

    tot = 0

    for i in range((npoints-1)//order):
        if(i==0):
            y = x[0:order+1]
        else:
            y = x[i*(order):i*(order)+order+1]
        y = np.matrix(y)
        #print(y)
        tot += (w*y.T)[0,0]*dx
    return tot

    
def e(x):
    return np.exp(x)

soln = np.e-1/np.e




