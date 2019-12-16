import numpy as np

# Order 2- Quadratic fitting >> Simpson's method

def poly_matrix(npt):
    x = np.linspace(-1,1,npt)
    mat = np.zeros([npt,npt])
    mat[:,0] = 1.0
    mat[:,1] = x

    if (npt>2):
        for i in range(1,npt-1):
            mat[:,i+1] = x*mat[:,i]
    return np.matrix(mat)

def area_coeffs(npt):
    vec = np.zeros([npt])
    for k in range(npt):
        vec[k] = (1-(-1)**(k+1))/(k+1)
    return np.matrix(vec)

def print_weights(orders):

    for odr in orders: 
        a = area_coeffs(odr+1)
        A_inv = np.linalg.inv(poly_matrix(odr+1))
        
        w =  a*A_inv
        print("Weights for order %d: "%(odr))
        print(w)
        print("Normalized Weights for order %d: "%(odr))
        print(w*(odr)/2)

        x = np.ones([odr+1])
        x = np.matrix(x)

        print("\nArea for unity input (validity check): %.1f"%(w*x.T)[0,0])
        print()

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
order = 2

'''
for i in range(1,10):
    print("Order %d, %2d points >> residual = %.16f"%(order,i*order+1,integrate(e,i*order+1,order)-soln))
'''

print("Order 2, 31 points >> residual = %.16f"%abs(integrate(e,31,2)-soln))
print("Order 4, 29 points >> residual = %.16f"%abs(integrate(e,29,4)-soln))
print("Order 6, 31 points >> residual = %.16f"%abs(integrate(e,31,6)-soln))
print("Order 8, 33 points >> residual = %.16f"%abs(integrate(e,33,8)-soln))




