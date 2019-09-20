#########################################
##### Assignment #1: Problem #3    ######
#########################################
import numpy as np

def optimized_integrater(func,a,b,tol,left_end = None,mid= None,right_end = None,recursion = False):

    # Break up the interval (b-a) into 4 pieces of equal size
    x=np.linspace(a,b,5)

    # Evaluate the function each of the 5 data points 
    y=func(x)

    if (recursion == False):
        # Keep track of function evaluations
        neval = len(x) 
        
        # Integral estimate using Simspson's rule with n = 2. 
        int_1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
        
        # Integral estimate using Simspson's rule with n = 4.
        int_2=((y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0)*(b-a)
    else:
        # Keep track of function evaluations
        neval = len(x)-3
        
        # Integral estimate using Simspson's rule with n = 2. 
        int_1=(left_end+4*mid+right_end)/6.0*(b-a)
        
        # Integral estimate using Simspson's rule with n = 4.
        int_2=((left_end+4*y[1]+2*mid+4*y[3]+right_end)/12.0)*(b-a)
        


    # Compute error estimation
    err = np.abs(int_2-int_1)

    # If error is less than specified tolerance return the computed integral to the user.
    if (err < tol):

        # Return the integral, error estimate, and number of function evaluations
        return (16.0*int_2-int_1)/15.0,err,neval

    # If not, split the interval 
    else:
        # Split the interval (b-a) in half
        mid=0.5*(b+a)

        # Evaluate the integral of the left hand side of the midpoint
        int_left,err_left,neval_left = optimized_integrater(func,a,mid,tol/2.0,y[0],y[1],y[2],True)

        # Evaluate the integral of the right hand side of the midpoint
        int_right,err_right,neval_right = optimized_integrater(func,mid,b,tol/2.0,y[2],y[3],y[4],True)

        ################################################################################################
        ### NOTE: since the interval is half as large as the original, we divide the tolerance by 2. ### 
        ################################################################################################

        # Keep track of the number function evaluations
        neval=neval+neval_left+neval_right

        # Sum the integrals
        int_tot =int_left + int_right

        # Sum the errors
        err=err_left+err_right

        # Return the integral, error estimate, and number of function evaluations
        return int_tot,err,neval


# Simple integrator from class lectures
def simple_integrater(fun,a,b,tol):
    x=np.linspace(a,b,5)
    #np.median(np.diff(x))
    y=fun(x)
    neval=len(x) #let's keep track of function evaluations
    f1=(y[0]+4*y[2]+y[4])/6.0*(b-a)
    f2=(y[0]+4*y[1]+2*y[2]+4*y[3]+y[4])/12.0*(b-a)
    myerr=np.abs(f2-f1)
    #print([a,b,f1,f2])
    if (myerr<tol):
        #return (f2)/1.0,myerr,neval
        return (16.0*f2-f1)/15.0,myerr,neval
    else:
        mid=0.5*(b+a)
        f_left,err_left,neval_left=simple_integrater(fun,a,mid,tol/2.0)
        f_right,err_right,neval_right=simple_integrater(fun,mid,b,tol/2.0)
        neval=neval+neval_left+neval_right
        f=f_left+f_right
        err=err_left+err_right
        return f,err,neval
