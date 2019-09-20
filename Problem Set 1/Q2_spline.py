###############################################
##### Assignment #1: Problem #2 part 1   ######
###############################################
import numpy as np
from matplotlib import pyplot as plt

def spline_interp(z,x,y):
    assert len(x) == len(y), "x and y must have the same dimensions, Dave."
    n = len(x)
    for i in range(n-1):
        if (x[i]<=z):
            if(z<=x[i+1]):
                t = [x[i-1],x[i],x[i+1],x[i+2]]
                g = [y[i-1],y[i],y[i+1],y[i+2]]

                interp = np.polyfit(t,g,3)
                p = np.poly1d(interp)
                p_prime = np.polyder(p)
                return p(z),p_prime(z)
    return 


# Import the data from .txt file
t,v,dv_dt = np.loadtxt("DT-670_Silicon_Diode_Specs.txt", unpack = True)
dv_dt = dv_dt/1000

func = lambda z: spline_interp(z,t,v)

# Ask the user for an input temperature and calculate the diode voltage
# at that temperature using the spline interpolation

c = True
while (c): 
    try:
        point = float(input("Enter a temperature in Kelvin: "))
        if(point > t[len(t)-1]):
            raise ValueError
        c = False
    except ValueError:
        print("\nPlease enter a temperature between %.1f K and %.1f K.\n" %(t[0],t[len(t)-1]))

# Get the interpolated value of the diode voltage at the temperature entered by the user.
interpolated_value = func(point)

# Print the result to the user.
print("The diode voltage at",point, "K is interpolated to be %.6f  V" %(interpolated_value[0]))

# Ask the user if they want to see a plot of the interpolation
c = True
while (c): 
    try:
        print("\nPlot the interpolation?")
        plot_interp = input("(y/n): ")
        if (plot_interp != 'y' and plot_interp != 'n'):
                raise ValueError
        c = False
    except ValueError:
        print("I'm sorry Dave, I cant do that.")

c = True
while (c): 
    try:
        print("\nPlot the derivative of the interpolation?")
        plot_derivative = input("(y/n): ")
        if (plot_derivative != 'y' and plot_derivative != 'n'):
                raise ValueError
        c = False
    except ValueError:
        print("I'm sorry Dave, I cant do that.")

if(plot_interp == 'y' or plot_derivative == 'y'):
    x = np.linspace(t[1],t[len(t)-2],100000)
    y = np.array([func(xi) for xi in x])
    y = y.transpose()
    
if(plot_interp == 'y'):
    
    # New figure
    plt.figure()

    # Plot the data
    plt.plot(t,v,'o',markersize = 3, color = "blue", label = "Measured Data", zorder = 1)

    plt.ylabel("Diode Voltage [V]")
    plt.xlabel("Temperature [K]")
    plt.title("Spline Interpolation")

    # Plot the spline interpolation
    plt.plot(x,y[0], color = "lawngreen", label = "Interpolation", zorder = 2)

    # Add the calculated point
    plt.plot(point,interpolated_value[0],'o', color = "Orange", markersize = 4, label = "T = %.1f K"%point)


    # Make the plot a little nicer 
    plt.xlim(1,501)
    plt.legend()
    plt.grid()
    

if(plot_derivative == 'y'):

    # New figure
    plt.figure()

    # Plot the measured derivative
    plt.plot(t,dv_dt*1000,'o',markersize = 3, color = "blue", label = "Measured derivative", zorder = 1)

    # Plot the derivative of the interpolation
    plt.plot(x,y[1]*1000, color = "red", label = "Interpolation derivative", zorder = 2)

    plt.ylabel("dV/dT [mV/K]")
    plt.xlabel("T [K]")

    plt.xlim(1,501)
    plt.legend()
    plt.grid()


# Show the plot(s) to user
plt.show()

