import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Animator
x_file = "x.txt"
y_file = "y.txt"
z_file = "z.txt"


npix = 64

x = np.loadtxt("Positional_data/"+x_file,delimiter=",");print("> x data loaded.")
y = np.loadtxt("Positional_data/"+y_file,delimiter=",");print("> y data loaded.")
z = np.loadtxt("Positional_data/"+z_file,delimiter=",");print("> z data loaded")

size = len(x)

plt.ion()
fig = plt.figure(1)
axs = plt.axes(projection='3d')
axs.set_axis_off()


for t in range(size):
    axs.clear()
    axs.set_axis_off()

    axs.scatter3D(x[t],y[t],z[t],s = 2)
    axs.set_xlim([0,npix]);axs.set_ylim([0,npix]);axs.set_zlim([0,npix])
    
    print("Frame: %d/%d\r"%(t,size), end = '')

    plt.pause(0.05)