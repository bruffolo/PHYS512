import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


try:
    sim_name = sys.argv[1]
except(IndexError): 
    print("Enter a one of the following command line arguments: ")
    print("> periodic\n> isolated")
    exit()


npix = 128
directory = "Data/"+sim_name

x = np.loadtxt(directory+"/x.txt",delimiter=",");print("> x data loaded.")
y = np.loadtxt(directory+"/y.txt",delimiter=",");print("> y data loaded.")
z = np.loadtxt(directory+"/z.txt",delimiter=",");print("> z data loaded.")
E = np.loadtxt(directory+"/E.txt",delimiter=",");print("> E data loaded.")

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
    text = 'Energy Conservation:  %.2f %%'%((E[t]/E[0])*100)
    axs.text2D(-0.12, 1.05, text,fontsize = 8, transform=axs.transAxes)
    print("Frame: %d/%d\r"%(t,size), end = '')

    plt.pause(0.001)
print()