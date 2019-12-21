import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os

wd = os.getcwd()

##########################################
###          Orbit Animator            ###
##########################################

try:
    file_ext = sys.argv[1]
except(IndexError): 
    print("Enter a one of the following command line arguments: ")
    print("> orbit\n> orbit_complex")
    exit()

x_file = "x_"+file_ext+".txt"
y_file = "y_"+file_ext+".txt"
z_file = "z_"+file_ext+".txt"


npix = 64

x = np.loadtxt(wd+"/Positional_data/"+x_file,delimiter=",");print("> x data loaded.")
y = np.loadtxt(wd+"/Positional_data/"+y_file,delimiter=",");print("> y data loaded.")
z = np.loadtxt(wd+"/Positional_data/"+z_file,delimiter=",");print("> z data loaded.")


npix = 64
size = len(x)

N_trajectories = len(x[0,:])

fig = plt.figure(1)
axs = plt.axes(projection='3d')
axs.set_axis_off()

colors = ["orange","limegreen","blue","red"]

for t in range(size):
    axs.clear()
    axs.set_axis_off()

    for i in range(1,N_trajectories):
        axs.plot3D(x[:t,i],y[:t,i],z[:t,i],color = colors[i-1])
        axs.scatter3D(x[t,i],y[t,i],z[t,i],color = colors[i-1])       

    axs.scatter3D(x[t,0],y[t,0],z[t,0],color = "Black")
    
    axs.set_xlim([0,npix]);axs.set_ylim([0,npix]);axs.set_zlim([0,npix])
    plt.pause(0.0001)
    print("Frame: %d/%d\r"%(t,size), end = '')
