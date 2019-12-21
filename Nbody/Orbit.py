import numpy as np
from numpy.fft import rfftn,irfftn
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import ctypes
from numba import jit
import sys
import os

wd = os.getcwd()

# Pull necessary functions from companion C library
c_lib = ctypes.cdll.LoadLibrary("PM_methods.so")

# Get the methods for the isolated BC's
to_grid = c_lib.to_grid
to_grid.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.c_long]

handle_isolated_boundaries = c_lib.handle_isolated_boundaries
handle_isolated_boundaries.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.c_long]

def to_grid_py(xy,grid,cell_loc):
    g_xy = np.asarray(np.round(xy),dtype='int')
    n = xy.shape[0]
    for i in range(n):
        grid[g_xy[i,0],g_xy[i,1],g_xy[i,2]] = grid[g_xy[i,0],g_xy[i,1],g_xy[i,2]] + 1
        cell_loc[i] = g_xy[i,0],g_xy[i,1],g_xy[i,2]

# System Configuration
npix = 2**6        # Pixel number
npt  = 5           # Number of particles
ndim = 3           # Number of dimenesions

# Particle mass parameters
m_equal  = True # Set all masses to common value 
m_random = False  # Set masses of particles randomly

m_max = 1        # Maximum particle mass 
m_min = 0.4      # Mininmum particle mass

# Number of iterations
niter = 3000

# Time step
dt = 0.001

# Check if user wants verbose printing to stdout
try:
    if(sys.argv[1] == '-b'): verbose = True
except(IndexError): 
    verbose = False

# Plotting parameters (during simulation)
show3D = False
nth = 1        # Plot every nth particle on screen
plt_iter = 1    # Plot particles every plt_iter iteration

# Data saving parameters
save_data = True
save_iter = 1   # Save position particle data every save_iter iteration
save_nth  = 1  # Save position particle data every save_nth particle

x_file = "x_orbit_complex.txt"
y_file = "y_orbit_complex.txt"
z_file = "z_orbit_complex.txt"

# Benchmark mode
benchmark = False

# Perferably if benchmark = True, you should set niter = 1, otherwise benchmarks will print out every iteration!

# Potential test mode
potential_test = False

########################################
###        Green's function          ###
########################################
green_func = np.zeros([2*npix+1,2*npix+1,2*npix+1])

# Method to fill the Green's function array
@jit(nopython=True,parallel = True)
def get_H(H):
    for i in range(npix+1):
        for j in range(npix+1):
            for k in range(npix+1):
                if(i == 0 and j == 0 and k == 0 ):
                    # This value describes the strength of the force between particles in the same cell
                    H[i,j,k] = 4/3  # Chose wisely
                else:
                    H[i,j,k] = 1/np.sqrt(i**2 + j**2 + k**2)
        #print("Progress: %.2f %%\r"%((i/npix)*100), end = '')
    #print()
    for i in range(npix+1):
        for j in range(npix+1):
            for k in range(npix+1):
                h = H[i,j,k]
                H[ 2*npix-i ,    j     ,     k    ] = h
                H[ 2*npix-i , 2*npix-j ,     k    ] = h
                H[ 2*npix-i ,    j     , 2*npix-k ] = h
                H[ 2*npix-i , 2*npix-j , 2*npix-k ] = h
                H[     i    , 2*npix-j ,      k   ] = h
                H[     i    , 2*npix-j , 2*npix-k ] = h
                H[     i    ,    j     , 2*npix-k ] = h
        #print("Progress: %.2f %%\r"%((i/npix)*100), end = '')
    #print()

get_H(green_func) # Make the Green's function array
green_ft = rfftn(green_func[:-1,:-1,:-1]) # Take the FT of Green's function array

#plt.imshow(green_func[:,:,0]);plt.colorbar();plt.show()

grid = np.zeros([2*npix,2*npix,2*npix]) # Make the grid array


########################################
###         Particle setup           ###
########################################

# Positions
x = np.zeros([npt])
y = np.zeros([npt])
z = np.zeros([npt])

# Masses

if(m_equal):
    m = np.ones([npt],dtype=float)*m_max
if(m_random):
    m = (np.random.rand(npt)*m_max) + m_min


# Velocities
vx = np.zeros([npt])
vy = np.zeros([npt])
vz = np.zeros([npt])

# Accelerations
ax = np.zeros([npt])
ay = np.zeros([npt])
az = np.zeros([npt])

# Starting positions
x[0] = npix/2
y[0] = npix/2    
z[0] = npix/2    
m[0] = 100000

x[1] = npix*(8.5/10)
y[1] = npix/2
z[1] = npix/2
vy[1] = np.sqrt(m[0]/(x[1]-x[0]))


x[2] = npix*(1/10)
y[2] = npix/2
z[2] = npix/2
vy[2] = -np.sqrt(m[0]/(x[0]-x[2]))

x[3] = npix*(7.5/10)
y[3] = npix/2
z[3] = npix/1.8
vy[3] = -np.sqrt(m[0]/(x[3]-x[0]))

x[4] = npix*(6.5/10)
y[4] = npix/2
z[4] = npix/2
vy[4]  = np.sqrt(m[0]/(x[4]-x[0]))


if(potential_test):
    npix = 64
    niter = 0
    show3D = False
    z[:] = npix//2


# Data saving arrays (will be printed to file)
if(save_data):
    x_save = np.zeros([niter//save_iter,npt//save_nth])
    y_save = x_save.copy()
    z_save = x_save.copy()

if(show3D): plt.ion()

# Print statments for user
if(verbose):
    print("\nPixel density: %d x %d x %d"%(npix,npix,npix))
    print("Number of particles: %d"%npt)
    print("Numer of iterations: ",niter)
    print("Time step: ",dt)
    if(show3D):
        print("\n3D Plot parameters:")
        print(">Number of plotted particles: %d"%(npt/nth))
        print(">Plot refreshed every %dth iteration"%(plt_iter))
    if(save_data):
        print("\nData save parameters: ")
        print(">Positional data of every %dth particle will be saved every %dth iteration"%(save_nth,save_iter))
        print(">Data will be saved in files: ",x_file,", ",y_file,", ",z_file)
    print()

if(benchmark):
    print("\n####### Benchmarks #######")

########################################
###         Initializtion            ###
########################################

t91 = time.time() # Simulation start time

lx = np.asarray(np.round(x),dtype='int')
ly = np.asarray(np.round(y),dtype='int')
lz = np.asarray(np.round(z),dtype='int')

# This array will track out of bounds particles
bc = np.zeros([npt],dtype = 'int')

# Assign particles to the grid
to_grid(lx.ctypes.data, ly.ctypes.data, lz.ctypes.data, grid.ctypes.data, m.ctypes.data, bc.ctypes.data, npt, npix)
l = np.array([lx,ly,lz])


# Potential calculation
density_ft = rfftn(grid)
potential = irfftn(density_ft*green_ft)

# Acceleration calculation in every grid cell
Fx_prev,Fy_prev,Fz_prev = np.gradient(potential)

ax_prev = Fx_prev[tuple(l)]/m
ay_prev = Fy_prev[tuple(l)]/m
az_prev = Fz_prev[tuple(l)]/m

# Reset grid
grid.fill(0)

if(potential_test):
    plt.figure(1)
    plt.imshow(potential[:,:,npix//2].T,origin='lower');plt.colorbar()
    plt.show()

########################################
###            Main Loop             ###
########################################

if(show3D):
    fig = plt.figure(1)
    axs = plt.axes(projection='3d')
    axs.set_axis_off()
    axs.view_init(45,0 )

for t in range(niter):
    t11 = time.time()

    lx = np.asarray(np.round(x),dtype='int')
    ly = np.asarray(np.round(y),dtype='int')
    lz = np.asarray(np.round(z),dtype='int')

    # Assign particles to the grid
    to_grid(lx.ctypes.data, ly.ctypes.data, lz.ctypes.data, grid.ctypes.data,m.ctypes.data, bc.ctypes.data, npt, npix)

    t2 = time.time()
    if(benchmark):print("Grid snap calc. time: ",t2-t11)

#-------------------------------------------------------------------------------------------------------------
    # Potential Calculation
    t1 = time.time()

    density_ft = rfftn(grid)
    potential = irfftn(density_ft*green_ft)
    t2 = time.time()
    if(benchmark):print("FFT calc. time: ",t2-t1)
#-------------------------------------------------------------------------------------------------------------
    # Force calculation
    t1 = time.time()
    Fx,Fy,Fz = np.gradient(potential[:npix,:npix,:npix])

    t2 = time.time()
    if(benchmark):print("Gradient Calc. time: ",t2-t1)
#-------------------------------------------------------------------------------------------------------------

    # Integration of particle positions and velocities (LEAPFROG method)
    t1 = time.time()    

    # Update positions
    x = x + vx*dt + (1/2)*ax_prev*dt**2
    y = y + vy*dt + (1/2)*ay_prev*dt**2
    z = z + vz*dt + (1/2)*az_prev*dt**2

    # Save the position data
    if(t%save_iter == 0):
        x_save[t//save_iter] = x[::save_nth].copy()
        y_save[t//save_iter] = y[::save_nth].copy()
        z_save[t//save_iter] = z[::save_nth].copy()

    lx = np.asarray(np.round(x),dtype='int')
    ly = np.asarray(np.round(y),dtype='int')
    lz = np.asarray(np.round(z),dtype='int')
    
    handle_isolated_boundaries(lx.ctypes.data,ly.ctypes.data,lz.ctypes.data,npt,npix)
    l = np.array([lx,ly,lz])

    # Compute current acceleration
    ax = Fx[tuple(l)]/m
    ay = Fy[tuple(l)]/m
    az = Fz[tuple(l)]/m

    # Handle particles outside the grid >>> they feel no force from the grid bound mass.
    bc_flag = bc.astype(bool)

    ax[bc_flag] = 0
    ay[bc_flag] = 0
    az[bc_flag] = 0

    # Update velocities
    vx = vx + (1/2)*( ax_prev + ax )*dt
    vy = vy + (1/2)*( ay_prev + ay )*dt
    vz = vz + (1/2)*( az_prev + az )*dt
        
    # Save current forces    
    ax_prev = ax.copy()
    ay_prev = ay.copy()
    az_prev = az.copy()

    t2 = time.time()
    if(benchmark):print("Integration time: ",t2-t1)
    else: print("Progress: %.2f %%, Rate: %.3f s/iter, Elasped time: %.2f s \r"%( (t/niter)*100,(t2-t11),(t2-t91) ), end = '')
#-------------------------------------------------------------------------------------------------------------
    # Plotting 

    if(show3D and t%plt_iter == 0): 
        axs.clear()
        axs.set_axis_off()

        axs.plot3D(x_save[:t,1],y_save[:t,1],z_save[:t,1])
        axs.scatter3D(x_save[t,1],y_save[t,1],z_save[t,1])
        '''
        axs.scatter3D(x_save[t,2],y_save[t,2],z_save[t,2])
        axs.plot3D(x_save[:t,2],y_save[:t,2],z_save[:t,2])

        axs.scatter3D(x_save[t,3],y_save[t,3],z_save[t,3])
        axs.plot3D(x_save[:t,3],y_save[:t,3],z_save[:t,3])

        axs.scatter3D(x_save[t,4],y_save[t,4],z_save[t,4],color = "red")
        axs.plot3D(x_save[:t,4],y_save[:t,4],z_save[:t,4],color = "red")
        '''
        #axs.plot3D(x_save[:t,0],y_save[:t,0],z_save[:t,0])
        axs.scatter3D(x_save[t,0],y_save[t,0],z_save[t,0])

        axs.set_xlim([0,npix]);axs.set_ylim([0,npix]);axs.set_zlim([0,npix])
        plt.pause(0.001)
#-------------------------------------------------------------------------------------------------------------
    # Reset the grid
    t1 = time.time()
    grid.fill(0)

    t2 = time.time()
    if(benchmark):print("Grid reset time: ",t2-t1)
#-------------------------------------------------------------------------------------------------------------
t99 = time.time(); print("\n\nTotal simulation time: %.2f s"%(t99-t91))

# Save x,y and z data to files
if(save_data):
    if(verbose): print("\nSaving Data ...")
    np.savetxt(wd+"/Positional_data/"+x_file,x_save,delimiter=",")
    np.savetxt(wd+"/Positional_data/"+y_file,y_save,delimiter=",")
    np.savetxt(wd+"/Positional_data/"+z_file,z_save,delimiter=",")
    if(verbose): print("Save complete!")