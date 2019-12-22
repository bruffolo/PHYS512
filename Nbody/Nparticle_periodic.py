import numpy as np
from numpy.fft import rfftn,irfftn
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import ctypes
from numba import jit
import sys

# Enable command line switches
switches = True

# Pull necessary functions from companion C library
c_lib = ctypes.cdll.LoadLibrary("PM_methods.so")

# Get the methods for the periodic BC's
to_grid_periodic = c_lib.to_grid_periodic
to_grid_periodic.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_long, ctypes.c_long]

handle_periodic_bc = c_lib.handle_periodic_bc
handle_periodic_bc.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,ctypes.c_long, ctypes.c_long]


def to_grid_py(xy,grid,cell_loc):
    g_xy = np.asarray(np.round(xy),dtype='int')
    n = xy.shape[0]
    for i in range(n):
        grid[g_xy[i,0],g_xy[i,1],g_xy[i,2]] = grid[g_xy[i,0],g_xy[i,1],g_xy[i,2]] + 1
        cell_loc[i] = g_xy[i,0],g_xy[i,1],g_xy[i,2]

########################################
###      System Configuration        ###
########################################
npix = 2**8        # Pixel number
npt  = 50      # Number of particles
ndim = 3           # Number of dimenesions

# Particle mass parameters
m_equal  = False # Set all masses to common value 
m_random = True  # Set masses of particles randomly

m_max = 1        # Maximum particle mass 
m_min = 0.4      # Mininmum particle mass

# Number of iterations
niter = 3000

# Time step
dt = 0.001

########################################
###             Modes                ###
########################################

########### Plotting parameters (live during simulation) ###########
show3D = False
nth    = 300        # Plot every nth particle on screen
plt_iter = 5  # Plot particles every plt_iter iteration

########### Data saving parameters ###########
save_data = True
save_iter = 10   # Save position particle data every save_iter iteration
save_nth  = 30   # Save position particle data every save_nth particle
filename = "periodic_l.txt"

########### Benchmark mode ###########
benchmark = False
# Perferably if benchmark = True, you should set niter = 1, otherwise benchmarks will print out every iteration!

########### Potential test mode ###########
potential_test = False

# Switches from user input
args = len(sys.argv) - 1
cmd_pars = ["-b","-v","-p","-s","-pt"]

if(switches):
    for i in range(args):
        if(sys.argv[i+1] == cmd_pars[0]):
            benchmark = True
        if(sys.argv[i+1] == cmd_pars[1]):
            verbose = True
        if(sys.argv[i+1] == cmd_pars[2]):
            show3D = True
        if(sys.argv[i+1] == cmd_pars[3]):
            save_data = True
        if(sys.argv[i+1] == cmd_pars[4]):
            potential_test = True

########################################
###        Green's function          ###
########################################
npix = npix//2
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

npix = 2*npix
grid = np.zeros([npix,npix,npix]) # Make the grid array


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

# Starting positions
x = np.random.rand(npt)*(npix-0.6)
y = np.random.rand(npt)*(npix-0.6)
z = np.random.rand(npt)*(npix-0.6)

if(potential_test):
    niter = 0
    show3D = False
    z[:] = npix//2
    x[0] = npix//2
    y[0] = npix//2


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

# Assign particles to the grid
to_grid_periodic(lx.ctypes.data, ly.ctypes.data, lz.ctypes.data, grid.ctypes.data,m.ctypes.data, npt, npix)

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
    plt.imshow(potential[:npix,:npix,npix//2].T,origin='lower');plt.colorbar()
    plt.savefig("Periodic_potential_map.png")
    plt.show()

########################################
###            Main Loop             ###
########################################

if(show3D):
    fig = plt.figure(1)
    axs = plt.axes(projection='3d')
    axs.set_axis_off()
    axs.view_init(60, 0)

for t in range(niter):
    t11 = time.time()
    to_grid_periodic(lx.ctypes.data, ly.ctypes.data, lz.ctypes.data, grid.ctypes.data,m.ctypes.data, npt, npix)

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

    # Handle particles encountering the boundary
    handle_periodic_bc(x.ctypes.data, y.ctypes.data, z.ctypes.data, npt, npix)
    
    # Save the position data
    if(save_data and t%save_iter == 0):
        x_save[t//save_iter] = x[::save_nth].copy()
        y_save[t//save_iter] = y[::save_nth].copy()
        z_save[t//save_iter] = z[::save_nth].copy()

    lx = np.asarray(np.round(x),dtype='int')
    ly = np.asarray(np.round(y),dtype='int')
    lz = np.asarray(np.round(z),dtype='int')
    l = np.array([lx,ly,lz])

    ax = Fx[tuple(l)]/m
    ay = Fy[tuple(l)]/m
    az = Fz[tuple(l)]/m

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
    else:print("Progress: %.2f %%, Rate: %.3f s/iter, Elasped time: %.2f s \r"%((t/niter)*100,(t2-t11),(t2-t91)), end = '')

#-------------------------------------------------------------------------------------------------------------
    # Plotting 
    if(0):
        plt.clf()
        #plt.imshow(potential[:npix,:npix,npix//2].T,origin='lower');plt.colorbar()
        plt.scatter(x,y,s=1)
        plt.xlim([0,npix]);plt.ylim([0,npix]);plt.axis('off')
        plt.pause(0.001)

    if(show3D and t%plt_iter == 0):

        axs.clear()
        axs.set_axis_off()
        axs.scatter3D(x[::nth],y[::nth],z[::nth],s=m[::nth])
        axs.set_xlim([0,npix]);axs.set_ylim([0,npix]);axs.set_zlim([0,npix])
        #text = '\nBound Mass: %.2f %%\n'%( (np.sum(grid)/np.sum(m))*100 )
        #axs.text2D(-0.01, 0.99, text,fontsize = 8, transform=axs.transAxes)
        plt.pause(.1)
#-------------------------------------------------------------------------------------------------------------
    
    # Reset the grid
    t1 = time.time()
    grid.fill(0)

    t2 = time.time()
    if(benchmark):print("Grid reset time: ",t2-t1)
#-------------------------------------------------------------------------------------------------------------
t99 = time.time(); print("\n\nTotal simulation time: %.2f s"%(t99-t91))

x_file = "Positional_data/"+"x_"+filename
y_file = "Positional_data/"+"y_"+filename
z_file = "Positional_data/"+"z_"+filename

# Save x,y and z data to files
if(save_data):
    if(verbose): print("\nSaving Data ...")
    np.savetxt(x_file,x_save,delimiter=",")
    np.savetxt(y_file,y_save,delimiter=",")
    np.savetxt(z_file,z_save,delimiter=",")
    if(verbose): print("Save complete!")