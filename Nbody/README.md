# N body simulation (BUG IN THIS CODE, WILL FIX ASAP) <img src= ./photo.png width="400"> <img src= ./Orbit_complex.png width="350">

The main N body code that can be used for large ranges of particle numbers are the programs:
- Nparticle_periodic.py
- Nparticle_isolated.py

These two programs are for handling the periodic and isolated systems, and are largely similar except for the green's
function calculation and handling of out of bounds particles.

### Dependencies
Novel requirements are c_types and numba packages.

A C library "C_methods.c" is called by the above code using [c_types] in python. The C library must be compiled with the command 

```gcc-9 -o PM_methods.so C_methods.c -O3 --shared```

in the command line. This compiles to a shared library called "PM_methods.so". The particle grid assignment and handling at the boundaries are operations that involve looping over particle positions, and are natively very slow in python. These calculations are instead passed off the the compiled c library for speed augmentation. 

### Simulation Configuration
To run a simulation, open either the Nparticle_periodic.py or Nparticle_isolated.py, and scroll down to the system configuration. Here you can input the number of particles, as well as the grid size and number of iterations and time step of the integration loop.


```
npix = 2**6        # Pixel number
npt  = 300000      # Number of particles
ndim = 3           # Number of dimensions

# Number of iterations
niter = 10000

# Time step
dt = 0.001
```
The particle masses, and initial conditions can be set in the "particle setup" portion of the code.

### Potential
A softened 1/r potential is used, and the grid based potential map for a single particle (w/isolated BC's)  is shown below with comparison to the analytic 1/r potential.

<img src= ./Isolated_potential_map.png width="350">  <img src= ./Isolated_potential_plot.png width="350"> 

A potential map showing the proper implementation of periodic boundary conditions is shown below:
![](./Periodic_potential_map.png)


### Integrator

This simulator uses the standard leapfrog integration technique.

### Animations

The Animation.py program is written to animate (really just rapidly plot) the outputs of the isolated and periodic Nbody simulation from the above scripts that have been pre-computed and saved in textfile form in the "Positional_data" folder. With the folder contents downloaded (and directory structure intact), you can run

```python Animation.py isolated```
or 
```python Animation.py periodic```

to see an animation of an nbody simulation with 300,000 randomly distributed particles, time step (dt) of 1e-3 and grid size of 128X128X128. For better handling of the figure only every 300 particles of the 300,000 are plotted. The figure is 3 dimensional and fully interactive, allowing the user to rotate figure at will during the animation. 

The Orbit_Animator.py program does a similar job as the above animator, but with extra lines drawn behind the moving particles illustrating their path in space, this is very handy for the visualization of a circular orbit. You can run

```python Orbit_Animation.py circular```
or 
```python Orbit_Animation.py complex```

to animate a pre-computed 2 body circular orbit, or more complex solar system like orbit, respectively. Here the orbiting mass(es) are much lighter than the fixed mass at the centre, the grid size is 64X64X64 and the integrator time step is 0.001. 

## Modes

### Verbose Mode
Adding the switch ```-v``` as a command line argument when calling the above will print all the information associated with a given simulation. 

### Benchmark Mode
Adding the switch ```-b``` as a command line argument, the speed of each component of the algorithm for any particle number/grid size will be printed to stdout. You should also set ```niter = 1``` (unless you want to see the output for many iterations of the integrator).

For example, for 100,000 particles @ 128X128X128 grid size, the benchmark mode returns:

```
####### Benchmarks #######
Grid snap calc. time:  0.004080295562744141 s
FFT calc. time:  1.4266560077667236 s
Gradient Calc. time:  0.06095528602600098 s
Integration time:  0.02000594139099121 s
Grid reset time:  0.0131072998046875 s
```
Showing that the FFT is largely the time limiting component of the calculation.

### Live Plotting Mode
Adding the switch ```-p``` as a command line argument, you can have the code do live animation during the computation. Note that for live plotting to not significantly slow down the code, you could set ```nth = npt/1000```, and set ```plt_iter ``` as is necessary. These commands control the number of particles that are live plotted (every "nth" particle is plotted) and how many iterations of the integrator are performed before the plot is updated.

### Data Saving Mode
Adding the switch ```-s``` as a command line argument allows you to save the calculated positional data of the particles. The data will be saved into a subdirectory of the "Data" folder, the name of this directory is given by the ```sim_name``` variable.

## Questions

### 1)

### 2) 
The algorithm employed has good energy conservation properties and does well at maintaining a circular orbit. 
The stability of a circular orbit with this algorithm can be seen by running the ```python Orbit_Animation.py circular``` to animate pre-computed data of the motion of a light mass in orbit around a much heavier mass. If instead you want to see an animation during the live computation of the particles motion, run ```python Orbit.py```.

### 3) 
Run either ```Nparticle_periodic.py``` or ```Nparticle_isolated.py``` with the parameters of your choice for grid size, dt, and the number of particles. You can have the code do live animation during the computation by setting ```show3D = True``` in the code.  

The periodic boundary conditions conserve energy relatively well until dense clusters of mass form, at which point the system energy can deflect up to 10% away from the initial value.
During the initial collapse of mass, the isolated and periodic boundary conditions have approximately equal energy conservation characteristics. As mass begins to be ejected from the grid edges, the isolated boundary conditions conserve energy much more poorly, although this is largely expected because of the way particles are handled after leaving the grid range. Specifically, after leaving the grid, (for the isolated BC's) the particles no longer feel any force from the grid bound mass and their gravitational potential potential energy is no longer added into the energy calculation. Energy conservation is plotted in the 3D animations as a percentage of current/initial energy.


