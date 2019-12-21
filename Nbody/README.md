# N body simulation (STILL A BUG IN THIS CODE, WILL FIX ASAP)

The main N body code that can be used for large ranges of particle numbers are the programs:
- Nparticle_periodic.py
- Nparticle_isolated.py

These two programs are for handling the periodic and isolated systems, and are largely similar except for the green's
function calculation and handling of out of bounds particles.

Along with these programs is a C library "C_methods.c" that is called by the python code using the c_types library. The C library must be compiled with the command "gcc-9 -o PM_methods.so C_methods.c -O3 --shared" in the command line. This compiles to a shared library called "PM_methods.so".

The Animation.py program is written to animate (really just rapidly plot) the outputs of the isolated and periodic Nbody simulation from the above scripts that have been pre-computed and saved in textfile form in the "Positional_data" folder. With the folder contents downloaded (and directory structure intact), you can run

```python Animation.py isolated```
or 
```python Animation.py periodic```

to see an animation of an nbody simulation with 100,000 particles, time step (dt) of 1e-3 and grid size of 64X64X64. For better handling of the figure only every 100 particles of the 100,000 are plotted. The figure is 3 dimensional and fully interactive, allowing the user to rotate figure at will during the animation. 