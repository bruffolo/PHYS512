# N body simulation

The main N body code that can be used for large ranges of particle numbers are the programs:
- Nparticle_periodic.py
- Nparticle_isolated.py

These two programs are for handling the periodic and isolated systems, and are largely similar except for the green's
function calculation and handling of out of bounds particles.

Along with these programs is a C library "C_methods.c" that is called by the python code using the c_types library. The C library must be compiled with the command "gcc-9 -o PM_methods.so C_methods.c -O3 --shared" in the command line. This compiles to a shared library called "PM_methods.so".