#include <stdio.h>
//Compile into a shared library from the command line using the command:
//gcc-9 -o PM_methods.so C_methods.c -O3 --shared

// Assign particles to grid w/ isolated BC's
void to_grid(long *lx,long *ly,long *lz, double *grid, double *mass,long *bounds, long npt,long npix)
{
  for (long i=0;i<npt;i++) 
  {
    if(lx[i] >= 0 && lx[i] <= npix-1){
      if(ly[i] >= 0 && ly[i] <= npix-1){
        if(lz[i] >= 0 && lz[i] <= npix-1)
        { 
          long myind = 4*lx[i]*npix*npix + 2*ly[i]*npix + lz[i];
          grid[myind] += mass[i];
        }
        else
        {bounds[i] = 1;}
      }
      else
      {bounds[i] = 1;} 
    }
    else
    {bounds[i] = 1;}
  } 
}

// Assign particles to grid w/ periodic BC's
void to_grid_periodic(long *lx,long *ly,long *lz, double *grid, double *mass, long npt,long npix)
{
  for (long i=0;i<npt;i++) 
  {
    long myind = lx[i]*npix*npix + ly[i]*npix + lz[i];
    grid[myind] += mass[i];
  } 
}

// Handle out of bound particles for isolated BC
void handle_isolated_boundaries(long *lx,long *ly,long *lz,long npt,long npix)
{
    for (long i=0;i<npt;i++) 
    {
        if(lx[i] < 0    )  lx[i] = 0;  
        if(lx[i] > npix-1) lx[i] = npix-1;

        if(ly[i] < 0     ) ly[i] = 0;
        if(ly[i] > npix-1) ly[i] = npix-1;
        
        if(lz[i] < 0     ) lz[i] = 0;
        if(lz[i] > npix-1) lz[i] = npix-1;
    }
}

// Handle out of bound particles for periodic bc
void handle_periodic_bc(double *x, double *y, double *z, long npt, long npix)
{
    for(int i = 0; i<npt; i++)
    {
      // Check x bc
      if(x[i] > npix-1)
      {x[i] = 0;} 
      else if(x[i] < 0)
      {x[i] = npix-1;}

      // Check y bc
      if(y[i] > npix-1)
      {y[i] = 0;} 
      else if(y[i] < 0)
      {y[i] = npix-1;}

      // Check z bc
      if(z[i] > npix-1)
      {z[i] = 0;} 
      else if(z[i] < 0)
      {z[i] = npix-1;}
    }
}