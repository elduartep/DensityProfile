# DensityProfile
Density profile computation

copile typing g++ Perfiles_2019_octantes_omp.c -lm -fopenmp

# You need to install
openmp

# IMPUTS
All files need to be in the working directory

## Dark matter catalog
ASCII file with x, y, z in Mpc/h

## Voids catalog
ASCII file with x, y, z, r in Mpc/h

## Halos catalog
ASCII file with x, y, z, m, r with xyzr in Mpc/h, m=mass in 10^10 M_\sun

## or Haloes catalog

# OUTPUTS




# parameters
parameters.h

you need so set here the cosmological parameters, the size of the box, number of particles, the name of the files, and if you want to work with halos or voids.
