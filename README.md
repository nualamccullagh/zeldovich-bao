# zeldovich-bao
The Zel'dovich approximation is a method used to set up the initial conditions of cosmological N-body simulations. It takes particles on a regular grid and moves them in a single time step so that the distribution of particles after they have been moved is close to Gaussian and has a specific power spectrum.

This code sets up an initial 3D grid of a specified number of particles in a cubic box of specified size, and applies the Zel'dovich approximation using the given power spectrum (supplied in a text file). The output is the final distribution of particles at a specified redshift. Code is also supplied to interpolate final particle distribution to a density grid through cloud-in-cell (CIC) and to compute the 2-point correlation function of the density field. The Zel'dovich approximation can be applied in real-space or redshift-space.

Files:

-zeldovich.py: This is the main code which runs the Zeldovich approximation for a given input power spectrum and output redshift and returns the final particle positions in either real or redshift space

-cic_dens_wrapper.py: This is a wrapper for the C-code (cicdensmodule.c) that computes the cloud-in-cell density on a grid from a set of particle positions

-spatial_stats.py: This module takes a density grid and computes the (1d) power spectrum, (1d and 2d) correlation function

-pk_indra7313.txt: An example linear power spectrum (generated using CAMB for WMAP7 cosmology)

-runZA.py: Script that shows how to run the above code and plot the resulting correlation function(s)

-dcicmod_script.sh: shell script that compiles cicdensmodule.c into a shared object library that can be called in python. Location of Python libraries may have to be changed for different systems
