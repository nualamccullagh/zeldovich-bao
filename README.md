# zeldovich-bao
Code for running a Zel'dovich simulation, interpolating particles to a grid, and computing the 2pt-correlation function (in real and redshift space)

Files:

-zeldovich.py: This is the main code which runs the Zeldovich approximation for a given input power spectrum and output redshift and returns the final particle positions in either real or redshift space
-cic_dens_wrapper.py: This is a wrapper for the C-code (cicdensmodule.c) that computes the cloud-in-cell density on a grid from a set of particle positions
-spatial_stats.py: This module takes a density grid and computes the (1d) power spectrum, (1d and 2d) correlation function
-pk_indra7313.txt: An example linear power spectrum (generated using CAMB for WMAP7 cosmology)
-runZA.py: Script that shows how to run the above code and plot the resulting correlation function(s)
