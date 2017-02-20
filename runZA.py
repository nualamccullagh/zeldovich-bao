import zeldovich as Z
import cic_dens_wrapper
import numpy as N
import matplotlib.pyplot as plt
import spatial_stats


# This shows how to use the code to run a Zeldovich simulation with the power spectrum in the directory
# It takes a redshift as an input parameter, and an optional smoothing window parameter (default = 1 Mpc/h)
# It uses a boxsize of 1Gpc/h and 128^3 particles and redshift-space distortion parameter f=0.5
# computes the 1d and 2d correlation functions and plots both

def test(redshift, smw=1.0):
	pkinit=N.loadtxt('pk_indra7313.txt')
	boxsize=1000.0
	ngrid=128
	
	dens0, x, y, z=Z.run(redshift, pkinit, 0.5, boxsize=boxsize, ngrid=ngrid, nparticles=ngrid, exactpk=True, smw=smw)

	dens=cic_dens_wrapper.get_dens(x, y, z, ngrid, boxsize)
	r, xi=spatial_stats.getXi(dens, nrbins=ngrid/2, boxsize=boxsize, get2d=False, deconvolve_cic=True, exp_smooth=0.0)

	rp, pi, xi2d=spatial_stats.getXi(dens, nrbins=ngrid/2, boxsize=boxsize, get2d=True, deconvolve_cic=True, exp_smooth=0.0)

	plt.plot(r, r**2*xi)
	plt.xlim([0, 200])
	plt.figure()

	plt.pcolormesh(rp, pi, N.arcsinh(xi2d*300.0), cmap='rainbow')
	plt.xlim([0, boxsize/2])
	plt.ylim([0, boxsize/2])
    
    plt.xlabel(r'$\xi_1(r)$')
	plt.colorbar()

	plt.show()
