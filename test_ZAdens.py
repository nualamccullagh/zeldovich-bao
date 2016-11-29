import numpy as N
import zeldovich as Z
import matplotlib.pyplot as plt
import cicdens

pkfile=''

def get_za_dens(redshift, boxsize=512.0, ngrid=512):
    pkinit=N.loadtxt(pkfile)
    x, y, z=Z.run(redshift, pkinit, boxsize=boxsize, ngrid=ngrid, nparticles=ngrid, exactpk=True)
    
    dens=N.zeros([ngrid, ngrid, ngrid], dtype=N.float)
    
    #compute the cloud-in-cell density on the grid "dens" using x, y, and z. Pointers are passed using x.ctypes.data and must be of type int64
    #particle positions and density grid must be of type numpy float (which are 64 bit, or double-precision floats)
    #the box size must also be a double precision float. ngrid and nparticles must be integers 
    cicdens.calcd_cic(N.int(ngrid), N.int(ngrid**3), N.float(boxsize), N.int64(x.ctypes.data), N.int64(y.ctypes.data), N.int64(z.ctypes.data), N.int64(dens.ctypes.data))
    
    #dens now contains the total mass (number) of particles in each cell, but we want the overdensity, which is:
    dens=dens/N.mean(dens)-1.0
    
    return dens
    
    
    
    