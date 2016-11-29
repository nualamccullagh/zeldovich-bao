import numpy as N
import scipy.interpolate as I


# redshift: float 
# pkinit: array of [k, pk] for the initial power spectrum
# boxsize: float (Mpc/h)
# ngrid: integer
def run(redshift, pkinit, f, boxsize=512.0, ngrid=512, nparticles=512, exactpk=True, smw=2.0):
    
    #make initial Gaussian density field
    seed=314159
    dens0=make_gauss_init(pkinit, boxsize=boxsize, ngrid=ngrid, seed=seed, exactpk=exactpk, smw=smw)
    #get displacements on the grid
    fx, fy, fz=get_disp(dens0, boxsize=boxsize, ngrid=ngrid)
    #get final positions of particles at a given redshift
    x, y, z=get_pos(fx, fy, fz*(1+f), redshift, boxsize=boxsize, ngrid=ngrid)
    
    return dens0, x, y, z
    
    
    
    
def make_gauss_init(pkinit, boxsize=512.0,ngrid=512,seed=314159,exactpk=True, smw=2.0):
    thirdim = ngrid/2+1
    kmin = 2*N.pi/N.float(boxsize)
    sk = (ngrid,ngrid,thirdim)
    
    
    a = N.fromfunction(lambda x,y,z:x, sk).astype(N.float)
    a[N.where(a > ngrid/2)] -= ngrid
    b = N.fromfunction(lambda x,y,z:y, sk).astype(N.float)
    b[N.where(b > ngrid/2)] -= ngrid
    c = N.fromfunction(lambda x,y,z:z, sk).astype(N.float)
    c[N.where(c > ngrid/2)] -= ngrid
    
    kgrid = kmin*N.sqrt(a**2+b**2+c**2).astype(N.float)
    
    rs = N.random.RandomState(seed)
    
    if (exactpk):
        dk = N.exp(2j*N.pi*rs.rand(ngrid*ngrid*(thirdim))).reshape((ngrid,ngrid,thirdim)).astype(N.complex64)
    else:
        dk = N.empty((ngrid,ngrid,thirdim),dtype=N.complex64)
        dk.real = rs.normal(size=ngrid*ngrid*(thirdim)).reshape((ngrid,ngrid,thirdim)).astype(N.float32)
        dk.imag = rs.normal(size=ngrid*ngrid*(thirdim)).reshape((ngrid,ngrid,thirdim)).astype(N.float32)
        dk /= N.sqrt(2.0)
    
    
    filt=N.exp(-kgrid*kgrid*smw*smw)
    
    pkinterp=I.interp1d(pkinit[:,0], pkinit[:,1])
    
    
    if (kgrid[0,0,0]==0):
        dk[0,0,0]=0
        wn0=N.where(kgrid!=0)
        
        dk[wn0] *= N.sqrt(filt[wn0]*pkinterp(kgrid[wn0]))*ngrid**3/boxsize**1.5
    else:
        dk *= N.sqrt(filt*pkinterp(kgrid.flatten())).reshape(sk)*ngrid**3/boxsize**1.5
    dk=nyquist(dk)
    dens = N.fft.irfftn(dk)
    return dens
    
    
def get_disp(dens, boxsize=512.0, ngrid=512):

    cell_len=N.float(boxsize)/N.float(ngrid)
    thirdim=ngrid/2+1
    kmin = 2*N.pi/N.float(boxsize)
    dens = N.fft.rfftn(dens)
    sk = (ngrid,ngrid,thirdim)
    
    a = N.fromfunction(lambda x,y,z:x, sk).astype(N.float)
    a[N.where(a > ngrid/2)] -= ngrid
    b = N.fromfunction(lambda x,y,z:y, sk).astype(N.float)
    b[N.where(b > ngrid/2)] -= ngrid
    c = N.fromfunction(lambda x,y,z:z, sk).astype(N.float)
    c[N.where(c > ngrid/2)] -= ngrid
    
    
    xp=N.zeros((ngrid, ngrid, ngrid/2+1), dtype=N.complex)
    yp=N.zeros((ngrid, ngrid, ngrid/2+1), dtype=N.complex)
    zp=N.zeros((ngrid, ngrid, ngrid/2+1), dtype=N.complex)
    
    kgrid = kmin*N.sqrt(a**2+b**2+c**2).astype(N.float)
    
    xp.real =-kmin*a*(dens.imag)/(kgrid*kgrid)
    xp.imag=kmin*a*(dens.real)/(kgrid*kgrid)
    xp[0,0,0]=0.
    yp.real =-kmin*b*(dens.imag)/(kgrid*kgrid)
    yp.imag=kmin*b*(dens.real)/(kgrid*kgrid)
    yp[0,0,0]=0.
    zp.real =-kmin*c*(dens.imag)/(kgrid*kgrid)
    zp.imag=kmin*c*(dens.real)/(kgrid*kgrid)
    zp[0,0,0]=0.
    
    a=0
    b=0
    c=0
    kgrid=0
    
    xp = nyquist(xp)
    yp = nyquist(yp)
    zp = nyquist(zp)
    
    xp=N.fft.irfftn(xp)
    yp=N.fft.irfftn(yp)
    zp=N.fft.irfftn(zp)
    
    return xp, yp, zp
    
    


def get_pos(fx, fy, fz, redshift, boxsize=512.0, ngrid=512):
    cell_len=N.float(boxsize)/N.float(ngrid)

    
    #setup particles on a uniform grid
    sk = (ngrid,ngrid,ngrid)
    a = N.fromfunction(lambda x,y,z:x+0.5, sk).astype(N.float)
    b = N.fromfunction(lambda x,y,z:y+0.5, sk).astype(N.float)
    c = N.fromfunction(lambda x,y,z:z+0.5, sk).astype(N.float)
    a=cell_len*a.flatten()
    b=cell_len*b.flatten()
    c=cell_len*c.flatten()
    
    #displacements, scaled by the growth function at the redshift we want
    d1=growthfunc(1./(1+redshift))/growthfunc(1.)
    x=fx*d1
    y=fy*d1
    z=fz*d1
    
    #assuming ngrid=nparticles, displace particles from the grid
    a+=x.flatten()
    b+=y.flatten()
    c+=z.flatten()
    
    #periodic boundary conditions
    a[N.where(a<0)]+=boxsize
    a[N.where(a>boxsize)]-=boxsize
    b[N.where(b<0)]+=boxsize
    b[N.where(b>boxsize)]-=boxsize
    c[N.where(c<0)]+=boxsize
    c[N.where(c>boxsize)]-=boxsize
    return a, b, c


def growthfunc(a, omega_m=0.289796, omega_l=0.710204):
    da=a/10000.
    integral = 0.
    for i in range(10000):
        aa=(i+1)*da
        integral+=da/((aa*N.sqrt(omega_m/(aa**3)+omega_l))**3)
    return 5*omega_m/2*N.sqrt(omega_m/a**3+omega_l)*integral

# ensuring arrays are hermitian
def nyquist(xp):
    ngrid=xp.shape[0]
    xp[ngrid/2+1:,1:,0]= N.conj(N.fliplr(N.flipud(xp[1:ngrid/2,1:,0])))
    xp[ngrid/2+1:,0,0] = N.conj(xp[ngrid/2-1:0:-1,0,0])
    xp[0,ngrid/2+1:,0] = N.conj(xp[0,ngrid/2-1:0:-1,0])
    xp[ngrid/2,ngrid/2+1:,0] = N.conj(xp[ngrid/2,ngrid/2-1:0:-1,0])
    
    xp[ngrid/2+1:,1:,ngrid/2]= N.conj(N.fliplr(N.flipud(xp[1:ngrid/2,1:,ngrid/2])))
    xp[ngrid/2+1:,0,ngrid/2] = N.conj(xp[ngrid/2-1:0:-1,0,ngrid/2])
    xp[0,ngrid/2+1:,ngrid/2] = N.conj(xp[0,ngrid/2-1:0:-1,ngrid/2])
    xp[ngrid/2,ngrid/2+1:,ngrid/2] = N.conj(xp[ngrid/2,ngrid/2-1:0:-1,ngrid/2])
    return xp