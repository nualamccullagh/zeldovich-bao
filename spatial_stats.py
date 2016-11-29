import numpy as N

sqrt2 = N.sqrt(2.)

# density is the configuration-space density grid (real-valued, shape = ngrid x ngrid x ngrid)
# nkbins is the number of bins to compute the power spectrum in. It computes it in log-spaced bins in the range (2*Pi/L to Pi*ngrid / L)
# deconvolve_cic if you want to deconvolve a cloud-in-cell window function
# exp_smooth is the size of an exponential smoothing window to deconvolve (if = 0.0, then it does nothing)
def getPk(density, nkbins=40, boxsize=512., deconvolve_cic=False, exp_smooth=0.0):
    #make sure the density has mean 0
    density=density-N.mean(density)
    ngrid=density.shape[0]
    
    
    #Fourier transform of density
    deltak=N.fft.rfftn(density)
    sk=deltak.shape
    
    #Square the density in Fourier space to get the 3D power, make sure k=0 mode is 0
    dk2=(deltak*N.conjugate(deltak)).astype(N.float)
    dk2[0,0,0]=0.0
    
    #set up k-grid
    kmin=2*N.pi/boxsize
    kny=N.pi*ngrid/boxsize
    
    a = N.fromfunction(lambda x,y,z:x, sk).astype(N.float)
    a[N.where(a > ngrid/2)] -= ngrid
    b = N.fromfunction(lambda x,y,z:y, sk).astype(N.float)
    b[N.where(b > ngrid/2)] -= ngrid
    c = N.fromfunction(lambda x,y,z:z, sk).astype(N.float)
    c[N.where(c > ngrid/2)] -= ngrid
    kgrid = kmin*N.sqrt(a**2+b**2+c**2).astype(N.float)
    if (deconvolve_cic):
        wx=N.sin(kmin*a*N.pi/(2*kny))**4/(kmin*a*N.pi/(2*kny))**4
        wx[N.where(a==0)]=1.0
        wy=N.sin(kmin*b*N.pi/(2*kny))**4/(kmin*b*N.pi/(2*kny))**4
        wy[N.where(b==0)]=1.0
        wz=N.sin(kmin*c*N.pi/(2*kny))**4/(kmin*c*N.pi/(2*kny))**4
        wz[N.where(c==0)]=1.0
        ww=wx*wy*wz
        dk2=dk2/ww
        
    if (exp_smooth!=0):
        filt=N.exp(-kgrid**2*exp_smooth**2)
        dk2=dk2/filt
        
        
    #Now we want to compute the 1-D power spectrum which involves averaging over shells in k-space
    
    #define the k-bins we want to compute the power spectrum in
    binedges=N.logspace(N.log10(kmin), N.log10(kny),nkbins)
    numinbin=N.zeros_like(binedges)
    pk=N.zeros_like(binedges)
    kmean=N.zeros_like(binedges)
    
    kgrid=kgrid.flatten()
    dk2 = dk2.flatten()
    index = N.argsort(kgrid)
    
    kgrid=kgrid[index]
    dk2=dk2[index]
    c0=0.*c.flatten()+1.
    c0[N.where(c.flatten()==0.)]-=0.5
    c0=c0[index]
    cuts = N.searchsorted(kgrid,binedges)
    
    
    for i in N.arange(0, nkbins-1):
        if (cuts[i+1]>cuts[i]):
            numinbin[i]=N.sum(c0[cuts[i]:cuts[i+1]])
            pk[i]=N.sum(c0[cuts[i]:cuts[i+1]]*dk2[cuts[i]:cuts[i+1]])
            kmean[i]=N.sum(c0[cuts[i]:cuts[i+1]]*kgrid[cuts[i]:cuts[i+1]])
    wn0=N.where(numinbin>0.)
    pk=pk[wn0]
    kmean=kmean[wn0]
    numinbin=numinbin[wn0]
    
    pk/=numinbin
    kmean/=numinbin
    
    pk*= boxsize**3/ngrid**6
    
    return kmean, pk
    
def getXi(density, nrbins=40, boxsize=512., get2d=False, deconvolve_cic=False, exp_smooth=0.0):
    #make sure the density has mean 0
    density=density-N.mean(density)
    ngrid=density.shape[0]
    
    
    #Fourier transform of density
    deltak=N.fft.rfftn(density)
    sk=deltak.shape
    
    #Square the density in Fourier space to get the 3D power, make sure k=0 mode is 0
    dk2=(deltak*N.conjugate(deltak)).astype(N.float)
    dk2[0,0,0]=0.0
    
    #set up k-grid
    kmin=2*N.pi/boxsize
    kny=N.pi*ngrid/boxsize
    
    a = N.fromfunction(lambda x,y,z:x, sk).astype(N.float)
    a[N.where(a > ngrid/2)] -= ngrid
    b = N.fromfunction(lambda x,y,z:y, sk).astype(N.float)
    b[N.where(b > ngrid/2)] -= ngrid
    c = N.fromfunction(lambda x,y,z:z, sk).astype(N.float)
    c[N.where(c > ngrid/2)] -= ngrid
    kgrid = kmin*N.sqrt(a**2+b**2+c**2).astype(N.float)
    if (deconvolve_cic):
        wx=N.sin(kmin*a*N.pi/(2*kny))**4/(kmin*a*N.pi/(2*kny))**4
        wx[N.where(a==0)]=1.0
        wy=N.sin(kmin*b*N.pi/(2*kny))**4/(kmin*b*N.pi/(2*kny))**4
        wy[N.where(b==0)]=1.0
        wz=N.sin(kmin*c*N.pi/(2*kny))**4/(kmin*c*N.pi/(2*kny))**4
        wz[N.where(c==0)]=1.0
        ww=wx*wy*wz
        dk2=dk2/ww
        
    if (exp_smooth!=0):
        filt=N.exp(-kgrid**2*exp_smooth**2)
        dk2=dk2/filt
        
    #Fourier transform back to config space
    xi3d = N.fft.irfftn(dk2)/ngrid**3
    
    #define coordinates in config space
    
    rmin=boxsize/N.float(ngrid)
    a = N.fromfunction(lambda x,y,z:x, (ngrid, ngrid, ngrid))
    a[N.where(a > ngrid/2)] -= ngrid
    b = N.fromfunction(lambda x,y,z:y, (ngrid, ngrid, ngrid))
    b[N.where(b > ngrid/2)] -= ngrid
    c = N.fromfunction(lambda x,y,z:z, (ngrid, ngrid, ngrid))
    c[N.where(c > ngrid/2)] -= ngrid
    
    
    #compute 2d correlation function, xi(rp, pi)
    #this only works for nrbins=ngrid/2, so we will just ignore nrbins here
    if (get2d):
        rp=rmin*N.sqrt(b**2+c**2)
        c0 = 0.*c + 1.
        c0[N.where(c == 0.)] -= 0.5
        binedges=rmin*N.linspace(0.0, ngrid/2+1,ngrid/2+1)
        xips = N.zeros((ngrid/2, ngrid/2), dtype=N.float)
        #line of sight separation, rz
        rz=rmin*N.sqrt(a[:,0,0]**2)
        iz=N.argsort(rz)
        rz=rz[iz]
        rzmean=N.zeros_like(binedges)
        cutsz=N.searchsorted(rz, binedges)
        c0_ind=c0[:,0,0][iz]
        # loop over planes perpendicular to pi direction
        for i in N.arange(ngrid/2):
            rzmean[i]=N.sum(c0_ind[cutsz[i]:cutsz[i+1]]*rz[cutsz[i]:cutsz[i+1]])/N.sum(c0_ind[cutsz[i]:cutsz[i+1]])
            index = N.argsort(rp[i,:,:].flatten())
            rp_sort = rp[i,:,:].flatten()[index]
            xi3d_sort = xi3d[i,:,:].flatten()[index]
            c0_indexed = c0[i,:,:].flatten()[index]
            cuts = N.searchsorted(rp_sort,binedges)
            numinbin = N.zeros_like(binedges)
            for ic in N.arange(0,ngrid/2):
                if (cuts[ic+1] > cuts[ic]):
                    xips[i,ic] = N.sum(c0_indexed[cuts[ic]:cuts[ic+1]]*xi3d_sort[cuts[ic]:cuts[ic+1]])/N.sum(c0_indexed[cuts[ic]:cuts[ic+1]])
	    index=N.argsort(rp[0, :, :].flatten())
        rp_sort=rp[0,:,:].flatten()[index]
        rmean=N.zeros_like(binedges)
        c0_indexed = c0[0,:,:].flatten()[index]
        cuts = N.searchsorted(rp_sort,binedges)
        for i in N.arange(0, ngrid/2):
            rmean[i]=N.sum(c0_indexed[cuts[i]:cuts[i+1]]*rp_sort[cuts[i]:cuts[i+1]])/N.sum(c0_indexed[cuts[i]:cuts[i+1]])
        return rmean[0:-1], rzmean[0:-1], xips
        
        
    else:
        rgrid=rmin*N.sqrt(a**2+b**2+c**2)
        #define the k-bins we want to compute the power spectrum in
        binedges=rmin*N.linspace(0.0, ngrid/2+1,nrbins)
        numinbin=N.zeros_like(binedges)
        xir=N.zeros_like(binedges)
        rmean=N.zeros_like(binedges)
        
        rgrid=rgrid.flatten()
        xi3d = xi3d.flatten()
        index = N.argsort(rgrid)
        
        rgrid=rgrid[index]
        xi3d=xi3d[index]
        c0=0.*c.flatten()+1.
        c0[N.where(c.flatten()==0.)]-=0.5
        c0=c0[index]
        cuts = N.searchsorted(rgrid,binedges)
        
        
        for i in N.arange(0, nrbins-1):
            if (cuts[i+1]>cuts[i]):
                numinbin[i]=N.sum(c0[cuts[i]:cuts[i+1]])
                xir[i]=N.sum(c0[cuts[i]:cuts[i+1]]*xi3d[cuts[i]:cuts[i+1]])
                rmean[i]=N.sum(c0[cuts[i]:cuts[i+1]]*rgrid[cuts[i]:cuts[i+1]])
        wn0=N.where(numinbin>0.)
        xir=xir[wn0]
        rmean=rmean[wn0]
        numinbin=numinbin[wn0]
        
        xir/=numinbin
        rmean/=numinbin
        
        return rmean, xir
    
    
    
    
    
    