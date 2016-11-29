import numpy as N
import scipy.interpolate as I
import scipy.special as S

#Returns the linear theory prediction for the real-space correlation function (z=0)
def xiLin(r, pklin, smw=2.0):
    return xinm(0, 0, r, pklin, smw=smw)



#Returns the linear theory prediction for the line-of-sight correlation function (z=0)
#the growth rate parameter (f) can be either a float or an array
#if a float, returns an array of the same length as r
#if an array, returns an array of shape (f.size, r.size)
def xi_LOS(r, f, pklin, smw=2.0):
    xi00=xinm(0, 0, r, pklin, smw=smw)
    xi20=xinm(2, 0, r, pklin, smw=smw)
    xi40=xinm(4, 0, r, pklin, smw=smw)


    if (N.asarray(f).size == 1):
        xiLOS=(1.0+2.0*f/3.0+f**2/5.0)*xi00 - (4.*f/3.+4.*f**2/7.)*xi20 + 8.*f**2/35.*xi40

    else:
        f=N.asarray(f)
        xiLOS=N.outer((1.0+2.0*f/3.0+f**2/5.0), xi00) - N.outer((4.*f/3.+4.*f**2/7.), xi20) + N.outer(8.*f**2/35., xi40)

    return xiLOS



def besph(n, x):
    z=N.array(x)
    if z.size == 1:
        if (n==0 & x==0.0):
            return 1.0
        elif x==0.0:
            return 0.0
    else:
        return S.jn(n+0.5, x)*N.sqrt(N.pi/(2*x))
    ans=N.zeros_like(z)
    nz=N.where(z!=0)
    if n==0:
        ans[N.where(z==0)]=1.0
    else:
        ans[N.where(z==0)]=0.0
    ans[nz]=S.jn(n+0.5, z[nz])*N.sqrt(N.pi/(2*z[nz]))
    return ans


def xinm(n, m, r, pklin, smw=2.0, L=20000, upper=4*N.pi):


    pkI=I.interp1d(pklin[:,0], pklin[:,1])

    deltak=upper/L
    kk=N.arange(pklin[:,0].min(), deltak*L, deltak)
    
    sk=r.shape
    
    jrk=besph(n, N.outer(r,kk))
    if (n==0):
        jrk[N.where(r==0)]=1.0
    else:
        jrk[N.where(r==0)]=0.0
    
    if (smw==0.0):
        w=1.0
    else:
        w=N.exp(-kk**2*smw**2/2.0)
    
    ans=N.sum((pkI(kk)*w**2*kk**(2+m)*jrk*deltak)/(2*N.pi**2), 1)

    ans=N.reshape(ans, sk)
    
    return ans
