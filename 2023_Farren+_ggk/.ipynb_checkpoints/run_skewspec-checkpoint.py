import basic, numpy as np, healpy as hp, os
from scipy.interpolate import InterpolatedUnivariateSpline as spline

D = '../data/2023_Farren+/'
cpmodel = 'modelw'
model   = 'RT'
gal     = 'green'

k, pk0 = np.loadtxt(D+'../'+cpmodel+'/Pk/Pklin.dat',unpack=True)
kn = np.size(k)

zcmb = 1088.69
zs   = [zcmb,zcmb]
zmin, zmax, zn = 0.01, 5, 30

# gal lrange
lmin, lmax = 20, 3000
# kappa lrange
Lmin, Lmax = lmin, 2048

Ln = Lmax-Lmin+1
L = np.linspace(Lmin,Lmax,Ln,dtype=np.int)

z, dz = basic.bispec.zpoints(zmin,zmax,zn)
zi, dndzi = np.loadtxt(D+'dNdz/'+gal+'.txt',unpack=True)
dNdz = spline(zi,dndzi)(z)
imax = np.argmax(z>np.max(zi))
dNdz[imax:] = 0

f = D+'Sl_'+gal+'_'+model+'_zn'+str(zn)+'_zmax'+str(zmax)+'_l'+str(lmin)+'-'+str(lmax)+'.dat'
skew = basic.bispec.skewspeclens(cpmodel,model,z,dz,zs,L,lmin,lmax,k,pk0,btype='kgg',dNdz=dNdz,theta=.1)
np.savetxt(f,np.array((L,skew[0,0,:],skew[0,1,:])).T)
