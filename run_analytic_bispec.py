
import basic
import numpy as np
import healpy as hp

cpmodel = 'modelp'
#model = 'RT'
model = 'GM'

#zmin, zmax = 0.0001, 1088.69
zmin, zmax = 0.0001, 1.0
zn = 100
zs = [zmax,zmax,zmax]

#calc = 'bispecsnr'
calc = 'bispecbin'

lmin = 1
lmax = 2000
olmin = lmin
olmax = lmax

bn = 20
lan = 0
kan = 0

D = 'data/'
L = np.linspace(0,olmax,olmax+1)
fpk = D+cpmodel+'/Pk/Pklin.dat'
k, pk0 = np.loadtxt(fpk,unpack=True)
kn = np.size(k)
z, dz = basic.bispec.zpoints(zmin,zmax,zn)

# load aps
nldd = np.zeros(lmax+1)
cldd = np.zeros(lmax+1)
nldd[2:] = np.loadtxt(D+'nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:lmax-1]
cldd[2:] = np.loadtxt(D+cpmodel+'/cl/fid.dat',unpack=True)[4][:lmax-1]
cldd[2:] = cldd[2:]/4.*(1.+1./L[2:])**2 + nldd[2:]*L[2:]**2*(L[2:]+1.)**2/4.

# bispectrum
if calc == 'bispec':
    bl0, pb0 = basic.bispec.bispeclens('equi',cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan)
    bl1, pb1 = basic.bispec.bispeclens('fold',cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan)
    bl2, pb2 = basic.bispec.bispeclens('sque',cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan)
    bl3, pb3 = basic.bispec.bispeclens('angl',cpmodel,model,z,dz,zn,zs,lmin,lmax,k,pk0,kn,lan,kan)
    np.savetxt('test_zs'+str(zs)+'.dat',np.array((L[1:],bl0,bl1,bl2,bl3,pb0,pb1,pb2,pb3)).T)

# binned bispectrum
if calc == 'bispecbin':
    shap = 'sque'
    bc, bl0, bl1 = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zn,zs,lmin,lmax,bn,k,pk0,kn,lan,kan)
    np.savetxt('bl_'+shap+'.dat',np.array((bc,bl0,bl1)).T)

# total bispectrum SNR
if calc == 'bispecsnr':
    snr = basic.bispec.bispeclens_snr(cpmodel,model,z,dz,zn,zs,2,lmax,cldd,k,pk0,kn)
    print(snr)
    #np.savetxt('bl_'+shap+'.dat',np.array((bc,bl0,bl1)).T)


