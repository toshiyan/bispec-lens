#!/usr/bin/env python
# coding: utf-8

import numpy as np
import basic
import cosmology
import misctools
import local


# experimental parameters
lmin, lmax = 1, 2048
L = np.linspace(0,lmax,lmax+1)
experiment = 's4'
#experiment = 'so'
lss_zbin = 6

# load P(k) at z=0
k, pk0 = np.loadtxt(local.root + 'modelw/Pk/Pklin_new.dat',unpack=True)

# source redshift
zs = [local.zcmb,local.zcmb,local.zcmb]

# redshift integration points
zn = 50
zmin, zmax = 0.0001, 20
z, dz = basic.bispec.zpoints(zmin,zmax,zn)
chi = basic.cosmofuncs.dist_comoving(z,**local.cps)
Hzi = basic.cosmofuncs.hubble(z,divc=True,**local.cps)

# gal coefficient
frac = local.galaxy_distribution(np.linspace(0,50,1000),zbn={'lss':lss_zbin})[3]
Scov = local.get_covariance_signal(lmax,add_lss=lss_zbin)
slss = local.get_covariance_signal(lmax,add_cmb=False,add_cib=False,add_lss=lss_zbin)
nlss = local.get_covariance_noise(lmax,frac=frac,add_cmb=False,add_cib=False,add_lss=lss_zbin)
vecs = Scov[0,2:2+lss_zbin,:]
ocov = slss + nlss
icov = ocov.copy()
icov[:,:,1:] = np.array( [ np.linalg.inv(ocov[:,:,l]) for l in range(1,lmax+1) ] ).T
coef = np.array( [ np.dot( icov[:,:,l], vecs[:,l] ) for l in range(lmax+1) ] )

# gal weight
zbin, dndzi, pz, __ = local.galaxy_distribution(z,zbn={'lss':lss_zbin})
wgal = np.zeros((lss_zbin,zn))
for zi in range(lss_zbin):
    wgal[zi,:] = Hzi*dndzi['lss']*pz['lss'][zi]*(1+.84*z)/frac['lss'][zi]

# modified kernel
wdel = np.zeros((zn,lmax+1))
for l in range(lmax+1):
    wdel[:,l] = 0.
    for zi in range(lss_zbin):
        wdel[:,l] += coef[l,zi]*wgal[zi,:]

# clkk
skk = basic.bispec.cl_flat(local.cpmodel,z,dz,zs[:2],lmax,k,pk0)
dkk = basic.bispec.cl_flat(local.cpmodel,z,dz,zs[:2],lmax,k,pk0,wdel=wdel)

# prepare arrays
Lmaxs = np.arange(100,1600,100)
oSNRt  = np.zeros(len(Lmaxs))
oSNRm  = np.zeros(len(Lmaxs))
dSNRt  = np.zeros(len(Lmaxs))
dSNRm  = np.zeros(len(Lmaxs))

# compute SNR
for i, Lmax in enumerate(Lmaxs):

    # kappa noise spectrum
    nlkk = np.zeros(Lmax+1)
    if experiment == 'actdr6':
        nlkk[3:] = np.loadtxt(local.root+'nldd/ACTDR6.dat',unpack=True)[1][:Lmax-2]
    if experiment == 'so':
        nlkk[2:] = np.loadtxt(local.root+'nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]
    if experiment == 's4':
        nlkk[2:] = np.loadtxt(local.root+'nldd/kappa_deproj0_sens0_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]
        #nlkk[2:] = np.loadtxt(local.root+'nldd/S4_s1_t3_rlmax4000.dat',unpack=True)[1][:Lmax-1]

    # observed kappa spectrum
    ockk = skk[:Lmax+1] + nlkk
    dckk = dkk[:Lmax+1] + nlkk
        
    # SNR computation
    oSNRt[i], oSNRm[i] = basic.bispec.bispeclens_snr(local.cpmodel,'RT',z,dz,zs,2,Lmax,ockk,k,pk0,btype='kkk')
    dSNRt[i], dSNRm[i] = basic.bispec.bispeclens_snr(local.cpmodel,'RT',z,dz,zs,2,Lmax,dckk,k,pk0,btype='kkk',wdel=wdel[:,:Lmax+1])
    print('Lmax=',Lmax,',snr=',np.around(oSNRt[i],decimals=4),',snr(del)=',np.around(dSNRt[i],decimals=4))
    print('Lmax=',Lmax,',snr=',np.around(oSNRm[i],decimals=4),',snr(del)=',np.around(dSNRm[i],decimals=4))

np.savetxt(local.root+local.cpmodel+'/snr/snr_kkk_zclean_lss_'+experiment+'_zn'+str(zn)+'.dat',np.array((Lmaxs,oSNRt,dSNRt,oSNRm,dSNRm)).T)


