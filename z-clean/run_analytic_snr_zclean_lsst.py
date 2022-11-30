#!/usr/bin/env python
# coding: utf-8

import numpy as np
import basic
import cosmology
import misctools
import local


lmin, lmax = 1, 2048
L = np.linspace(0,lmax,lmax+1)
experiment = 's4'
#experiment = 'so'

k, pk0 = np.loadtxt(local.root + 'modelw/Pk/Pklin_new.dat',unpack=True)

# source redshift
zs = [local.zcmb,local.zcmb,local.zcmb]

# mass spectra
ckk, ckI, cII, nlII = np.loadtxt(local.root+'modelw/cl/mass.dat',unpack=True,usecols=(1,2,3,4))
fdel = ckI/(cII+nlII+1e-30)

zn = 50

zmin, zmax = 0.0001, 20
z, dz = basic.bispec.zpoints(zmin,zmax,zn)
chi = basic.cosmofuncs.dist_comoving(z,**local.cps)


# CIB weight
wcib = local.cib_weight(z,local.cps,local.nu)

# modified kernel
wdel = np.zeros((zn,lmax+1))
for l in range(lmax+1):
    if l<100: continue #low-ell cut
    wdel[:,l] = fdel[l]*wcib

# clkk
skk = basic.bispec.cl_flat(local.cpmodel,z,dz,zs[:2],lmax,k,pk0)
dkk = basic.bispec.cl_flat(local.cpmodel,z,dz,zs[:2],lmax,k,pk0,wdel=wdel)

Lmaxs = np.arange(100,2100,100)
oSNRt  = np.zeros(len(Lmaxs))
oSNRm  = np.zeros(len(Lmaxs))
dSNRt  = np.zeros(len(Lmaxs))
dSNRm  = np.zeros(len(Lmaxs))

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

np.savetxt(local.root+local.cpmodel+'/snr/snr_kkk_zclean_cib_'+experiment+'_zn'+str(zn)+'.dat',np.array((Lmaxs,oSNRt,dSNRt,oSNRm,dSNRm)).T)


