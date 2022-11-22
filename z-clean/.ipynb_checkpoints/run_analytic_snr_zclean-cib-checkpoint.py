#!/usr/bin/env python
# coding: utf-8

import numpy as np
import basic
import camb
import os
import cosmology
import misctools


lmin, lmax = 1, 2048
L = np.linspace(0,lmax,lmax+1)

kk_noise = 's4'

# Set matter power spectrum at z: P(k,z)
fname = local.root + 'modelw/Pk/Pklin_new.dat'
if misctools.check_path(fname,verbose=True):
    k, pk0 = np.loadtxt(fname,unpack=True)
else:
    k, pk0 = local.compute_matter_pk(local.H0,local.ombh2,local.omch2,local.As,local.ns)
    np.savetxt(fname,np.array((k,pk0[0])).T)

# source redshift
zs = [zcmb,zcmb,zcmb]

# mass spectra
ckI, cII, nlII = np.loadtxt(local.root+'modelw/cl/mass.dat',unpack=True,usecols=(1,2))
fdel = ckI/(cII+nlII+1e-30)

#for zn, zimax in [(10,.5),(20,1.),(20,1.5),(25,2.)]:
for zn, zimax in [(10,.5)]:
#for zn, zimax in [(25,2.5),(30,3.),(30,4.),(100,zcmb)]:

    zmin, zmax = 0.0001, min(40.,zimax)
    z, dz = basic.bispec.zpoints(zmin,zmax,zn)
    chi = basic.cosmofuncs.dist_comoving(z,**cps)
    
    # CIB weight
    wcib = local.cib_weight(z,local.cps,local.nu)
    
    # modified kernel
    wdel = np.dot((fdel,wcib))
    #wdel = np.array([ fdel[zi]*wcib[lmax] ])

    # clkk
    clkk = basic.bispec.cl_flat(cpmodel,z,dz,zs[:2],lmax,k,pk0)

    #Lmaxs = np.arange(100,300,100)
    Lmaxs = np.arange(100,2100,100)
    SNR = np.zeros(len(Lmaxs))

    for i, Lmax in enumerate(Lmaxs):

        # kappa noise spectrum
        nlkk = np.zeros(Lmax+1)
        if kk_noise == 'advact':
            nlkk[2:] = np.loadtxt(local.root+'nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:Lmax-1]
        if kk_noise == 'so':
            nlkk[2:] = np.loadtxt(local.root+'nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]
        if kk_noise == 's4':
            nlkk[2:] = np.loadtxt(local.root+'nldd/S4_s1_t3_rlmax4000.dat',unpack=True)[1][:Lmax-1]

        # observed kappa spectrum
        ckk = clkk[:Lmax+1] + nlkk
        
        # SNR computation
        SNR[i] = basic.bispec.bispeclens_snr(cpmodel,'RT',z,dz,zs,2,Lmax,ckk,k,pk0,btype='kkk')
        print('zs=',zsi,',Lmax=',Lmax,',snr=',np.around(SNR[i],decimals=4))

    snr = np.array( [ np.sqrt(SNR[i]**2-SNR[i-1]**2) for i in range(1,len(Lmaxs)) ] )
    
    np.savetxt(root+cpmodel+'/snr/snr_kkk_zclean_'+str(zsi)+'_'+kk_noise+'.dat',np.array((Lmaxs,SNR)).T)

