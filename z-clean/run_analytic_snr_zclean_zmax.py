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

for zmax in [.5,1.,1.5,2.,2.5]:

    zmin = 0.0001
    zn = int(20*(zmax-zmin))

    z, dz = basic.bispec.zpoints(zmin,zmax,zn)
    chi = basic.cosmofuncs.dist_comoving(z,**local.cps)

    # clkk
    skk = basic.bispec.cl_flat(local.cpmodel,z,dz,zs[:2],lmax,k,pk0)

    Lmaxs = np.arange(100,1600,100)
    SNRt  = np.zeros(len(Lmaxs))
    SNRm  = np.zeros(len(Lmaxs))

    for i, Lmax in enumerate(Lmaxs):
        
        nlkk = local.kappa_noise(Lmax,experiment)

        # observed kappa spectrum
        ockk = skk[:Lmax+1] + nlkk
        
        # SNR computation
        SNRt[i], SNRm[i] = basic.bispec.bispeclens_snr(local.cpmodel,'RT',z,dz,zs,2,Lmax,ockk,k,pk0,btype='kkk')
        print('Lmax=',Lmax,',snr=',np.around(SNRt[i],decimals=4),',snr(lss)=',np.around(SNRm[i],decimals=4))

    np.savetxt(local.root+local.cpmodel+'/snr/snr_kkk_zclean_'+experiment+'_zn'+str(zn)+'_zmax'+str(zmax)+'.dat',np.array((Lmaxs,SNRt,SNRm)).T)


