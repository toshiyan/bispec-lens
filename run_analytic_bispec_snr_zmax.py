#!/usr/bin/env python
# coding: utf-8

import numpy as np
import basic
import camb
import os
import cosmology
import misctools
from camb import model, initialpower
from camb.sources import GaussianSourceWindow, SplinedSourceWindow
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))


cpmodel = 'modelw'
H0    = 70.0
ombh2 = 0.02254
omch2 = 0.11417
As    = 2.41e-09
ns    = 0.97
h0    = H0/100.
Om    = (omch2+ombh2)/h0**2
cps   = {'H0':H0,'Om':Om,'Ov':1-Om,'w0':-1,'wa':0.}


lmin, lmax = 1, 2048
L = np.linspace(0,lmax,lmax+1)


# Set matter power spectrum at z: P(k,z)
fname = 'data/modelw/Pk/Pklin_new.dat'
if misctools.check_path(fname,verbose=True):
    k, pk0 = np.loadtxt(fname,unpack=True)
else:
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(lmax, lens_potential_accuracy=5)
    pars.NonLinear = model.NonLinear_both
    pars.Accuracy.AccuracyBoost = 1
    pars.Accuracy.lSampleBoost = 1
    pars.set_matter_power(redshifts=[0.],kmax=1e2/h0*1.1,k_per_logint=200)
    results = camb.get_results(pars)
    k, _, pk0 = results.get_linear_matter_power_spectrum(have_power_spectra=True,nonlinear=False)
    np.savetxt(fname,np.array((k,pk0[0])).T)

#k, pk0 = np.loadtxt( 'data/'+cpmodel+'/Pk/Pklin.dat', unpack=True )

zcmb = 1088.69

kk_noise = 'so'

for zn, zsi in [(10,.5),(20,1.),(20,1.5),(25,2.),(30,4.),(100,zcmb)]:
#for zn, zsi in [(20,1.),(25,2.),(30,4.),(100,zcmb)]:
    zs = [zsi,zsi,zsi]
    zmin, zmax = 0.0001, min(40.,zsi)
    z, dz = basic.bispec.zpoints(zmin,zmax,zn)
    chi = basic.cosmofuncs.dist_comoving(z,**cps)
    
    # clkk
    clkk = basic.bispec.cl_flat(cpmodel,z,dz,zs[:2],lmax,k,pk0)

    #Lmaxs = np.arange(100,300,100)
    Lmaxs = np.arange(100,2100,100)
    SNR = np.zeros(len(Lmaxs))

    for i, Lmax in enumerate(Lmaxs):

        # kappa noise spectrum
        nlkk = np.zeros(Lmax+1)
        if kk_noise == 'advact':
            nlkk[2:] = np.loadtxt('data/nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:Lmax-1]
        if kk_noise == 'so':
            nlkk[2:] = np.loadtxt('data/nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]

        # observed kappa spectrum
        ckk = clkk[:Lmax+1] + nlkk
        
        # SNR computation
        SNR[i] = basic.bispec.bispeclens_snr(cpmodel,'RT',z,dz,zs,2,Lmax,ckk,k,pk0,btype='kkk')
        print('zs=',zsi,',Lmax=',Lmax,',snr=',np.around(SNR[i],decimals=4))

    snr = np.array( [ np.sqrt(SNR[i]**2-SNR[i-1]**2) for i in range(1,len(Lmaxs)) ] )
    
    np.savetxt('data/'+cpmodel+'/snr/snr_kkk_zclean_'+str(zsi)+'_'+kk_noise+'.dat',np.array((Lmaxs,SNR)).T)

