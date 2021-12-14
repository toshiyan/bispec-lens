#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import basic
import camb
import os
import cosmology
from matplotlib.pyplot import *
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


lmin = 1
lmax = 2048
L = np.linspace(0,lmax,lmax+1)


# Set matter power spectrum at z: P(k,z)

pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
pars.InitPower.set_params(As=As, ns=ns)
pars.set_for_lmax(lmax, lens_potential_accuracy=5)
pars.NonLinear = model.NonLinear_both
pars.Accuracy.AccuracyBoost = 1
pars.Accuracy.lSampleBoost = 1


pars.set_matter_power(redshifts=[0.],kmax=5e1/h0*1.1,k_per_logint=200)
results = camb.get_results(pars)
k, _, pk0 = results.get_linear_matter_power_spectrum(have_power_spectra=True,nonlinear=False)


# In[7]:
#k, pk0 = np.loadtxt( 'data/'+cpmodel+'/Pk/Pklin.dat', unpack=True )


zcmb = 1088.69


for zn, zsi in [(10,.5),(20,1.),(25,2.),(30,4.),(100,zcmb)]:
    zs = [zsi,zsi,zsi]
    zmin, zmax = 0.0001, min(40.,zsi)
    z, dz = basic.bispec.zpoints(zmin,zmax,zn)
    chi = basic.cosmofuncs.dist_comoving(z,**cps)
    
    # clkk
    clkk = basic.bispec.cl_flat(cpmodel,z,dz,zs[:2],lmax,k,pk0)

    #Lmaxs = np.array([100,500,1000,2000,3000])
    Lmaxs = np.array([100,200,300,400,500])
    snr = np.zeros(len(Lmaxs))

    for i, Lmax in enumerate(Lmaxs):
    
        nlkk = np.zeros(Lmax+1)
        #nlkk[2:] = np.loadtxt(D+'nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:lmax-1]
        nlkk[2:] = np.loadtxt('data/nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]
        ckk = clkk[:Lmax+1] + nlkk
        
        snr[i] = basic.bispec.bispeclens_snr(cpmodel,'RT',z,dz,zs,2,Lmax,ckk,k,pk0,btype='kkk')
        print('zs=',zsi,'Lmax=',Lmax,'snr=',np.around(snr[i],decimals=4))
    np.savetxt('data/'+cpmodel+'/snr/snr_kkk_zclean_'+str(zsi)+'_kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',np.array((Lmaxs,snr)).T)




