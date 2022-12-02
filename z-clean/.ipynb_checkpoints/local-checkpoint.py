#///////////////////////////////////////////////////////////////////////////////////////////////////#
# This file is intended to provide parametes, functions, etc, affecting the delensing code globally #
# Set up analysis parameters, filenames, arrays, functions                                          #
#///////////////////////////////////////////////////////////////////////////////////////////////////#

import numpy as np
import healpy as hp
import sys

# for camb
import os
import camb
from camb import model, initialpower
from camb.sources import GaussianSourceWindow, SplinedSourceWindow
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))

# from cmblensplus/wrap/
import basic

# from cmblensplus/utils/
import cosmology
import curvedsky as cs
import constant as c
import cmb


#////////// Define fixed values //////////#
# cosmological parameters
cpmodel = 'modelw'
root  = '../data/'
H0    = 70.0
h0    = H0/100.
Om    = 0.279
As    = 2.41e-09
ns    = 0.9645
ombh2 = 0.02254
omch2 = Om*h0**2 - ombh2
cps   = {'H0':H0,'Om':Om,'Ov':1-Om,'w0':-1,'wa':0.}
zcmb = 1088.69


# CIB model
nu = 353.


# CIB noise spectrum
def nlII(nu,lmax):
    l  = np.linspace(0,lmax,lmax+1)
    Jysr = c.MJysr2uK(nu)/c.Tcmb
    nI = 2.256e-10
    nl = ( nI + .00029989393 * (1./(l[:lmax+1]+1e-30))**(2.17) ) * Jysr**2
    nl[:100] = nl[100]
    return nl


def compute_matter_pk(H0,ombh2,omch2,As,ns):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(lmax, lens_potential_accuracy=5)
    pars.NonLinear = model.NonLinear_both
    pars.Accuracy.AccuracyBoost = 1
    pars.Accuracy.lSampleBoost = 1
    pars.set_matter_power(redshifts=[0.],kmax=1e2/(H0*.01)*1.1,k_per_logint=200)
    results = camb.get_results(pars)
    k, _, pk0 = results.get_linear_matter_power_spectrum(have_power_spectra=True,nonlinear=False)
    return k, pk0


def cib_weight(zi,cps,nu):
    rzi = basic.cosmofuncs.dist_comoving(zi,**cps)
    return cosmology.window_cib(rzi,zi,nu)


def kappa_noise(Lmax,experiment):
    
    # kappa noise spectrum
    nlkk = np.zeros(Lmax+1)
    if experiment == 'actdr6':
        nlkk[3:] = np.loadtxt(root+'nldd/ACTDR6.dat',unpack=True)[1][:Lmax-2]
    if experiment == 'so':
        nlkk[2:] = np.loadtxt(root+'nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]
    if experiment == 's4':
        nlkk[2:] = np.loadtxt(root+'nldd/kappa_deproj0_sens0_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:Lmax-1]

    return nlkk


# galaxy survey parameters
def galaxy_distribution( zi, survey=['lss'], zbn={'lss':4}, z0={'lss':.311}, nz_b={'lss':1.}, sig={'lss':.05}):
    
    zbin, dndzi, pz = {}, {}, {}
    
    if zbn['lss']==1:
        zbin['lss'] = np.array([2.5,7.])
    if zbn['lss']==2:
        zbin['lss'] = np.array([2.,2.5,7.])
    if zbn['lss']==3:
        zbin['lss'] = np.array([1.7,2.,2.5,7.])
    if zbn['lss']==4:
        zbin['lss'] = np.array([1.5,1.7,2.,2.5,7.])
    if zbn['lss']==5:
        zbin['lss'] = np.array([1.3,1.5,1.7,2.,2.5,7.])
    if zbn['lss']==6:
        zbin['lss'] = np.array([1.1,1.3,1.5,1.7,2.,2.5,7.])
        #zbin['lss'] = np.array([0.,.5,1.,2.,3.,4.,7.])
    if zbn['lss']==7:
        zbin['lss'] = np.array([0.9,1.1,1.3,1.5,1.7,2.,2.5,7.])

    for s in survey:
        dndzi[s] = basic.galaxy.dndz_sf(zi,2.,nz_b[s],z0=z0[s])
        pz[s]    = {zid: basic.galaxy.photoz_error(zi,[zbin[s][zid],zbin[s][zid+1]],sigma=sig[s],zbias=0.) for zid in range(zbn[s])}

    # fractional number density
    frac = {}
    for s in survey:
        frac[s] = {zid: np.sum(dndzi[s]*pz[s][zid])/np.sum(dndzi[s]) for zid in range(zbn[s]) }
    
    return zbin, dndzi, pz, frac


def tracer_list(add_cmb=True, add_lss=5, add_cib=True):
    
    # construct list of mass tracers to be combined
    klist = {}

    # store id for cmb lensing maps
    kid = 0
    
    if add_cmb:
        klist[kid] = 'cmb'
        kid += 1

    # store id for cib maps
    if add_cib: 
        klist[kid] = 'cib'
        kid += 1
        
    # store id for Euclid galaxy maps
    for z in range(add_lss):
        klist[kid] = 'lss'+str(z+1)+'n'+str(add_lss)
        kid += 1

    return klist

        
#//// Load analytic spectra and covariance ////#

def tracer_filename(m0,m1):

    return root + cpmodel + '/cl/'+m0+'-'+m1+'.dat'


def read_camb_cls(lmax=2048,lminI=100,return_klist=False,**kwargs):

    klist = tracer_list(**kwargs)
    
    # load cl of mass tracers
    cl = {}    
    for I, m0 in klist.items():
        for J, m1 in klist.items():
            if J<I: continue
            l, cl[m0+m1] = np.loadtxt( tracer_filename(m0,m1) )[:,:lmax+1]

            # remove low-ell CIB
            if m0=='cib' or m1=='cib':
                cl[m0+m1][:lminI] = 1e-20

    if return_klist:
        return l, cl, klist
    else:
        return l, cl
        

def get_covariance_signal(lmax,lmin=1,lminI=100,**kwargs): 
        # signal covariance matrix

        # read camb cls
        l, camb_cls, klist = read_camb_cls(lminI=lminI,return_klist=True,**kwargs)
        nkap = len(klist.keys())

        # form covariance
        Cov = np.zeros((nkap,nkap,lmax+1))
        
        for I, m0 in klist.items():
            for J, m1 in klist.items():
                if J<I: continue
                Cov[I,J,lmin:] = camb_cls[m0+m1][lmin:lmax+1]
                
        # symmetrize
        Cov = np.array( [ Cov[:,:,l] + Cov[:,:,l].T - np.diag(Cov[:,:,l].diagonal()) for l in range(lmax+1) ] ).T
        
        return Cov


def get_spectrum_noise(lmax,lminI=100,nu=353.,return_klist=False,frac=None,**kwargs):
    
    klist = tracer_list(**kwargs)

    l  = np.linspace(0,lmax,lmax+1)    
    nl = {}
    
    #//// prepare reconstruction noise of LB and S4 ////#
    if 'cib' in klist.values():
        Jysr = c.MJysr2uK(nu)/c.Tcmb
        nI = 2.256e-10
        nl['cib'] = ( nI + .00029989393 * (1./(l[:lmax+1]+1e-30))**(2.17) ) * Jysr**2
        nl['cib'][:lminI] = nl['cib'][lminI]

    for m in klist.values():
        if 'lss' in m:
            if frac is None:
                f = 1./kwargs['add_lss']
            else:
                f = frac['lss'][int(m[3])-1]
            nl[m] = np.ones(lmax+1)*c.ac2rad**2/(40.*f)

    for m in nl.keys():
        nl[m][0] = 0.
    
    if return_klist:
        return nl, klist
    else:
        return nl


def get_covariance_noise(lmax,lminI=100,frac=None,**kwargs):
    
    nl, klist = get_spectrum_noise(lmax,lminI=lminI,return_klist=True,frac=frac,**kwargs)
    nkap = len(klist.keys())

    Ncov = np.zeros((nkap,nkap,lmax+1))

    for I, m in enumerate(nl.keys()):
        Ncov[I,I,:] = nl[m]
 
    return Ncov

