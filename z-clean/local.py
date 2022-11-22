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

