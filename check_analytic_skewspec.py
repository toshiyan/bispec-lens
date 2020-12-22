#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np, camb, basic
from matplotlib.pyplot import *


# In[2]:


# Function to get the sigma_8 and P(k) at z = 0.0 from CAMB
def get_pklin0(h0,wb,Om,As,ns,maxkh=10.,k_per_logint=200):

    # Reset CAMB
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=1e2*h0, ombh2=wb, omch2=Om*h0**2-wb,omk=0,mnu=0)
    pars.set_dark_energy() #LCDM (default)
    pars.InitPower.set_params(ns=ns, r=0, As=As)

    # Set matter power spectrum at z=0: P(k,z)
    pars.set_matter_power(redshifts=[0.],kmax=maxkh/h0*1.1,k_per_logint=k_per_logint)
    pars.NonLinear = camb.model.NonLinear_none
    
    # Calculate the intermediate results for these parameters
    results = camb.get_results(pars)
    #results.calc_power_spectra(pars)
    
    print(results.get_sigma8())
    
    # Calculate the CAMB power spectra to be interpolated
    kh, _z, pk = results.get_linear_matter_power_spectrum(have_power_spectra=True,nonlinear=False)

    return kh, pk[0]


# In[3]:


h0 = 0.7
wb = 0.022
Om = 0.3
As = 2e-9
ns = 0.96


# In[4]:


kh, pklin0 = get_pklin0(h0,wb,Om,As,ns)


# In[5]:


#loglog(kh,pklin0)


# In[6]:


lmin = 100
lmax = 2000
bn = lmax-lmin
ols = np.linspace(lmin,lmax,bn+1,dtype=np.int)
zss = [1.0334]
zn = 30
D = '/global/homes/t/toshiyan/scratch/bispec/skewspec/'
#zss = [0.5078]

for zs0 in zss:
    zs = [zs0,zs0]
    zmin, zmax = 0.0001, min(zs0,zs0)
    f = D+'Sl_Om'+str(Om)+'_h0'+str(h0)+'_zs'+str(zs0)+'_zn'+str(zn)+'_l'+str(lmin)+'-'+str(lmax)+'.dat'
    print(zs0)
    skew = basic.bispec.skewspeclens('input','RT',zmin,zmax,zn,zs,ols,lmin,lmax,kh,pklin0,pb=True,Om=Om,H0=1e2*h0,mnu=0.,ns=ns)
    np.savetxt(f,np.array((ols,skew[0,:],skew[1,:],skew[2,:])).T)

