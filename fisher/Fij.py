#//////////////////////////////////////////////////#
# * Fisher forecast combining 2pt, 3pt and DESI
#//////////////////////////////////////////////////#

import numpy as np, h5py
from matplotlib.pyplot import *

# read cl
def READ_CL(f,ln,l): 
  CL  = (np.loadtxt(f)).T
  fac = (2.*np.pi)/(l**2+l)
  cl  = np.zeros((3,3,ln))
  cl[0,0,:] = CL[1,:ln]*fac
  cl[0,1,:] = CL[3,:ln]*fac
  cl[0,2,:] = CL[5,:ln]/l**3
  cl[1,1,:] = CL[2,:ln]*fac
  cl[2,2,:] = CL[4,:ln]/l**4
  return cl

# return fisher matrix
def FisherMatrix(dlnCdp,fsky):
  s1, s2, ln, pn = dlnCdp.shape
  F = np.zeros((pn,pn,ln))
  for i in range(pn):
    for j in range(i,pn):
      for l in range(ln):
        F[i,j,l] = fsky*(l+2.+.5)*np.trace(np.dot(dlnCdp[:,:,l,i],dlnCdp[:,:,l,j]))
        F[j,i,l] = F[i,j,l]
  return F

# loading Fisher matrix of lensing bispectrum
def loadbispFij(fr,cp,dp,f):
  cn = len(cp)
  F2 = np.zeros((cn,cn))
  with h5py.File(f,'r') as hf:
    for p in range(cn):
      for q in range(p,cn):
        if cp[p]=='tau' or cp[q]=='tau':
          F2[p,q] = 0.
        else:
          if p==q:
            tag = cp[p]+'_dp'+str(fr[cp[p]][dp[p]])
          else:
            tag = cp[p]+'_dp'+str(fr[cp[p]][dp[p]])+'_'+cp[q]+'_dp'+str(fr[cp[q]][dp[q]])
          F2[p,q] = np.array(hf.get(tag))
  return F2

# values of fractional delta p
fr = {}
fr['Obh'] = np.array([.1])
fr['Omh'] = np.array([.007,.005,.003])
fr['O_L'] = np.array([.1,.07,.05])
fr['A_s'] = np.array([.1])
fr['n_s'] = np.array([.005])
#fr['w_0'] = np.array([.05,.03])
#fr['mnu'] = np.array([.3,.2,.1])
fr['tau'] = np.array([.4])


lmin=2; lmax=4000; ln=lmax-lmin+1; l=np.linspace(lmin,lmax,ln)
llmax = 3000 # phi
blmax = 1000 # bispectrum phi
rlmax = lmax # CMB
cp = ['w_0','mnu','Omh','Obh','O_L','n_s','A_s','tau']
cn = len(cp)
root = '../../../data/forecast/model0/'

# signal covariance 
CL = READ_CL(root+'cl/fid.dat',ln,l)

# experiments: fsky and noise covariance 
ac2rad = np.pi/180./60.
fsky = .5
nlpp = (1.*ac2rad/2.72e6)**2*np.exp(-l*(l+1)*(3.*ac2rad)**2/(8.*np.log(2.)))
nldd = (np.loadtxt(root+'../nldd/S4_s1_t3_rlmax'+str(rlmax)+'.dat')).T[1]
#nldd = (np.loadtxt(root+'../nldd/S5_rlmax'+str(rlmax)+'.dat')).T[1]

# noise covariance
NL = np.zeros((3,3,ln))
NL[0,0,:] = nlpp/2.
NL[1,1,:] = nlpp
NL[2,2,:] = nldd[:ln]
NL[2,2,llmax:ln] = NL[2,2,llmax:ln]*1000000
NL[0,0,:49] = NL[0,0,:49]*1000000
NL[1,1,:49] = NL[1,1,:49]*1000000

# DESI
sig = np.array([4.1,1.7,.88,.55,.38,.28,.21,.18,.18,.17,.16,.14,.15,.16,.19,.28,.41,.52])*1e-3
fp = np.zeros((cn,18))
#fp[3,:], fp[2,:], fp[4,:], fp[0,:] = np.loadtxt(root+'rsdv.dat',unpack=True,usecols=(1,2,3,4))
Fd = np.zeros((cn,cn))
for i, p in enumerate(cp):
  for j, q in enumerate(cp):
    Fd[i,j] = np.sum(fp[i,:]*fp[j,:]/sig**2)

# primary CMB + 2pt
cr = {'Obh':.1,'Omh':.005,'O_L':.03,'A_s':.1,'n_s':.005,'w_0':.005,'mnu':.2,'tau':.4}
iC = np.array([np.linalg.inv(CL[:,:,L]+NL[:,:,L]) for L in range(ln)])
dlnCdp = np.zeros((3,3,ln,cn))
for i, p in enumerate(cp):
  dC = READ_CL(root+'cl/'+p+'_d'+str(cr[p])+'.dat',ln,l)  
  for L in range(ln):
    dlnCdp[:,:,L,i] = np.dot(iC[L,:,:],dC[:,:,L])
dlnCdp[:2,:2,:5,:] = 0. # remove large scale

for p in cp:
  for j, d in enumerate(fr[p]):

    print p, d
    DP = np.zeros(cn)
    for i, pp in enumerate(cp):
      if pp==p:
        DP[i] = j
        dC = READ_CL(root+'cl/'+pp+'_d'+str(d)+'.dat',ln,l)  
      else:
        di = 0
        if pp=='w_0': di = 1
        if pp=='O_L': di = 2
        if pp=='mnu': di = 1
        DP[i] = di

    # bispectrum Fij
    f  = root+'/l'+str(blmax)+'_rl'+str(rlmax)+'_S4_MG.h5'
    F2 = fsky*loadbispFij(fr,cp,DP,f)
    F2 = F2 + F2.T - np.diag(F2.diagonal())

    # Combined Fisher matrix
    Fl = FisherMatrix(dlnCdp,fsky)
    F1 = np.sum(Fl[:,:,:lmax],axis=2)
    Fl = FisherMatrix(dlnCdp[:2,:2,:,:],fsky)
    F0 = np.sum(Fl[:,:,:lmax],axis=2)
    for i in range(3):
      if i==0: F = F1
      if i==1: F = F0 + F2
      if i==2: F = F1 + F2 #+ Fd
      sig = np.sqrt(np.linalg.inv(F).diagonal())
      print sig

