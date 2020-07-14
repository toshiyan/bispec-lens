#//////////////////////////////////////////////////#
# * plot ellipse
#//////////////////////////////////////////////////#

import numpy as np, matplotlib as mpl, h5py
from matplotlib.pyplot import *
from matplotlib.patches import Ellipse

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
fr['O_L'] = np.array([.5,.3,.1])
fr['A_s'] = np.array([.1])
fr['n_s'] = np.array([.005])
fr['w_0'] = np.array([.05,.03])
fr['mnu'] = np.array([.3,.2,.1])
fr['tau'] = np.array([.5,.4,.3])

# cosmological parameters
cp = ['w_0','mnu','Omh','Obh','O_L','n_s','A_s','tau']
cn = len(cp)
root = '../../../data/forecast/model0/'

# multipole
lmin=2; lmax=4000; ln=lmax-lmin+1; l=np.linspace(lmin,lmax,ln)
llmax = 2000
rlmax = lmax

# signal covariance 
CL = READ_CL(root+'cl/fid.dat',ln,l)

# experiments: fsky and noise covariance 
ac2rad = np.pi/180./60.
fsky = 1.
nlpp = (1.*ac2rad/2.72e6)**2*np.exp(-l*(l+1)*(3.*ac2rad)**2/(8.*np.log(2.)))
nldd = (np.loadtxt(root+'../nldd/S4_s1_t3_rlmax'+str(rlmax)+'.dat')).T[1]

NL = np.zeros((3,3,ln))
NL[0,0,:] = nlpp/2.
NL[1,1,:] = nlpp
NL[2,2,:] = nldd[:ln]
#NL[2,2,llmax:ln] = nldd[llmax:ln]*100000.

# primary CMB + 2pt
cr = {'Obh':.1,'Omh':.005,'O_L':.03,'A_s':.1,'n_s':.005,'w_0':.005,'mnu':.2,'tau':.4}
iC = np.array([np.linalg.inv(CL[:,:,L]+NL[:,:,L]) for L in range(ln)])
dlnCdp = np.zeros((3,3,ln,cn))
for i, p in enumerate(cp):
  dC = READ_CL(root+'cl/'+p+'_d'+str(cr[p])+'.dat',ln,l)  
  for L in range(ln):
    dlnCdp[:,:,L,i] = np.dot(iC[L,:,:],dC[:,:,L])
dlnCdp[:2,:2,:5,:] = 0.

DP = np.zeros(cn)
for i, pp in enumerate(cp):
  di = 0
  if pp=='w_0': di = 0
  if pp=='Omh': di = 1
  if pp=='mnu': di = 1
  if pp=='O_L': di = 2
  DP[i] = di

# Fisher matrix (F_l = (2l+1)/2 * Tr(dlnC/dpi dlnC/dpj))
Fl = FisherMatrix(dlnCdp,fsky)
F1 = np.sum(Fl[:,:,:lmax],axis=2)
Fl = FisherMatrix(dlnCdp[:2,:2,:,:],fsky)
F0 = np.sum(Fl[:,:,:lmax],axis=2)

# bispectrum Fij
F2 = np.zeros((cn,cn))
f2 = fsky*loadbispFij(fr,cp,DP,root+'/l'+str(llmax)+'_rl'+str(rlmax)+'_zn100_all.h5')
f2 = f2 + f2.T - np.diag(f2.diagonal())
F2[:,:] = f2

# ellipse
a = 0
b = 4
xlim(-1-.25,-1+.25)
ylim(0.72-0.2,.72+.2)
sig = np.zeros((5,cn))
ax = subplot(111)
cr = ['k','c','b','g','b']
ls = ['solid','dashed','dotted','solid','solid']
circ = Circle((-1,0.72),.0005,color='k')
ax.add_artist(circ)
lab = [r'primary CMB',r'2pt',r'3pt',r'2pt+3pt']
xlabel(r'$w$',fontsize=20)
ylabel(r'$\Omega_\Lambda$',fontsize=20)

# DESI
#mu = np.array([4.1,1.7,.88,.55,.38,.28,.21,.18,.18,.17,.16,.14,.15,.16,.19,.28,.41,.52])*1e-3
#fp = np.zeros((cn,18))
#fp[3,:], fp[2,:], fp[4,:], fp[0,:] = np.loadtxt(root+'fz.dat',unpack=True,usecols=(1,2,3,4))
#Fd = np.zeros((cn,cn))
#for i, p in enumerate(cp):
#  for j, q in enumerate(cp):
#    Fd[i,j] = np.sum(fp[i,:]*fp[j,:]/mu**2)

for i in range(4):
  if i==0: F = F0
  if i==1: F = F1
  if i==2: F = F0 + F2
  if i==3: F = F1 + F2
  iF = np.linalg.inv(F)
  if i>0:  sig[i,:] = np.sqrt(iF.diagonal())

  sF = np.zeros((2,2))
  sF[0,0] = iF[a,a]
  sF[0,1] = iF[a,b]
  sF[1,1] = iF[b,b]
  sF[1,0] = sF[0,1]

  lam, v = np.linalg.eig(sF)
  lam = np.sqrt(lam)*1.516575089
  ell = Ellipse(xy=(-1,.72),width=2*lam[0],height=2*lam[1],angle=np.rad2deg(np.arccos(v[0,0])))
  ell.set_facecolor('none')
  ell.set_edgecolor(cr[i])
  ell.set_linestyle(ls[i])
  ell.set_linewidth(2)
  ax.add_patch(ell)
  ell.set(label=lab[i])

legend(loc=2,frameon=False)
print sig[1,:]/sig[3,:], sig[1,0:2], sig[3,0:2]

#np.savetxt("sigma.txt",sig[1:,:].T,fmt=4*'%f1 ')
#savefig('ellipse.png')
show()

