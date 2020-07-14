# pipeline submit
import os, numpy as np
from quicksub import *

def set_pfile(tag,cp,d,model='',lmax=500,nldd='S4',zs=1088.69):
  f = 'sub_'+tag+'.sh'
  p = 'sub_p_'+tag+'.ini'
  set_sbatch_params(f,tag,mem='8G',t='0-10:00')
  add('./exe '+p,f)
  add('lmin = 2',p,True)
  add('lmax = '+str(lmax),p)
  if nldd=='S4': add('nldd = ../../../data/forecast/nldd/S4_s1_t3_rlmax4000.dat',p)
  add('cldd = ../../../data/forecast/model0/cl/fid.dat',p)
  add('pklin = ../../../data/forecast/model0/Pk/Pklin.dat',p)
  add('dPk = ../../../data/forecast/model0/dPk/',p)
  add('zs = '+str(zs),p)
  add('zmin = 0.001',p)
  add('zmax = '+str(np.amin([zs,30.])),p)
  add('zn = 200',p)
  add('model = '+model,p)
  add('pb = T',p)
  add('Obh = 0.02226',p)
  add('Omh = 0.14220419154',p)
  add('O_L = 0.6879845',p)
  add('A_s = 2.13e-9',p)
  add('n_s = 0.96530',p)
  add('w_0 = -1.0000',p)
  add('w_a = 0.0',p)
  add('mnu = 0.06',p)
  add('p1 = '+cp[0],p)
  add('p2 = '+cp[1],p)
  add('der1 = '+str(d[0]),p)
  add('der2 = '+str(d[1]),p)
  add('tag = '+tag,p)
  #os.system('sbatch '+f)
  os.system('sh '+f)
  os.system('rm -rf '+f+' '+p)


# values of delta p
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

n = 0
#for nldd in ['','S4']:
for nldd in ['']:
  #for lmax in [100,300,500,700,1000,1500]:
  for lmax in [500,1500]:
    for M in ['TR','GM']:
      for zs in [1088.69]: 
          print lmax, M, nldd, zs

          # usual Fij matrix
          for p in range(len(cp)):
            print p
            for q in range(p,len(cp)):
              if p==q:
                for dp1 in fr[cp[p]]:
                  tag = cp[p]+'_dp'+str(dp1)+M+str(lmax)+nldd+str(ztra)+str(zs)
                  set_pfile(tag,[cp[p],cp[p]],[dp1,dp1],model=M,lmax=lmax,nldd=nldd,zs=zs,ztra=ztra)
                  n = n + 1
              else:
                for dp1 in fr[cp[p]]:
                  for dp2 in fr[cp[q]]:
                    tag = cp[p]+'_dp'+str(dp1)+'_'+cp[q]+'_dp'+str(dp2)+M+str(lmax)+nldd+str(ztra)+str(zs)
                    set_pfile(tag,[cp[p],cp[q]],[dp1,dp2],model=M,lmax=lmax,nldd=nldd,zs=zs,ztra=ztra)
                    n = n + 1

          print n

