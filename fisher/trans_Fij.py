#///////////////////////////////////////////////////////////////////////////////////////////////////#
# Save F_ij elements to .h5
#///////////////////////////////////////////////////////////////////////////////////////////////////#
import numpy as np, h5py

# values of delta p
fr = {}
fr['Obh'] = np.array([.1])
fr['Omh'] = np.array([.007,.005,.003])
fr['O_L'] = np.array([.2,.1,.07,.05])
fr['A_s'] = np.array([.1])
fr['n_s'] = np.array([.005])
fr['w_0'] = np.array([.05,.03])
fr['mnu'] = np.array([.3,.2,.1])

# maximum multipole of CMB
rlmax = 4000

# zs
zs = 1088.69

#cp = ['w_0','mnu','Omh','Obh','O_L','n_s','A_s']

for M in ['TR','GM']:
  for nldd in ['','S4']:
    for lmax in [100,300,500,700,1000,1500]:

      # filename to be saved
      if nldd=='': f ='l'+str(lmax)+'_rl'+str(rlmax)+'_CV_'+Md+'.h5'
      if nldd!='': f ='l'+str(lmax)+'_rl'+str(rlmax)+'_'+nldd+'_'+Md+'.h5'

      # tag of output data from f90 code
      etag = M+str(lmax)+str(nldd)+str(ztra)+str(zs)

      # making h5
      with h5py.File(f,'w') as hf:
        for p in range(len(cp)):
          print p
          for q in range(p,len(cp)):
            if p==q:
              for dp in fr[cp[p]]:
                tag = cp[p]+'_dp'+str(dp)
                l, F = np.loadtxt('data_'+tag+etag+'.dat',unpack=True)
                hf.create_dataset(tag,data=np.sum(F))
            else:
              for dp in fr[cp[p]]:
                for dq in fr[cp[q]]:
                  tag = cp[p]+'_dp'+str(dp)+'_'+cp[q]+'_dp'+str(dq)
                  l, F = np.loadtxt('data_'+tag+etag+'.dat',unpack=True)
                  hf.create_dataset(tag,data=np.sum(F))

      # check
      with h5py.File(f,'r') as hf:
        print len(hf.keys())
        for p in range(len(cp)):
          for q in range(p,len(cp)):
            if p==q:
              for dp in fr[cp[p]]:
                print cp[p], np.array(hf.get(cp[p]+'_dp'+str(dp)))
            else:
              for dp in fr[cp[p]]:
                for dq in fr[cp[q]]:
                  tag = cp[p]+'_dp'+str(dp)+'_'+cp[q]+'_dp'+str(dq)
                  print cp[p], cp[q], np.array(hf.get(tag))

