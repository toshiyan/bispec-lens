# compute shear bispectrum

import numpy as np
import healpy as hp
import curvedsky
import sim_takahashi

# set parameters
Lmax = 2048
bn   = 20
sn   = 108
nres = 13

# choose multipoles within a multipole bin
print('define multipole bins')
bp = np.array([np.int(Lmax*(i/bn)) for i in range(bn+1)])
bc = (bp[1:]+bp[:-1])*.5

# compute binned bispectrum
cl = np.zeros((sn,bn))

for i in range(sn):

    #/// read map ///#
    # Takahashi sims
    filename = 'data/inp/kappa/allskymap_nres'+str(nres)+'r'+str(i).zfill(3)+'.zs16.mag.dat'
    kmap = sim_takahashi.read_jacobimap(filename)

    # get alm
    kalm = curvedsky.utils.hp_map2alm(2**nres,Lmax,Lmax,kmap)
    cl[i,:] = curvedsky.utils.alm2bcl(bn,Lmax,kalm)
    
# output
np.savetxt('data/sim/sim_aps_nres'+str(nres)+'_zs16.dat',(np.concatenate((bc[None,:],np.mean(cl,axis=0),np.std(cl,axis=0)))).T,fmt='%.5e')

