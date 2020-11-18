# compute shear bispectrum

import numpy as np
import healpy as hp
import curvedsky
from scipy.io import FortranFile
import sim_takahashi

# set parameters
Lmax = 2048
bn   = 20
sn   = 108
bst  = 2
nres = 13

# choose multipoles within a multipole bin
print('define multipole bins')
bp = np.array([np.int(Lmax*(i/bn)) for i in range(bn+1)])
bc = (bp[1:]+bp[:-1])*.5
sL = bp[:2]

# compute binned bispectrum
bl = np.zeros((sn,4,bn))
print('norm equi')
hequi = curvedsky.bispec.bispec_norm(bn,bp,'equi',bst)
print('norm fold')
hfold = curvedsky.bispec.bispec_norm(bn,bp,'fold',bst)
print('norm sque')
hsque = curvedsky.bispec.bispec_norm(bn,bp,'sque',bst,sL)
print('norm isos')
hisos = curvedsky.bispec.bispec_norm(bn,bp,'isos',bst)

for i in range(sn):

    #/// read map ///#
    # Euclid map
    #kmap = hp.fitsfunc.read_map('../data/fullsky/shear/TQUmappy.flagship.n8192.step'+I+'.degrade.n2048.fits')

    # Takahashi sims
    filename = 'data/inp/kappa/allskymap_nres'+str(nres)+'r'+str(i).zfill(3)+'.zs16.mag.dat'
    kmap = sim_takahashi.read_jacobimap(filename)

    # get alm
    print('get kalm',2**nres)
    kalm = curvedsky.utils.hp_map2alm(2**nres,Lmax,Lmax,kmap)

    print('equi')
    bl[i,0,:] = curvedsky.bispec.bispec_bin(bn,bp,Lmax,kalm,'equi',bst) * np.sqrt(4*np.pi)/hequi

    print('fold')
    bl[i,1,:] = curvedsky.bispec.bispec_bin(bn,bp,Lmax,kalm,'fold',bst) * np.sqrt(4*np.pi)/hfold

    print('sque')
    bl[i,2,:] = curvedsky.bispec.bispec_bin(bn,bp,Lmax,kalm,'sque',bst,sL) * np.sqrt(4*np.pi)/hsque

    print('isos')
    bl[i,3,:] = curvedsky.bispec.bispec_bin(bn,bp,Lmax,kalm,'isos',bst) * np.sqrt(4*np.pi)/hisos

# output
np.savetxt('data/sim/sim_bispec_nres'+str(nres)+'_zs16.dat',(np.concatenate((bc[None,:],np.mean(bl,axis=0),np.std(bl,axis=0)))).T,fmt='%.5e')

