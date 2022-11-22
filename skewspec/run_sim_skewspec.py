# compute shear bispectrum
import numpy as np
import curvedsky
import sim_takahashi as simrt

# set parameters
D = '/global/u1/t/toshiyan/scratch/RTsim/'
Lmin = 100
Lmax = 2000
sn = 1
nr = 13

# compute binned bispectrum
nside = 2**nr
npix  = 12*nside**2

zid = [9,16,21,25,66]

for zn0 in zid:
    zs0 = simrt.zid2zs(zn0)

    for zn1 in zid:
        if zn0==zn1: continue
        zs1 = simrt.zid2zs(zn1)

        cl = np.zeros((sn,Lmax+1))

        for theta in [10.,5.,2.,0.]:

            ll = np.linspace(0,Lmax,Lmax+1)
            wl = np.exp(-ll*(ll+1.)*(theta*np.pi/10800.)**2/16./np.log(2))

            for i in range(sn):

                # Takahashi sims
                kmap0 = simrt.read_jacobimap(D+'/allskymap_nres'+str(nr)+'r'+str(i).zfill(3)+'.zs'+str(zn0)+'.mag.dat')

                # get alm
                print('kalm',nside)
                kalm0 = curvedsky.utils.hp_map2alm(nside,Lmax,Lmax,kmap0)
                kalm0[:Lmin,:] = 0.
                kalm0 *= wl[:,None]
                kmap0 = curvedsky.utils.hp_alm2map(npix,Lmax,Lmax,kalm0)

                if zn0==zn1:
                    salm = curvedsky.utils.hp_map2alm(nside,Lmax,Lmax,kmap0**2)
                    cl[i,:] = curvedsky.utils.alm2cl(Lmax,salm,kalm0)
                else:
                    kmap1 = simrt.read_jacobimap(D+'allskymap_nres'+str(nr)+'r'+str(i).zfill(3)+'.zs'+str(zn1)+'.mag.dat')
                    kalm1 = curvedsky.utils.hp_map2alm(nside,Lmax,Lmax,kmap1)
                    kalm1[:Lmin,:] = 0.
                    kalm1 *= wl[:,None]
                    kmap1 = curvedsky.utils.hp_alm2map(npix,Lmax,Lmax,kalm1)
                    salm = curvedsky.utils.hp_map2alm(nside,Lmax,Lmax,kmap1**2)
                    cl[i,:] = curvedsky.utils.alm2cl(Lmax,kalm0,salm)

            # output
            np.savetxt('sim_nres'+str(nr)+'_zs'+str(zs0)+'x'+str(zs1)+'_l'+str(Lmin)+'-'+str(Lmax)+'_t'+str(theta)+'.dat',np.array((np.linspace(0,Lmax,Lmax+1),np.mean(cl,axis=0),np.std(cl,axis=0))).T,fmt='%.5e')

