
import basic
import numpy as np
import healpy as hp

cpmodel = 'modelw'
model = 'RT'
#model = 'GM'
#model = 'SC'

zmin, zmax = 0.0001, 1088.69
#zmin, zmax = 0.0001, 1.0334
zn = 150
zs = [zmax,zmax,zmax]

calc = 'bispecsnr'
#calc = 'bispecbin'

btype = 'ggk'
#btype = 'kkk'

lmin = 1
lmax = 2048
olmin = lmin
olmax = lmax

bn = 20

D = 'data/'
L = np.linspace(0,olmax,olmax+1)
k, pk0 = np.loadtxt( D+cpmodel+'/Pk/Pklin.dat', unpack=True )
kn = np.size(k)

z, dz = basic.bispec.zpoints(zmin,zmax,zn)

if btype == 'kkk': 
    dNdz = None
else:
    dNdz = basic.galaxy.dndz_sf(z,2.,1.,1.)

# bispectrum
if calc == 'bispec':
    bl0, pb0 = basic.bispec.bispeclens('equi',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl1, pb1 = basic.bispec.bispeclens('fold',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl2, pb2 = basic.bispec.bispeclens('sque',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    bl3, pb3 = basic.bispec.bispeclens('angl',cpmodel,model,z,dz,zs,lmin,lmax,k,pk0,btype=btype,dNdz=dNdz)
    np.savetxt('test_zs'+str(zs)+'.dat',np.array((L[1:],bl0,bl1,bl2,bl3,pb0,pb1,pb2,pb3)).T)

# binned bispectrum
if calc == 'bispecbin':
    bl, pb = {}, {}
    for shap in ['equi','fold','sque','angl']:
        bc, bl[shap], pb[shap] = basic.bispec.bispeclens_bin(shap,cpmodel,model,z,dz,zs,lmin,lmax,bn,k,pk0,btype=btype,dNdz=dNdz)
    np.savetxt('data/modelw/bl/bl_'+model+'_b'+str(bn)+'_zs1.dat',np.array((bc,bl['equi'],bl['fold'],bl['sque'],bl['angl'],pb['equi'],pb['fold'],pb['sque'],pb['angl'])).T)

# total bispectrum SNR
if calc == 'bispecsnr':

    # load aps
    lmaxs = np.array([100,500,1000,2000,3000])
    #lmaxs = np.array([100,500,1000])
    snr = np.zeros(len(lmaxs))

    for i, lmax in enumerate(lmaxs):
        
        L = np.linspace(0,lmax,lmax+1)
        cldd = np.zeros(lmax+1)
        cldd[2:] = np.loadtxt(D+cpmodel+'/cl/fid.dat',unpack=True)[4][:lmax-1]
        nldd = np.zeros(lmax+1)
        #nldd[2:] = np.loadtxt(D+'nldd/advact_s6_t1_rlmax4000.dat',unpack=True)[1][:lmax-1]
        nldd[2:] = np.loadtxt(D+'nldd/so/kappa_deproj0_sens2_16000_lT30-3000_lP30-5000.dat',unpack=True)[7][:lmax-1]
        cldd[2:] = cldd[2:]/4.*(1.+1./L[2:])**2 + nldd[2:]#*L[2:]**2*(L[2:]+1.)**2/4.
        
        if btype=='kkk':
            clgg = None
        else:
            clgg = basic.bispec.cl_flat(cpmodel,z,dz,[zmax,zmax],lmax,k,pk0,pktype='T12',cltype='gg',dNdz=dNdz)
            #cgg = sum(dNdz**2*Pl(:,l)*dz*kernel_cgg(z,cp)) 
            clgg += (np.pi/10800.)**2/30.

        snr[i] = basic.bispec.bispeclens_snr(cpmodel,model,z,dz,zs,2,lmax,cldd,k,pk0,btype=btype,dNdz=dNdz,cgg=clgg)
        print(lmax,snr[i])

    np.savetxt('snr_'+btype+'_kappa_deproj0_sens2_16000_lT30-3000_lP30-5000_'+str(zn)+'.dat',np.array((lmaxs,snr)).T)


