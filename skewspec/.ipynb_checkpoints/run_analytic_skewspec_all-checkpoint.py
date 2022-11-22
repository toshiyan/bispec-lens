import basic, numpy as np, healpy as hp, os


def compute_skewspec(f,cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pb):
    #if os.path.exists(f): return
    skew = basic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,bn+1,ols,lmin,lmax,k,pk0,kn,theta=theta,pb=pb)
    #skew = basic.bispec.skewspeclens(cpmodel,model,z,dz,zn,zs,bn+1,ols,lmin,lmax,k,pk0,kn,pb=pb)
    np.savetxt(f,np.array((ols,skew[0,:],skew[1,:],skew[2,:])).T)
    

cpmodel = 'modelw'
model = 'RT'
#model = 'GM'

D = 'data/'
k, pk0 = np.loadtxt(D+cpmodel+'/Pk/Pklin.dat',unpack=True)
kn = np.size(k)

lmin = 100
lmax = 2000
#lmax = 1024
bn = 25
ols = np.linspace(lmin,lmax,bn+1,dtype=np.int)
zss = [0.5078,1.0334,1.5345,2.0548,1100]
#zss = [0.5078]

# auto skew spec
for zs0 in zss:

    for zs1 in zss:
        
        zs = [zs0,zs1]
        zmin, zmax = 0.0001, min(zs0,zs1)
        zn = 100
        z, dz = basic.bispec.zpoints(zmin,zmax,zn)

        #for theta in [10.,5.,2.,0.]:
        for theta in [10.,5.,2.]:

            for pb, PB in [(False,''),(True,'_pb')]:
            
                if zs0==zs1: #auto skew spec
                
                    f = D+'skewspec/Sl_'+model+'_zs'+str(zs0)+'_zn'+str(zn)+'_l'+str(lmin)+'-'+str(lmax)+'_t'+str(theta)+PB+'.dat'
                    print(zs0,theta,pb)
                    compute_skewspec(f,cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pb)

                else: # cross
                    f = D+'skewspec/Sl_'+model+'_zs'+str(zs0)+'x'+str(zs1)+'_zn'+str(zn)+'_l'+str(lmin)+'-'+str(lmax)+'_t'+str(theta)+PB+'.dat'
                    print(zs0,zs1,theta,pb)
                    print(zs)

                compute_skewspec(f,cpmodel,model,z,dz,zn,zs,bn,ols,lmin,lmax,k,pk0,kn,theta,pb)

