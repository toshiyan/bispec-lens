!///////////////////////////////////////////////////////////////////////!
! * CMB Lensing Bispectrum Fisher matrix
!///////////////////////////////////////////////////////////////////////!

program main
  use readfile, only: set_params_file, read_prm, read_str, read_int, read_dbl, read_val, read_log
  use myconst,  only: pi
  use myutils,  only: GLdxs, GLpoints, loadtxt, linspace, savetxt, str, GAUSS_LEGENDRE_PARAMS, gl_initialize, gl_finalize
  use myfunc,   only: cosmoparams
  use cmblbisp
  use local
  implicit none
  character(LEN=128) :: pkfile, cps(2), dstr(2)
  logical :: pb
  integer :: i, j, l, zn, kn, lmax, lmin, pm, p
  double precision :: zs, zmin, zmax, h, Om, Omh, Ov, w0, wa, mnu, ders(2)
  double precision, allocatable :: F(:),a(:),z(:),dz(:),cl(:),datl(:,:),nldd(:,:),ell(:),fac(:,:,:),pl(:,:,:,:),wp(:,:,:,:),k(:,:,:,:),ck(:,:,:,:),abc(:,:,:,:,:), ckk(:)
  type FiducialParams
    double precision :: Obh, Omh, O_L, n_s, A_s, w_0, w_a ,mnu
  end type FiducialParams
  type(fiducialparams) :: fp
  type(cosmoparams) :: cp
  type(gauss_legendre_params) :: gl

  call set_params_file

  ! numerical parameters
  zn   = read_int('zn')
  zs   = read_dbl('zs')
  kn   = 1380
  lmin = read_int('lmin')
  lmax = read_int('lmax')
  zmin = read_dbl('zmin')
  zmax = read_dbl('zmax')
  ders = [read_dbl('der1'),read_dbl('der2')]
  dstr = [read_str('der1'),read_str('der2')]
  cps  = [read_str('p1'),read_str('p2')]
  pb   = read_log('pb')

  ! check errors
  if (cps(2)=='') stop 'should specify parameter'

  ! set fiducial cosmology
  fp%Obh = read_dbl('Obh')
  fp%Omh = read_dbl('Omh')
  fp%O_L = read_dbl('O_L')
  fp%n_s = read_dbl('n_s')
  fp%A_s = read_dbl('A_s')
  fp%w_0 = read_dbl('w_0')
  fp%w_a = read_dbl('w_a')
  fp%mnu = read_dbl('mnu')

  ! derivatives
  do p = 1, 2
    if (cps(p)=='O_L') ders(p) = ders(p)*fp%Omh
    if (cps(p)=='n_s') ders(p) = ders(p)*fp%n_s
    if (cps(p)=='mnu') ders(p) = ders(p)*fp%mnu
  end do

  !* precomputing interpolation points for z
  allocate(z(zn),a(zn),dz(zn))
  call gl_initialize(GL,zn,1d-15)
  z  = GLpoints([zmin,zmax],GL%z) ! points
  dz = GLdxs([zmin,zmax],GL%w)    ! width
  a  = 1d0/(1d0+z)                ! scale factor

  allocate(abc(3,zn,lmax,2,2),k(zn,lmax,2,2),Pl(zn,lmax,2,2),fac(zn,2,2),wp(zn,lmax,2,2),ck(zn,lmax,2,2),ckk(lmax))
  Pl=0d0; wp=0d0; fac=0d0; k=0d0; ck=0d0; abc=1d0; ckk=0d0

  do p = 1, 2

    do pm = 1, 2

      !* fiducial cosmological parameters
      Omh = fp%Omh
      Ov  = fp%O_L
      Om  = 1d0 - Ov
      w0  = fp%w_0
      wa  = fp%w_a
      mnu = fp%mnu

      !* add diff
      if(cps(p)=='Omh')  Omh = Omh + (3-2*pm)*ders(p)*Omh  !ders is in ln Omh
      if(cps(p)=='O_L')  Ov  = Ov  + (3-2*pm)*ders(p)
      if(cps(p)=='w_0')  w0  = w0  + (3-2*pm)*ders(p)
      if(cps(p)=='w_a')  wa  = wa  + (3-2*pm)*ders(p)
      if(cps(p)=='mnu')  mnu = mnu + (3-2*pm)*ders(p)

      !* set of derived parameters
      cp%Om  = 1d0 - Ov
      cp%Ov  = Ov
      h      = dsqrt(Omh/(1d0-Ov))
      cp%H0  = h*1d2
      cp%nu  = mnu/(93.14d0*Omh)
      cp%w0  = w0
      cp%wa  = wa

      !* read template k and Pk (if an explicit derivative is available, fiducial P(k) is loaded)
      if (pm==1) pkfile = trim(read_str('dPk'))//'/pk_'//trim(cps(p))//'_d'//trim(dstr(p))//'_+.dat'
      if (pm==2) pkfile = trim(read_str('dPk'))//'/pk_'//trim(cps(p))//'_d'//trim(dstr(p))//'_-.dat'

      !* load linear Pk at z=0
      allocate(datl(zn+1,kn))
      call loadtxt(pkfile,datl(1:2,:),fsize=[2,kn])
      ! k/h -> k ( e.g. 1h Mpc^-1 -> 0.72 Mpc^-1 )
      datl(1,:) = datl(1,:)*h
      datl(2,:) = datl(2,:)/h**3

      !* prep bisp quantities
      call prep_lens_bispectrum(z,dz,zs,cp,datl(1,:),datl(2,:),read_str('model'),k(:,:,p,pm),pl(:,:,p,pm),fac(:,p,pm),abc(:,:,:,p,pm),wp(:,:,p,pm),ck(:,:,p,pm))
      call prep_lens_aps(z,dz,zs,cp,datl(1,:),datl(2,:),ckk)

      deallocate(datl)

    end do
  end do

  !* clkk + noise
  write(*,*) 'cl'
  allocate(datl(6,lmax),nldd(2,lmax),cl(lmax),ell(lmax)); datl=0d0; nldd=0d0; ell=1d0
  ell = linspace(1,lmax)
  if (read_val('nldd'))  call loadtxt(read_str('nldd'),nldd(1:2,2:lmax),fsize=[2,lmax-1])
  if (zs<=1000d0) then
    cl = ckk
  else
    call loadtxt(read_str('cldd'),datl(1:6,2:lmax),fsize=[6,lmax-1])
    cl = datl(5,:)/4d0*(1d0+1d0/ell)**2 + nldd(2,:)*ell**2*(1d0+ell)**2/4d0
  end if
  deallocate(nldd,datl)

  !* bispectrum Fisher matrix
  allocate(F(lmax));  F=0d0

  write(*,*) 'Fisher'
  call fisher_bisp([lmin,lmax],k,Pl,cl,fac,abc,wp,ck,ders,F,pb)

  ! save
  call savetxt('data_'//trim(read_str('tag'))//'.dat',linspace(1,lmax),F,ow=.true.)
  write(*,*) sum(F)

  ! free up memory
  call gl_finalize(GL)
  deallocate(z,a,dz,F,cl,ckk,Pl,wp,ck,k,fac,abc,ell)


end program main

