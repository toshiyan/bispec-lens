!///////////////////////////////////////////////////////////////////////!
! * CMB Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

program main
  use myconst
  use readfile
  use myfunc
  use myutils
  use cmblbisp
  use utils_galaxy
  implicit none
  character(4) :: bisptype
  integer :: i, l, zn, kn, lmax, imax
  double precision :: h, zran(2), alpha, beta, zm
  double precision, allocatable, dimension(:)     :: z, dz, ckk, cgg, fac, dNdz
  double precision, allocatable, dimension(:,:)   :: blll, rcl, datl, Pl, wp, k, ck, nldd, snr, bl
  double precision, allocatable, dimension(:,:,:) :: abc
  type(gauss_legendre_params) :: gl
  type(cosmoparams) :: cp

  call set_params_file

  cp%Om = read_dbl('Om')
  cp%H0 = read_dbl('H0')
  cp%Ov = 1d0 - cp%Om
  cp%w0 = read_dbl('w0')
  cp%wa = read_dbl('wa')
  cp%nu = read_dbl('mnu')/(93.14d0*(cp%H0/100d0)**2*cp%Om)
  h     = cp%H0/100d0

  zn = read_int('zn')
  call set_prm('zran',zran,[0.001d0,40d0])
  kn = 1380

  imax = read_int('imax')

  alpha = read_dbl('alpha')
  beta  = read_dbl('beta')
  zm    = read_dbl('zm')

  bisptype = read_str('bisptype')

  !* load linear Pk at z=0
  allocate(datl(2,kn))
  call loadtxt(read_str('pklin'),datl(1:2,:),fsize=[2,kn])
  datl(1,:) = datl(1,:)*h    ! k
  datl(2,:) = datl(2,:)/h**3 ! linear Pk at z=0

  !* precomputing interpolation points for z
  allocate(z(zn),dz(zn))
  select case (read_int('zspace')) 
  case(0)
    z  = linspace(zran,zn)
    dz = z(2)-z(1)
  case(1)
    call gl_initialize(gl,zn,1d-15)
    z  = glpoints(zran,gl%z) ! points
    dz = gldxs(zran,gl%w)    ! width
  end select

  allocate(dNdz(zn))
  dNdz = nz_SF_arr(z,alpha,beta,z0_sF(alpha,beta,zm))  !; call savetxt('dNdz.dat',z,dNdz,ow=.true.)

  allocate(snr(3,imax)); snr=0d0

  do i = 1, imax

    lmax = 100 + (i-1)*200

    !if (lmax/=900) cycle

    allocate(Pl(zn,lmax),k(zn,lmax),fac(zn),abc(3,zn,lmax),wp(zn,lmax),ck(zn,lmax));  Pl=0d0; k=0d0; fac=0d0; abc=1d0; wp=0d0; ck=0d0
    call prep_lens_bispectrum(z,dz,read_dbl('zs'),cp,datl(1,:),datl(2,:),read_str('model'),k,pl,fac,abc,wp,ck,btype=bisptype)

    !* z-kernel
    if (bisptype=='gkk') fac  = fac*dNdz
    if (bisptype=='ggk') fac  = fac*dNdz**2

    !allocate(blll(6,zn),bl(3,lmax)); blll=0d0
    !call bisp_equi([1,lmax],k,Pl,fac,abc,wp,ck,bl(2,:),btype='LSS',blll=blll(2:6,:))
    !blll(1,:) = z
    !call savetxt('bl.dat',blll,ow=.true.)
    !stop

    !* clkk + noise
    write(*,*) 'cl'
    allocate(rcl(6,lmax),nldd(4,lmax),cgg(lmax),ckk(lmax),bl(3,lmax)); rcl=0d0; nldd=1d0; cgg=0d0; bl=0d0
    call loadtxt(read_str('nldd'),nldd(1:4,2:lmax),fsize=[4,lmax-1])
    call loadtxt(read_str('clfile'),rcl(1:6,2:lmax),fsize=[6,lmax-1])
    ckk = rcl(5,:)/4d0*(1d0+1d0/rcl(1,:))**2 + nldd(2,:)*rcl(1,:)**2*(1d0+rcl(1,:))**2/4d0
    do l = 2, lmax
      cgg(l) = sum(dNdz**2*Pl(:,l)*dz*kernel_cgg(z,cp)) + ac2rad**2/read_dbl('ngal')
    end do
    !call savetxt('clgg.dat',rcl(1,:),cgg,ckk,ow=.true.)
    deallocate(rcl,nldd)

    !* bispectrum
    !bl(1,:) = linspace(1,lmax)
    !call bisp_equi([1,lmax],k,Pl,fac,abc,wp,ck,bl(2,:),btype='LSS')
    !call bisp_fold([1,lmax/2],k,Pl,fac,abc,wp,ck,bl(3,:),btype='LSS')
    !call savetxt('bl.dat',bl,ow=.true.)

    !* SNR
    snr(1,i)   = lmax
    snr(2:3,i) = snr_xbisp([2,lmax],zn,k,Pl,cgg,ckk,fac,abc,wp,ck,bisptype)
    snr(2,i)   = dsqrt(snr(2,i))
    write(*,*) snr(1,i), snr(2,i)

    deallocate(k,Pl,cgg,ckk,fac,abc,wp,ck,bl)

  end do

  call savetxt('snr_'//trim(bisptype)//'_zn'//str(zn)//'.dat',snr,ow=.true.)

  if (read_int('zspace')==1) call gl_finalize(gl)
  deallocate(z,dz,datl,dNdz,snr)

end program main

