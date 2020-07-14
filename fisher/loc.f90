
module local
  use myconst, only: pi
  use cmblbisp
  implicit none
  private pi

contains


subroutine fisher_bisp(eL,k,Pk,Cl,fac,abc,wp,ck,der,F,pb)
  implicit none
  !I/O
  logical, intent(in) :: pb
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: k(:,:,:,:), Pk(:,:,:,:), Cl(:), fac(:,:,:), wp(:,:,:,:), ck(:,:,:,:), abc(:,:,:,:,:), der(:)
  double precision, intent(out) :: F(:)
  !internal
  integer :: zn, lmax, l1, l2, l3, i, p, q, pm
  double precision :: cov, Del, F2(1:3), l1l2, l2l3, l3l1, al1, al2, al3
  double precision, allocatable :: bisp(:,:), dbisp(:)

  zn   = size(k,dim=1)
  lmax = size(k,dim=2)

  do l1 = eL(1), eL(2)

    if (mod(l1,100)==0) write(*,*) l1

    do l2 = l1, lmax
      do l3 = l2, lmax

        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        if (mod(l1+l2+l3,2)==1) cycle
        Del = 1d0
        if (l1==l2.and.l2/=l3) Del = 2d0
        if (l1/=l2.and.l2==l3) Del = 2d0
        if (l1==l2.and.l2==l3) Del = 6d0
        cov  = Del*cl(l1)*cl(l2)*cl(l3)

        ! computing bispectrum for each cosmological parameter shifts
        allocate(bisp(2,2),dbisp(2));  bisp=0d0; dbisp=0d0
        do p = 1, 2 !cosmological parameters
          do pm = 1, 2 !p+dp, p-dp

            !1) post-Born
            if (pb) call bisp_postborn(l1,l2,l3,wp(:,:,p,pm),ck(:,:,p,pm),bisp(p,pm))

            !2) LSS
            do i = 1, zn

              ! F2 kernel
              call F2_Kernel([k(i,l1,p,pm),k(i,l2,p,pm),k(i,l3,p,pm)],abc(:,i,l1,p,pm),abc(:,i,l2,p,pm),F2(1))
              call F2_Kernel([k(i,l2,p,pm),k(i,l3,p,pm),k(i,l1,p,pm)],abc(:,i,l2,p,pm),abc(:,i,l3,p,pm),F2(2))
              call F2_Kernel([k(i,l3,p,pm),k(i,l1,p,pm),k(i,l2,p,pm)],abc(:,i,l3,p,pm),abc(:,i,l1,p,pm),F2(3))
              bisp(p,pm) = bisp(p,pm) + fac(i,p,pm) * (F2(1)*Pk(i,l1,p,pm)*Pk(i,l2,p,pm) + F2(2)*Pk(i,l2,p,pm)*Pk(i,l3,p,pm) + F2(3)*Pk(i,l3,p,pm)*Pk(i,l1,p,pm))

            end do
          end do
        end do

        ! Fisher matrix
        bisp = bisp * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
        do p = 1, 2
          dbisp(p) = (bisp(p,1)-bisp(p,2))/(2d0*der(p))
        end do
        F(l1) = F(l1) + dbisp(1)*dbisp(2) / cov

        deallocate(bisp,dbisp)
      end do
    end do
  end do

end subroutine fisher_bisp


end module local

