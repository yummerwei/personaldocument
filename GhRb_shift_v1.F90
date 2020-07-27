
program main
    use GrapheneRibbons
    use Math_Operation, only : FermiF
    use GaussianDelta
    implicit none
#include "../common/ASSERT.h"     
    type(RibbonInfo_data) :: rbi
    type(kdata) :: kd
    complex*16, pointer :: HT(:, :)
    real*8, pointer :: EV(:), fk(:), EVout(:, :)
    real*8, pointer :: inw(:)
    integer :: NumK, NumW
    real*8 :: k, dk, dw
    integer :: i, j, fn, s1, s2
    real*8 :: tempt, beta, We, Ws
    integer, pointer :: bnd_idx(:)
    complex*16, pointer :: sht(:), inj(:), sht_bndp(:, :, :), inj_bndp(:, :, :)
    
    character*100 :: rbifile, chtmp
    real*8 :: ChemMu, GaussBroadEn
    real*8 :: pre_eps
    logical :: sgn, ifrefb
    integer :: refine_bnd, bnds, bnde
    
    fn = 22
    open(fn, file="input_GhRb_shift_v1.dat")
    read(fn, *)chtmp, rbifile
    call NASSERT(chtmp .eq. "RibbonInfoFile:")        
    
    read(fn, *)chtmp, NumK
    call NASSERT(chtmp .eq. "DivisionNumk:")
    
    dk = 2.d0*Pi/dble(NumK)
        
    
    read(fn, *)chtmp, ChemMu, Tempt
    call NASSERT(chtmp .eq. "InitialDopingTem:")
    
    beta = e_charge / (kB * Tempt)
    
    read(fn, *)chtmp, NumW, Ws, We
    call NASSERT(chtmp .eq. "EnergyDivision:")

    dw = (We - Ws)/dble(NumW-1)
    allocate(inw(1:numw))
    do i = 1, numw, 1
      inw(i) = Ws + dble(i-1)*dw
    end do

    read(fn, *)chtmp, GaussBroadEn
    call NASSERT(chtmp .eq. "GaussBroadEn:")
    
    read(fn, *)chtmp, pre_eps
    call NASSERT(chtmp .eq. "PRE_EPS:")
    pre_eps = max(pre_eps, 1.d-12)
    
    read(fn, *)chtmp, ifrefb, refine_bnd
    call NASSERT(chtmp .eq. "Bands_Contrib:")

    close(fn)
    
    call Initialize_GaussianDelta(GaussBroadEn)
    
    call Read_Ribbon_info(rbi, trim(rbifile))
    call initialize_kdata(kd)
    
    
    kd%dim = rbi%dim
    
    allocate(sht(1:numw)); sht = 0.d0
    allocate(inj(1:numw)); inj = 0.d0
    
    if(ifrefb)then
      bnds = rbi%dim/2 - refine_bnd
      bnde = rbi%dim/2 + 1 + refine_bnd
      allocate(bnd_idx(bnds:bnde))
      do i = bnds, rbi%dim/2, 1
        bnd_idx(i) = i - rbi%dim/2 - 1
      end do
      do i = rbi%dim/2 + 1, bnde, 1
        bnd_idx(i) = i - rbi%dim/2
      end do
      
      allocate(inj_bndp(1:numw, bnds:bnde, bnds:bnde)) ; inj_bndp = 0.d0
      allocate(sht_bndp(1:numw, bnds:bnde, bnds:bnde)) ; sht_bndp = 0.d0
    end if
    
    allocate(HT(rbi%Dim, rbi%Dim))
    allocate(EV(rbi%Dim))
    
    allocate(EVout(rbi%Dim, 0:NumK-1))
    
    allocate(fk(rbi%Dim))
    
    do i = 0, NumK - 1, 1
      k = dk * dble(i)
      call construct_ribbon_HT(rbi, k, HT)
      call MZHEEV(HT, rbi%Dim, 'N', EV(:))
      EVout(:, i) = EV
      
      do j = 1, rbi%dim, 1
        fk(j) = FermiF(beta*(EV(j) - ChemMu))
      end do
      
      call CheckIfCalk(kd%dim, ev, fk, numw, inw, pre_eps, sgn)

      if(sgn) then
        call Get_kAlldata(rbi, k, kd)
        
        call get_results(kd, fk, numw, inw, inj, sht, ifrefb, bnds, bnde, inj_bndp, sht_bndp)
        
!        call get_injection(kd, fk, numw, inw, inj)
!        call get_shift(kd, fk, numw, inw, sht)
          
      end if
!      if(i .eq. NumK/2) write(22, *)
    end do
    
    open(22, file="rbiinfo.txt")
    call Write_Ribbon_info(rbi, 22)
    do i = 0, Numk - 1, 1
      write(22, '(1000(e13.6e3,1x))')dble(i)/dble(Numk), EVout(:, i)
    end do    
    close(22)
        
    inj = inj * (pi * e_charge**2/hbar)*(dk/(2.d0*Pi)/(rbi%LatCons*1.d-10)) * 2.d0 ! spin
    
    sht = sht * dcmplx(0.d0, -1.d0)*(pi * e_charge**2/hbar)*(dk/(2.d0*Pi)/(rbi%LatCons*1.d-10)) * 2.d0
    
    if(ifrefb) then
      inj_bndp = inj_bndp * (pi * e_charge**2/hbar)*(dk/(2.d0*Pi)/(rbi%LatCons*1.d-10)) * 2.d0 ! spin
      sht_bndp = sht_bndp * dcmplx(0.d0, -1.d0)*(pi * e_charge**2/hbar)*(dk/(2.d0*Pi)/(rbi%LatCons*1.d-10)) * 2.d0
    end if
    
    open(33, file="shift.txt")
    
    write(33, *)"#(m) ", rbi%width*1.d-10
    write(33, *)"#NN ", rbi%Num + 1
    write(33, *)"#Temperature(K) ", Tempt
    write(33, *)"#Chemical Potential(eV) ", ChemMu
    write(33, *)"#GaussianBroadening(eV) ", GaussBroadEn
    write(33, *)"#Ef(real), w(real), current injection coefficient(complex), shift conductivity(complex)"
    do i = 1, numw, 1
      write(33, '(100(e13.6e3,1x))')(rbi%phi(2, 1)-rbi%phi(1, 1))/(rbi%a2(2)*dble(1)*1.d-10), inw(i),  inj(i), sht(i)
    end do
    write(33, *)
    close(33)
    
    if(ifrefb) then
      open(33, file="shift_refine.txt")
    
      write(33, *)"#(m) ", rbi%width*1.d-10
      write(33, *)"#NN ", rbi%Num + 1
      write(33, *)"#Temperature(K) ", Tempt
      write(33, *)"#Chemical Potential(eV) ", ChemMu
      write(33, *)"#GaussianBroadening(eV) ", GaussBroadEn
      write(33, *)"#Column Meaning:"
      write(33, *)"#Ef(real) ", 1
      write(33, *)"#w(real)  ", 2
      write(33, *)"#INJECTION COEFFICIENTS: " 
      j = 3
      do s1 = bnds+1, bnde, 1
        write(33, '("# ", $)')
        do s2 = bnds, s1-1, 1
          write(33, '("[",I3,"=>",I3,":(",I3,",",I3,")]", $)')bnd_idx(s2), bnd_idx(s1), j, j+1
          j = j+2
        end do
        write(33, *)
      end do
      
      write(33, *)"#SHIFT CONDUCTIVITY: "
      do s1 = bnds+1, bnde, 1
        write(33, '("# ", $)')
        do s2 = bnds, s1-1, 1
          write(33, '("[",I3,"=>",I3,":(",I3,",",I3,")]", $)')bnd_idx(s2), bnd_idx(s1), j, j+1
          j = j+2
        end do
        write(33, *)
      end do
            
      do i = 1, numw, 1
        write(33, '(100000(e13.6e3,1x))')(rbi%phi(2, 1)-rbi%phi(1, 1))/(rbi%a2(2)*dble(1)*1.d-10), inw(i),  &
          & ((inj_bndp(i, s1, s2), s2 = bnds, s1-1), s1=bnds+1,bnde),&
          & ((sht_bndp(i, s1, s2), s2 = bnds, s1-1), s1=bnds+1,bnde)
      end do
    end if
    write(33, *)
    close(33)
    
end program main

subroutine CheckIfCalk(dim, ev, fk, numw, inw, EPS, sgn)
    use GrapheneRibbons
    use GaussianDelta
    implicit none
    integer, intent(in) :: dim, numw    
    real*8, intent(in) :: ev(dim), fk(dim), inw(numw), EPS
    logical, intent(out) :: sgn
    
    integer :: s1, s2, iw
    real*8 :: dfe, df


    
    do s1 = 1, dim, 1
      do s2 = 1, dim, 1
        df = fk(s2) - fk(s1)
        if(abs(df) .lt. EPS) cycle
        
        do iw = 1, numw, 1
          dfe = df * GaussianDelta_f (eV(s1) - ev(s2) - inw(iw))
          if(abs(dfe) .gt. EPS) then
            sgn = .true.
            return
          endif
          
        end do
        
      end do
    end do
    
    sgn = .false.
end subroutine CheckIfCalk

subroutine get_results(kd, fk, numw, inw, inj, sht, ifrefb, bs, be, inj_bndp, sht_bndp)
    use GrapheneRibbons
    use GaussianDelta
    implicit none
    type(kdata), intent(in) :: kd
    integer, intent(in) :: numw, bs, be
    real*8, intent(in) :: fk(kd%dim), inw(numw)
    logical, intent(in) :: ifrefb
    complex*16, intent(inout) :: inj(numw)
    complex*16, intent(inout) :: inj_bndp(numw, bs:be, bs:be)
    complex*16, intent(inout) :: sht(numw)
    complex*16, intent(inout) :: sht_bndp(numw, bs:be, bs:be)
    
    integer :: s1, s2, iw, j
    real*8 :: df, dfe
    
    complex*16 :: ctmp1, ctmp2
    
    do s1 = 1, kd%dim, 1
      do s2 = 1, kd%dim, 1
        df = fk(s2) - fk(s1)
        if(abs(df) .lt. 1.d-8) cycle
        
        do iw = 1, numw, 1
          dfe = df * GaussianDelta_f (kd%eV(s1) - kd%ev(s2) - inw(iw))
          
          if(abs(dfe) .lt. 1.d-9) cycle
          
          ctmp1 =  dfe*(kd%vT(s1, s1, 1) - kd%vT(s2, s2, 1)) &
            & * ( kd%rT(s2,s1,2)*kd%rT(s1,s2,1) &
            & - kd%rT(s2,s1,1)*kd%rT(s1, s2, 2) )
          
          inj(iw) = inj(iw) + ctmp1
          
          ctmp2 =  dfe * ( kd%rT(s1,s2,1)*kd%rT_derivk(s2,s1,2) &
            & - kd%rT(s1,s2,2)*kd%rT(s2, s1, 1) )
          
          sht(iw) = sht(iw) + ctmp2
          
          if(ifrefb) then
            if(s1 .ge. bs .and. s1 .le. be) then
              if(s2 .ge. bs .and. s2 .le. be) then
                inj_bndp(iw, s1, s2) = inj_bndp(iw, s1, s2) + ctmp1
                sht_bndp(iw, s1, s2) = sht_bndp(iw, s1, s2) + ctmp2
              endif
            end if
          end if
          
        end do
      end do
    end do
    
    
end subroutine get_results

! subroutine get_injection(kd, fk, numw, inw, inj, inj_bndp)
!     use GrapheneRibbons
!     use GaussianDelta
!     implicit none
!     type(kdata), intent(in) :: kd
!     integer, intent(in) :: numw    
!     real*8, intent(in) :: fk(kd%dim), inw(numw)
!     complex*16, intent(inout) :: inj(numw)
!     complex*16, optional, intent(inout) :: inj_bndp(numw, kd%dim/2-3:kd%dim/2+4, kd%dim/2-3:kd%dim/2+4)
    
!     integer :: s1, s2, iw, j, bs, be
!     real*8 :: df, dfe
    
!     logical :: sbp, sbp1, sbp2
!     complex*16 :: ctmp
    
!     if(present(inj_bndp)) then
!       sbp = .true.
!       bs = kd%dim/2 - 3
!       be = kd%dim/2 + 4
!     else
!       sbp = .false.
!     end if
    
!     do s1 = 1, kd%dim, 1
!       do s2 = 1, kd%dim, 1
!         df = fk(s2) - fk(s1)
!         if(abs(df) .lt. 1.d-8) cycle
        
!         do iw = 1, numw, 1
!           dfe = df * GaussianDelta_f (kd%eV(s1) - kd%ev(s2) - inw(iw))
          
!           if(abs(dfe) .lt. 1.d-9) cycle
!           ctmp =  dfe*(kd%vT(s1, s1, 1) - kd%vT(s2, s2, 1)) &
!             & * ( kd%rT(s2,s1,2)*kd%rT(s1,s2,1) &
!             & - kd%rT(s2,s1,1)*kd%rT(s1, s2, 2) )
          
!           inj(iw) = inj(iw) + ctmp
          
!           if( sbp ) then
!             if(s1 .ge. bs .and. s1 .le. be .and. s2 .ge. bs .and. s2 .le. be) then
!               inj_bndp(iw, s1, s2) = inj_bndp(iw, s1, s2) + ctmp
!             endif
            
!           end if
          
!         end do
!       end do
!     end do
    
    
! end subroutine get_injection



! subroutine get_shift(kd, fk, numw, inw, sht, sht_bndp)
!     use GrapheneRibbons
!     use GaussianDelta
!     implicit none
!     type(kdata), intent(in) :: kd
!     integer, intent(in) :: numw    
!     real*8, intent(in) :: fk(kd%dim), inw(numw)
!     complex*16, intent(inout) :: sht(numw)
!     complex*16, optional, intent(inout) :: sht_bndp(numw, kd%dim/2-3:kd%dim/2+4, kd%dim/2-3:kd%dim/2+4)
    
!     integer :: s1, s2, iw, j,  bs, be
!     real*8 :: df, dfe
    
!     logical :: sbp
!     complex*16 :: ctmp
    
!     if(present(sht_bndp)) then
!       sbp = .true.
!       bs = kd%dim/2 - 3
!       be = kd%dim/2 + 4
!     else
!       sbp = .false.
!     end if
    
!     do s1 = 1, kd%dim, 1
!       do s2 = 1, kd%dim, 1
!         df = fk(s2) - fk(s1)
!         if(abs(df) .lt. 1.d-8) cycle
        
!         do iw = 1, numw, 1
!           dfe = df * GaussianDelta_f (kd%eV(s1) - kd%ev(s2) - inw(iw))
          
!           if(abs(dfe) .lt. 1.d-9) cycle
!           ctmp =  dfe * ( kd%rT(s1,s2,1)*kd%rT_derivk(s2,s1,2) &
!             & - kd%rT(s1,s2,2)*kd%rT(s2, s1, 1) )
          
!           sht(iw) = sht(iw) + ctmp
          
!           if(sbp) then
!             if(s1 .ge. bs .and. s1 .le. be .and. s2 .ge. bs .and. s2 .le. be) then
!               sht_bndp(iw, s1, s2) = sht_bndp(iw, s1, s2) + ctmp
!             end if
!           endif
          
!         end do
!       end do
!     end do
    
    
! end subroutine get_shift
