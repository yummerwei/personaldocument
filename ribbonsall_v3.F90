
!!! v2: (1) all definition are following H_{nm}(R) = < n R | H | m 0 > 
!!!         H_{nm}(k) = e^{-ik R} H_{nm}(R)
!!!     (2) gamma0 = -2.7
!!! v3: add mass energy in RibbonInfo_data by MassEn
module GrapheneRibbons
    use constant
    use Utilities, only : First_Available_File_Pointer
    implicit none
#include "../common/ASSERT.h"            
    
    real*8, parameter :: latcons_gh = 2.46d0    !graphene lattice constant Unit: AA
    real*8, parameter :: ZERO_EN = 1.d-10       ! 
    type RibbonInfo_data
      integer :: type   !type = 1: zigzag, type = 2: armchair
      
      real*8 :: latcons  !Unit: AA
      real*8 :: a(2)     !primitive lattice vector
      real*8 :: a2(2)    ! a2 vector to generate the ribbon
      real*8 :: tau(2)   ! B atom bias
      integer :: Num
      real*8 :: gamma0   !nearest hopping energy. Unit: eV
      
      real*8 :: MassEn
      
      real*8, pointer :: phi(:, :)   ! scalar potential !eV
      
      integer :: Dim
      integer, pointer :: idx(:, :), map(:, :)
      
      real*8 :: width    ! ribbon width between two outmost atoms, Unit: A      
    end type RibbonInfo_data
    
    type kdata
      integer :: dim
      real*8, pointer :: eV(:)
      complex*16, pointer :: HT(:, :)
      complex*16, pointer :: vT(:, :, :)
      complex*16, pointer :: mT(:, :, :)
      complex*16, pointer :: rT(:, :, :)
      
      complex*16, pointer :: vT_derivk(:, :, :)
      complex*16, pointer :: rT_derivk(:, :, :)
    end type kdata
    
contains
    
    subroutine initialize_kdata(kd)
        implicit none
        type(kdata), intent(inout) :: kd
        nullify(kd%eV)
        nullify(kd%HT)
        nullify(kd%vT)
        nullify(kd%mT)
        nullify(kd%rT)
        nullify(kd%vT_derivk)
        nullify(kd%rT_derivk)
        
    end subroutine initialize_kdata
    
    subroutine Read_Ribbon_info(rbi, filename)
        implicit none
        type(RibbonInfo_data), intent(inout) :: rbi
        character*(*), intent(in) :: filename
        
        integer :: fn
        character*100 :: chtmp
        character*20 :: ctype, cphi
        real*8 :: Ef
        
        integer :: i, j, m
        
        call  First_Available_File_Pointer(fn)
        open(fn, file=trim(filename), status="old")
        read(fn, *)chtmp, ctype
        call NASSERT(chtmp .eq. "RibbonType:") 
        call NASSERT(ctype .eq. "ZigZag" .or. ctype .eq. "ArmChair")
        read(fn, *)chtmp, rbi%Num
        call NASSERT(chtmp .eq. "RibbonWidthNum:")        
        read(fn, *)chtmp, cphi
        call NASSERT(chtmp .eq. "Phiy:")
        if(cphi .eq. "HOMOElf") then
          read(fn, *) Ef
        endif
        close(fn)        
        
        rbi%gamma0 = - 2.7d0 !
        if(ctype .eq. "ZigZag")then
          rbi%type = 1
          rbi%latcons = latcons_gh
          rbi%a2 = (/0.5d0, sqrt(3.d0)/2.d0/)*latcons_gh
          rbi%tau = (rbi%a2 + (/1.d0, 0.d0/)*rbi%latcons) / 3.d0
        else if(ctype .eq. "ArmChair")then
          rbi%type = 2
          rbi%latcons = latcons_gh*sqrt(3.d0)
          rbi%a2 = (/sqrt(3.d0)/2.d0, 0.5d0/)*latcons_gh
          rbi%tau = (/sqrt(3.d0)/3.d0, 0.d0/)*latcons_gh          
        endif
        rbi%a = (/1.d0, 0.d0/)*rbi%latcons
        
        if(cphi .eq. "HOMOElf") then
          allocate(rbi%phi(0:rbi%Num, 2))
          do i = 0, rbi%Num, 1
            rbi%phi(i, 1) = Ef*rbi%a2(2)*dble(i)*1.d-10
            rbi%phi(i, 2) = Ef*(rbi%a2(2)*dble(i) + rbi%tau(2) )*1.d-10
          end do
          rbi%phi = rbi%phi - (rbi%phi(rbi%Num, 2) + rbi%phi(0, 1))/2.d0
        end if
        
        rbi%width = rbi%a2(2)*rbi%Num + rbi%tau(2)
        
        rbi%Dim = (rbi%Num + 1)*2
        allocate(rbi%idx(2, rbi%Dim))
        allocate(rbi%map(0:rbi%Num, 2))
        
        m = 0
        do j = 1, 2, 1
          do i = 0, rbi%Num, 1
            m = m + 1
            rbi%idx(:, m) = (/i, j/)
            rbi%map(i, j) = m
          end do
        end do
        
        rbi%MassEn = 0.d0
    end subroutine Read_Ribbon_info
    
    subroutine Write_Ribbon_info(rbi, fn)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        integer, intent(in) :: fn
        integer :: i
        
        if(rbi%type .eq. 1) then
          write(fn, '(A)')"#Ribbon Information: ZigZag"
        else
          write(fn, '(A)')"#Ribbon Information: ArmChair"
        end if
        write(fn, '(A, I10)')"#Ribbon Atom number: ", (rbi%Num + 1)*2
        write(fn, '(A, f12.6)')"#Ribbon width (A): ", rbi%width
        do i = 0, rbi%Num, 1
          write(fn, '(A, 2(e13.6e3,1x))')"#Ribbon_Phi: ", rbi%phi(i, :)
        end do
    end subroutine Write_Ribbon_info
    
    subroutine construct_ribbon_HT(rbi, k, HT)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: HT(rbi%Dim, rbi%Dim)
        
        integer :: i, n, m
        complex*16 :: ct1
        
        HT = 0.d0
        
        if(rbi%type .eq. 1)then
          ct1 = (exp(dcmplx(0.d0, -k)) + 1.d0)*rbi%gamma0

          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            HT(n, n) = rbi%phi(i, 1) + rbi%MassEn
            HT(m, m) = rbi%phi(i, 2) - rbi%MassEn
            
            HT(n, m) = ct1
            if(i .gt. 0) HT(n, rbi%map(i - 1, 2)) = rbi%gamma0
            
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))*rbi%gamma0
          
          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            HT(n, n) = rbi%phi(i, 1) + rbi%MassEn
            HT(m, m) = rbi%phi(i, 2) - rbi%MassEn
            
            HT(n, m) = rbi%gamma0
            
            if(i .gt. 0) HT(n, rbi%map(i - 1, 2)) = rbi%gamma0
            if(i .lt. rbi%Num) HT(n, rbi%map(i + 1, 2)) = ct1
            
          end do
        endif
        
        do n = 1, rbi%Dim, 1
          do m = n + 1, rbi%Dim, 1            
            HT(m, n) = dconjg(HT(n, m))
          end do
        end do
        
!!$        do n = 1, rbi%Dim, 1
!!$          write(*, '(100(e12.6,1x))')abs(HT(n, :))
!!$        end do
!!$        stop
    end subroutine construct_ribbon_HT
    
    
!!! Direct Fourier transformation v from wannier representaton to Bloch representation
    subroutine construct_ribbon_velocity_w(rbi, k, vT)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: vT(rbi%Dim, rbi%Dim, 2)
        
        integer :: i, n, m
        complex*16 :: ct1
        complex*16 :: cft

        cft = rbi%gamma0 * e_charge * 1.d-10 / dcmplx(0.d0, hbar)
        
        vT = 0.d0
        
        if(rbi%type .eq. 1)then
          ct1 = exp(dcmplx(0.d0, -k)) 

          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            vT(n, m, :) = ( ct1 * (rbi%a - rbi%tau) - rbi%tau) * cft
            if(i .gt. 0) vT(n, rbi%map(i - 1, 2), :) = cft * (rbi%a2 - rbi%tau)
            
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))
          
          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            vT(n, m, :) = -rbi%tau * cft
            
            if(i .gt. 0) vT(n, rbi%map(i - 1, 2), :) = cft * (rbi%a2 - rbi%tau)
            if(i .lt. rbi%Num) vT(n, rbi%map(i + 1, 2), :) = ct1 * cft * (rbi%a - rbi%a2 - rbi%tau)
            
          end do
        endif
        
        do n = 1, rbi%Dim, 1
          do m = n + 1, rbi%Dim, 1            
            vT(m, n, :) = dconjg(vT(n, m, :))
          end do
        end do
        
!!$        do n = 1, rbi%Dim, 1
!!$          write(*, '(100(e12.6,1x))')abs(HT(n, :))
!!$        end do
!!$        stop
    end subroutine construct_ribbon_velocity_w
    
    
    
!!! Direct Fourier transformation v from wannier representaton to Bloch representation
    subroutine construct_ribbon_MassTerm_w(rbi, k, mT)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: mT(rbi%Dim, rbi%Dim, 3) ! (xx, xy, yy)
        
        integer :: i, n, m
        complex*16 :: ct1
        complex*16 :: cft
        
        real*8 :: mtau(3), ma(3), ma2(3)
        
        cft = rbi%gamma0 * e_charge * (1.d-10 / dcmplx(0.d0, hbar))**2
        
        mT = 0.d0
        
        mtau = (/rbi%tau(1) * rbi%tau(1), rbi%tau(1) * rbi%tau(2), rbi%tau(2) * rbi%tau(2)/)
        ma = (/rbi%a(1) * rbi%a(1), rbi%a(1) * rbi%a(2), rbi%a(2) * rbi%a(2)/)
        ma2 = (/rbi%a2(1) * rbi%a2(1), rbi%a2(1) * rbi%a2(2), rbi%a2(2) * rbi%a2(2)/)
        
!        write(*, *)cft, mtau, ma, ma2
        if(rbi%type .eq. 1)then
          ct1 = exp(dcmplx(0.d0, -k)) 

          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            mT(n, m, :) = ( ct1 * (ma - mtau) - mtau) * cft
            if(i .gt. 0) mT(n, rbi%map(i - 1, 2), :) = cft * (ma2 - mtau)
            
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))
          
          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            mT(n, m, :) = -mtau * cft
            
            if(i .gt. 0) mT(n, rbi%map(i - 1, 2), :) = cft * (ma2 - mtau)
            if(i .lt. rbi%Num) mT(n, rbi%map(i + 1, 2), :) = ct1 * cft * (ma - ma2 - mtau)
            
          end do
        endif
        
        do n = 1, rbi%Dim, 1
          do m = n + 1, rbi%Dim, 1            
            mT(m, n, :) = dconjg(mT(n, m, :))
          end do
        end do
        
!!$        do n = 1, rbi%Dim, 1
!!$          write(*, '(100(e12.6,1x))')abs(HT(n, :))
!!$        end do
!!$        stop
    end subroutine construct_ribbon_MassTerm_w
    
    subroutine Get_rk_from_vk(dim, vk, eV, rk)
        implicit none
        integer, intent(in) :: dim
        complex*16, intent(in) :: vk(dim, dim, 2)
        real*8, intent(in) :: eV(dim)
        complex*16, intent(out) :: rk(dim, dim, 2)
        
        integer :: i, j
        complex*16, parameter :: cft = hbar/dcmplx(0.d0, e_charge)
        
        
        rk = 0.d0
        do i = 1, dim, 1
          do j = i+1, dim, 1
            if(eV(j) - eV(i) .gt. ZERO_EN) then
              rk(i, j, :) = vk(i, j, :) * (cft / (eV(i) - eV(j)))
              rk(j, i, :) = dconjg(rk(i, j, :))
            end if
          end do
        end do
        
    end subroutine Get_rk_from_vk
    
!!! dv_k^b/dk : b=(x, y);  k===> x
    subroutine Get_vkderivk(dim, vk, rk, Mk, vkdk)
        implicit none
        integer, intent(in) :: dim
        complex*16, intent(in) :: vk(dim, dim, 2), rk(dim, dim, 2), Mk(dim, dim, 3)
        complex*16, intent(out) :: vkdk(dim, dim, 2)
        
        integer :: i, j
        
        vkdk(:, :, 1) = hbar * Mk(:, :, 1) + CI*(matMul(rk(:, :, 1), vk(:, :, 1)) - matMul(vk(:, :, 1), rk(:, :, 1)))
        vkdk(:, :, 2) = hbar * Mk(:, :, 2) + CI*(matMul(rk(:, :, 1), vk(:, :, 2)) - matMul(vk(:, :, 2), rk(:, :, 1)))
!        vkdk(:, :, 3) = hbar * Mk(:, :, 2) + CI*(matMul(rk(:, :, 2), vk(:, :, 1)) - matMul(vk(:, :, 1), rk(:, :, 2)))
!        vkdk(:, :, 4) = hbar * Mk(:, :, 3) + CI*(matMul(rk(:, :, 2), vk(:, :, 2)) - matMul(vk(:, :, 2), rk(:, :, 2)))
        
    end subroutine Get_vkderivk
    
    
!!! dr_k^b/dk : b=(x, y); k====> x
    subroutine Get_rkderivk(dim, ev, vk, vkdk, rkdk)
        implicit none
        integer, intent(in) :: dim
        real*8, intent(in) :: ev(dim)
        complex*16, intent(in) :: vk(dim, dim, 2), vkdk(dim, dim, 2)
        complex*16, intent(out) :: rkdk(dim, dim, 2)
        
        integer :: i, j
        complex*16, parameter :: cft1 = hbar/dcmplx(0.d0, e_charge), cft2 = hbar**2/dcmplx(0.d0, e_charge**2)
        complex*16 :: Dk(2)
        
        
        rkdk = 0.d0
        
        do i = 1, dim, 1
          do j = i+1, dim, 1
            if(eV(j) - eV(i) .gt. ZERO_EN) then
              Dk = vk(i, i, :) - vk(j, j, :)
              
              rkdk(i, j, :) = vkdk(i, j, :) * (cft1 / (eV(j) - eV(i))) - cft2/(eV(j) - ev(i))**2 * (/&
                & vk(i, j, 1) * Dk(1), &
                & vk(i, j, 2) * Dk(1)/)! , &
                ! & vk(i, j, 1) * Dk(2), &
                ! & vk(i, j, 2) * Dk(2)/)
              
              rkdk(j, i, :) = dconjg(rkdk(i, j, :))
            end if
          end do
        end do
    end subroutine Get_rkderivk
    
    subroutine Get_kAlldata(rbi, k, kd)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi        
        real*8, intent(in) :: k
        type(kdata), intent(inout) :: kd
        
        complex*16 :: HT(rbi%dim, rbi%dim)
        integer :: i
        
        if(.not. associated(kd%HT)) allocate(kd%HT(rbi%dim, rbi%dim))
        
        call construct_ribbon_HT(rbi, k, kd%HT)          
        
        if(.not. associated(kd%eV)) allocate(kd%eV(rbi%dim))
        call MZHEEV(kd%HT, rbi%Dim, 'V', kd%EV)
        
        if(.not. associated(kd%vT)) allocate(kd%vT(rbi%dim, rbi%dim, 2))
        call construct_ribbon_velocity_w(rbi, k, kd%vT)          

        do i = 1, 2, 1
          kd%vT(:, :, i) = MatMul(dconjg(transpose(kd%HT)), MatMul(kd%vT(:, :, i), kd%HT))
        end do
        
        if(.not. associated(kd%mT)) allocate(kd%mT(rbi%dim, rbi%dim, 3))
        call construct_ribbon_MassTerm_w(rbi, k, kd%mT)          

        do i = 1, 3, 1
          kd%mT(:, :, i) = MatMul(dconjg(transpose(kd%HT)), MatMul(kd%mT(:, :, i), kd%HT))
        end do
        
        
        if(.not. associated(kd%rT)) allocate(kd%rT(rbi%dim, rbi%dim, 2))
        call Get_rk_from_vk(rbi%dim, kd%vT, kd%eV, kd%rT)
        
        
        if(.not. associated(kd%vT_derivk)) allocate(kd%vT_derivk(rbi%dim, rbi%dim, 2))
        call Get_vkderivk(rbi%dim, kd%vT, kd%rT, kd%mT, kd%vT_derivk)
        
        
        if(.not. associated(kd%rT_derivk)) allocate(kd%rT_derivk(rbi%dim, rbi%dim, 2))
        call Get_rkderivk(rbi%dim, kd%eV, kd%vT, kd%vT_derivk, kd%rT_derivk)        

    end subroutine Get_KAlldata
    
    subroutine construct_ribbon_velocity(rbi, k, vT)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: vT(rbi%Dim, rbi%Dim, 2)
        
        complex*16 :: HT(rbi%Dim, rbi%Dim)
        complex*16 :: X0(rbi%Dim, 2)
        
        complex*16, parameter :: cft = e_charge*1.d-10/(hbar*CI)
        
        complex*16 :: ct(2)
        
        integer :: i, n, m
        
        
        vT = 0.d0
!!!!  dHk/dk        
        if(rbi%type .eq. 1)then
          ct = exp(dcmplx(0.d0, -k))*rbi%gamma0*cft*rbi%a
          
          do i = 0, rbi%Num, 1
            vT(rbi%map(i, 1), rbi%map(i, 2), :) = ct
            vT(rbi%map(i, 2), rbi%map(i, 1), :) = dconjg(ct)
          end do
          
        else
          
          ct = exp(dcmplx(0.d0, -k))*rbi%gamma0*cft*rbi%a
          
          do i = 0, rbi%Num, 1
            if(i .lt. rbi%Num)then
              vT(rbi%map(i, 1), rbi%map(i + 1, 2), :) = ct           
              vT(rbi%map(i + 1, 2), rbi%map(i, 1), :) = dconjg(ct)
            end if
          end do
          
        endif
!        write(*, *)"### ",(maxval(abs(vT(:, :, i) - dconjg(transpose(vT(:, :, i))))), i = 1, 2)
!!!!   [X, H]
        call construct_ribbon_HT(rbi, k, HT)
        
        do n = 1, rbi%Dim, 1
          X0(n, :) = rbi%idx(1, n) * rbi%a2 
          if(rbi%idx(2, n) .eq. 2) X0(n, :) = x0(n, :) + rbi%tau          
          X0(n, :) = X0(n, :) * cft
        end do
        
        do n = 1, rbi%Dim, 1
          do m = 1, rbi%Dim, 1
            vT(n, m, :) = vT(n, m, :) + (X0(n, :) - X0(m, :))*HT(n, m)
          end do
        end do
                
    end subroutine construct_ribbon_velocity
    
    
    subroutine construct_ribbon_HT0(rbi, k, HAA, HAB, HBB)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: HAA(0:rbi%Num, 0:rbi%Num), HAB(0:rbi%Num, 0:rbi%Num), HBB(0:rbi%Num, 0:rbi%Num)
!        complex*16, intent(out) :: HT(2*(rbi%Num + 1), 2*(rbi%Num + 1))
        
        integer :: i
        complex*16 :: ct1
        
        HAA = 0.d0; HAB = 0.d0; HBB = 0.d0
        
        if(rbi%type .eq. 1)then
          ct1 = (exp(dcmplx(0.d0, -k)) + 1.d0)*rbi%gamma0

          do i = 0, rbi%Num, 1
            HAB(i, i) = ct1
            if(i .ge. 1) HAB(i, i - 1) = rbi%gamma0
            
            HAA(i, i) = rbi%phi(i, 1) + rbi%MassEn
            HBB(i, i) = rbi%phi(i, 2) - rbi%MassEn
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))*rbi%gamma0
          
          do i = 0, rbi%Num, 1
            HAB(i, i) = rbi%gamma0
            
            if(i .gt. 0) HAB(i, i - 1) = rbi%gamma0
            if(i .lt. rbi%Num) HAB(i, i + 1) = ct1
            
            HAA(i, i) = rbi%phi(i, 1) + rbi%MassEn
            HBB(i, i) = rbi%phi(i, 2) - rbi%MassEn
          end do
        endif
        
        
    end subroutine construct_ribbon_HT0
    
    subroutine test()
        implicit none
        type(RibbonInfo_data) :: rbi
        complex*16, pointer :: HT(:, :), vT(:, :, :), vT2(:, :, :), mT(:, :, :)
        real*8, pointer :: EV(:)
        real*8 :: k, dk
        integer :: i, NumK, j, fn
        
        call Read_Ribbon_info(rbi, "input_ribbonsall_v1")
        allocate(HT(rbi%Dim, rbi%Dim))
        allocate(vT(rbi%Dim, rbi%Dim, 2))
        allocate(vT2(rbi%Dim, rbi%Dim, 2))
        allocate(mT(rbi%Dim, rbi%Dim, 3))
        allocate(EV(rbi%Dim))
        
        NumK = 300
        dk = 2.d0*Pi/dble(NumK)
        
        call First_Available_File_Pointer(fn)
        open(fn, file=trim("ap.txt"))
        open(33, file=trim("en.txt"))
        call Write_Ribbon_info(rbi, fn)
        do i = 0, NumK, 1
          k = dk*dble(i)
          call construct_ribbon_HT(rbi, k, HT)          
          call MZHEEV(HT, rbi%Dim, 'V', EV)
          
          write(33, '(1000(e13.6e3,1x))')k, EV
          
          call construct_ribbon_velocity(rbi, k, vT)
          call construct_ribbon_velocity_w(rbi, k, vT2)          
!          write(*, *)maxval(abs(vT-vT2)), maxval(abs(vT2))
          
          call construct_ribbon_MassTerm_w(rbi, k, mT)          
          
          
          vT(:, :, 1) = MatMul(dconjg(transpose(HT)), MatMul(vT(:, :, 1), HT))

          vT(:, :, 2) = MatMul(dconjg(transpose(HT)), MatMul(vT(:, :, 2), HT))
          
!          write(fn, '(1000(e13.6e3,1x))')k, abs(vT(95, 99, 1)), abs(vT(95,99,2)), EV(99)-EV(95)
          
        end do
        close(fn)
        close(fn)
        deallocate(HT)
        deallocate(EV)
        deallocate(vT)
    end subroutine test
    
    
    
    subroutine test1()
        implicit none
        type(RibbonInfo_data) :: rbi
        type(kdata) :: kd
        real*8 :: k, dk
        integer :: i, NumK, j, fn, n, m, nn
        real*8, pointer :: nc(:)
        
        call Read_Ribbon_info(rbi, "input_ribbonsall_v1")
        
        call initialize_kdata(kd)
        
        allocate(nc(rbi%Dim))
        
        NumK = 300
        dk = 2.d0*Pi/dble(NumK)
        
        call First_Available_File_Pointer(fn)
        open(fn, file=trim("ap.txt"))
        call Write_Ribbon_info(rbi, fn)
        n = rbi%num+1; m = rbi%num+4
        
        nc = 0.d0
        
        do i = 0, NumK, 1
          k = dk*dble(i)

          call Get_kAlldata(rbi, k, kd)
          
          do j = 1, rbi%Dim, 1
            if(kd%eV(j) .lt. 1.d-3) then
              nc(:) = nc(:) + abs(kd%HT(:, j))**2
            endif
          end do
          
          write(fn, '(100(e12.6,1x))')k, kd%eV(n), kd%ev(m), ((kd%rT_derivk(n, m, nn)*kd%rT(m, n, j), nn = 1, 2), j= 1,2)
          ! write(fn, '(I5,1x,100(e12.6,1x))')1, k, kd%eV(:)
          ! write(fn, '(I5,1x,100(e12.6,1x))')2, k, kd%vT(n, m, :)
          ! write(fn, '(I5,1x,100(e12.6,1x))')3, k, kd%mT(n, m, :)
          ! write(fn, '(I5,1x,100(e12.6,1x))')4, k, kd%rT(n, m, :)
          ! write(fn, '(I5,1x,100(e12.6,1x))')5, k, kd%vT_derivk(n, m, :)
          ! write(fn, '(I5,1x,100(e12.6,1x))')6, k, kd%rT_derivk(n, m, :)
          
        end do
        close(fn)
        
        do i = 0, rbi%num, 1
          write(fn, *)i, nc(rbi%map(i, 1)), nc(rbi%map(i, 2))
        end do
    end subroutine test1
    
    
    subroutine test_elf()
        implicit none
        type(RibbonInfo_data) :: rbi
        complex*16, pointer :: HT(:, :)
        real*8, pointer :: EV(:, :), EV0(:, :), EVx(:, :)
        real*8 :: k, dk, elf
        integer :: i, NumK, j, fn, m
        
        call Read_Ribbon_info(rbi, "input_ribbonsall_v1")
        allocate(HT(rbi%Dim, rbi%Dim))
        
        NumK = 1000        
        allocate(EV(rbi%Dim, 0:NumK))
        allocate(EV0(rbi%Dim, 0:NumK))
        EV = 0.d0
        EV0 = 0.d0
        
        dk = Pi/dble(NumK)
        
        call First_Available_File_Pointer(fn)
        open(fn, file=trim("en.txt"))
        call Write_Ribbon_info(rbi, fn)
        do j = 0, 300, 1
          elf = dble(j) * 40.d0 
          do m = 0, rbi%Num, 1
            rbi%phi(m, 1) = elf*rbi%a2(2)*dble(m)*1.d-10
            rbi%phi(m, 2) = elf*(rbi%a2(2)*dble(m) + rbi%tau(2) )*1.d-10
          end do
          rbi%phi = rbi%phi - (rbi%phi(rbi%Num, 2) + rbi%phi(0, 1))/2.d0
          
          if(j .eq. 0) then
            EVx => EV0
          else
            EVx => EV
          end if
          
          write(*, *)j
!          do i = 0, NumK, 1
          do i = 0, NumK, 1
            k = dk * dble(i)
!            i = 0
!            k = 2.51d0 
            call construct_ribbon_HT(rbi, k, HT)          
            call MZHEEV(HT, rbi%Dim, 'N', EVx(:, i))
            
            if(j .gt. 0) then
              write(fn, '(1000(e19.12e3,1x))')k, elf, EV(26, i), EV0(26, i)
            end if
                      
          end do
          write(fn, *)
        end do
        
        close(fn)
        deallocate(HT)
        deallocate(EV)
        deallocate(EV0)
    end subroutine test_elf
    
    
end module GrapheneRibbons
