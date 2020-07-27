    
module GrapheneRibbons
    use constant
    use Utilities, only : First_Available_File_Pointer
    implicit none
#include "../common/ASSERT.h"            
    
    real*8, parameter :: latcons_gh = 2.46d0    !graphene lattice constant Unit: AA
    
    type Latt
    	real*8 :: a  ! Unit:AA, zigzag direction
    	real*8 :: b  ! Armchair direction
    	real*8 :: c  ! z direction
    	real*8 :: p  ! distortion parameter along y direction
    	! The point coordinates to the following Atomic Positions
    	!A(0,p,c/2) B(a/2,b/2-p,c/2) C(0,-p,-c/2) D(a/2,b/2+p,-c/2)
    	real*8 :: t(4)  !Hopping Parameters
    	!t(1) = t12; t(2)=t13; t(3)=t14; t(4)=t23 ;t(5) = p
    end type Latt


    type RibbonInfo_data
      integer :: type   !type = 1: zigzag, type = 2: armchair
      
      real*8 :: latcons  !Unit: AA
      real*8 :: a(2)     !primitive lattice vector
      real*8 :: a2(2)    ! a2 vector to generate the ribbon
      real*8 :: tau(2)   ! B atom bias
      integer :: Num
      real*8 :: gamma0   !nearest hopping energy. Unit: eV
      
      real*8, pointer :: phi(:, :)   ! scalar potential !eV
      
      integer :: Dim
      integer, pointer :: idx(:, :), map(:, :)
      
      real*8 :: width    ! ribbon width between two outmost atoms, Unit:A      
    end type RibbonInfo_data
    
contains
    
    subroutine Read_Lattice_info(Latt,filename)    !从外部读取晶格信息：abc,
		implicit none
		type(Latt), intent(inout) :: Latt
		character*(*), intent(in) :: filename
 		integer :: fn
        character*100 :: chtmp
        character*20 :: ctype, cphi
        real*8 :: Ef
		
		call  First_Available_File_Pointer(fn)


    end subroutine Read_Lattice_info
    
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
        
        rbi%gamma0 = 2.7d0 !
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
            
            HT(n, n) = rbi%phi(i, 1)
            HT(m, m) = rbi%phi(i, 2)
            
            HT(n, m) = ct1
            if(i .gt. 0) HT(n, rbi%map(i - 1, 2)) = rbi%gamma0
            
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))*rbi%gamma0
          
          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            
            HT(n, n) = rbi%phi(i, 1)
            HT(m, m) = rbi%phi(i, 2)            
            
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

subroutine construct_ribbon_H0T(rbi, k, H0T)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: H0T(rbi%Dim, rbi%Dim)
        
        integer :: i, n, m
        complex*16 :: ct1
        
        H0T = 0.d0
        
        if(rbi%type .eq. 1)then
          ct1 = (exp(dcmplx(0.d0, -k)) + 1.d0)*rbi%gamma0

          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            H0T(n, m) = ct1
            if(i .gt. 0) H0T(n, rbi%map(i - 1, 2)) = rbi%gamma0
            
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))*rbi%gamma0
          
          do i = 0, rbi%Num, 1
            n = rbi%map(i, 1); m = rbi%map(i, 2)
            H0T(n, m) = rbi%gamma0
            
            if(i .gt. 0) H0T(n, rbi%map(i - 1, 2)) = rbi%gamma0
            if(i .lt. rbi%Num) H0T(n, rbi%map(i + 1, 2)) = ct1
            
          end do
        endif
        
        do n = 1, rbi%Dim, 1
          do m = n + 1, rbi%Dim, 1            
            H0T(m, n) = dconjg(H0T(n, m))
          end do
        end do
end subroutine construct_ribbon_H0T

subroutine construct_ribbon_HTS(Latt,EP,kx,ky,kz,HT)
implicit none
	real*8, intent(in) :: kx,ky,kz
	type(Latt), intent(in) :: Latt
	type(Hoppings), intent(in) :: t
	complex*16, intent(out) :: HT(4, 4)
	! t12 = ;t13 = ;t14 = ;t23 = ; p = 
	HT(1,2) = t(1) * exp(dcmplx(0.d0,kx*Latt%a/2+ky*Latt%b/2-2*Latt%p*ky))
	HT(1,3) = t(2) * exp(dcmplx(0.d0,-kz*Latt%c-2*Latt%p*ky))	
	HT(1,4) = t(3) * exp(dcmplx(0.d0,kx*Latt%a/2+ky*Latt%b/2-kz*Latt%c))
	HT(2,3) = t(4) * exp(dcmplx(0.d0,-kx*Latt%a/2-ky*Latt%b/2-kz*Latt%c))
	HT(2,4) = t(2) * exp(dcmplx(0.d0,-kz*Latt%c+2*Latt%p*ky))
	HT(3,4) = t(1) * exp(dcmplx(0.d0,kx*Latt%a/2-ky*Latt%b/2+2*Latt%p*ky))

	do n = 1, 4, 1
          do m = n + 1, 4, 1            
            HT(m, n) = dconjg(HT(n, m))
          end do
        end do

end subroutine construct_ribbon_HTS

    subroutine construct_ribbon_velocity(rbi, k, vT)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: vT(rbi%Dim, rbi%Dim, 2)
        
        complex*16 :: H0T(rbi%Dim, rbi%Dim)
        real*8 :: X0(rbi%Dim, 2)
        
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
              vT(rbi%map(i + 1, 2), rbi%map(i, 1), :) = ct            
            end if
          end do
          
        endif
        
!!!!   [X, H]
        call construct_ribbon_H0T(rbi, k, H0T)
        
        do n = 1, rbi%Dim, 1
          X0(n, :) = rbi%idx(1, n) * rbi%a2 
          if(rbi%idx(2, n) .eq. 2) X0(n, :) = x0(n, :) + rbi%tau          
          X0(n, :) = X0(n, :) * cft
        end do
        
        do n = 1, rbi%Dim, 1
          do m = 1, rbi%Dim, 1
            vT(n, m, :) = vT(n, m, :) + (X0(n, :) - X0(m, :))*H0T(n, m)
          end do
        end do
                
    end subroutine construct_ribbon_velocity
    
    
    subroutine construct_ribbon_HT0(rbi, k, HAA, HAB, HBB)
        implicit none
        type(RibbonInfo_data), intent(in) :: rbi
        real*8, intent(in) :: k
        complex*16, intent(out) :: HAA(0:rbi%Num, 0:rbi%Num),HAB(0:rbi%Num, 0:rbi%Num), HBB(0:rbi%Num, 0:rbi%Num)
!        complex*16, intent(out) :: HT(2*(rbi%Num + 1), 2*(rbi%Num + 1))
        
        integer :: i
        complex*16 :: ct1
        
        HAA = 0.d0; HAB = 0.d0; HBB = 0.d0
        
        if(rbi%type .eq. 1)then
          ct1 = (exp(dcmplx(0.d0, -k)) + 1.d0)*rbi%gamma0

          do i = 0, rbi%Num, 1
            HAB(i, i) = ct1
            if(i .ge. 1) HAB(i, i - 1) = rbi%gamma0
            
            HAA(i, i) = rbi%phi(i, 1)
            HBB(i, i) = rbi%phi(i, 2)
          end do
          
        else
          
          ct1 = exp(dcmplx(0.d0, -k))*rbi%gamma0
          
          do i = 0, rbi%Num, 1
            HAB(i, i) = rbi%gamma0
            
            if(i .gt. 0) HAB(i, i - 1) = rbi%gamma0
            if(i .lt. rbi%Num) HAB(i, i + 1) = ct1
            
            HAA(i, i) = rbi%phi(i, 1)
            HBB(i, i) = rbi%phi(i, 2)            
          end do
        endif
        
        
    end subroutine construct_ribbon_HT0

subroutine construct_ribbon_vp(rbi, k, HT, EV, H0T, vT, vP,vP2)
implicit none
type(RibbonInfo_data),intent(in) :: rbi
real*8, intent(in) :: k,EV(:)
complex*16, intent(in) :: HT(:, :), vT(:, :, :),H0T(:,:)
complex*16, intent(out) :: vP(rbi%Dim, rbi%Dim, 2),vP2(rbi%Dim,rbi%Dim,2)

integer :: i,n,m,j

do i = 1,rbi%Dim,1
   do m = 1,rbi%Dim,1
      do n = 1,rbi%Dim,1
        if ((m .ne. n) .and. (m .ne. i) .and. ( n .ne. i)) vP(m,n,:) = vP(m,n,:) + (HT(m,i))/(EV(m)-EV(i))*vT(i,n,:) + vT(m,i,:)*(HT(i,n))/(EV(i)-EV(n))
       do j = 1,rbi%Dim,1
      if ((m .ne. i) .and. (m .ne. j) .and. (m .ne. n) .and. (n .ne. i) .and.(n .ne. j) .and. (i .ne. j)) vP2(m,n,:) = vP2(m,n,:) +(HT(m,i))/(EV(m)-EV(i))* vT(i,j,:) *(HT(j,n))/(EV(j)-EV(n)) 
       end do
        end do
     end do
  end do

!do i = 1,rbi%Dim,1
!  do j = 1,rbi%Dim,1
!   do m = 1,rbi%Dim,1
!      do n = 1,rbi%Dim,1
!        if ((m .ne. i) .and. (m .ne. j) .and. (m .ne. n) .and. (n .ne. i) &
!&.and.(n .ne. j) .and. (i .ne. j)) &
!& vP2(m,n,:) = vP2(m,n,:) + (HT(m,i)-H0T(m,i))/(EV(m)-EV(i))* vT(i,j,:)*&
!& (HT(j,n)-H0T(j,n))/(EV(j)-EV(n))
!        end do
!     end do
!  end do
!end do


end subroutine construct_ribbon_vp

!subroutine construct_ribbon_vp2(rbi, k, HT, EV, H0T, vT, vP2)
!implicit none
!type(RibbonInfo_data),intent(in) :: rbi
!real*8, intent(in) :: k,EV(:)
!complex*16, intent(in) :: HT(:, :), vT(:, :, :),H0T(:,:)
!complex*16, intent(out) :: vP2(rbi%Dim, rbi%Dim, 2)

!integer :: i,n,m,j

!do i = 1,rbi%Dim,1
!  do j = 1,rbi%Dim,1
!   do m = 1,rbi%Dim,1
!      do n = 1,rbi%Dim,1
!        if ((m .ne. i) .and. (m .ne. j) .and. (m .ne. n) .and. (n .ne. i) .and. (n .ne. j) .and. (i .ne. j))  vP2(m,n,:) = vP2(m,n,:) + (HT(m,i)-H0T(m,i))/(EV(m)-EV(i))*vT(i,j,:)*(HT(j,n)-H0T(j,n))/(EV(j)-EV(n))
!        end do
!     end do
!  end do
!end do

!end subroutine construct_ribbon_vp2

    subroutine test()
        implicit none
        type(RibbonInfo_data) :: rbi
        complex*16, pointer :: HT(:, :), vT(:, :,:),H0T(:,:),vP(:,:,:),vP2(:,:,:)
        real*8, pointer :: EV(:)
        real*8 :: k, dk
        integer :: i, NumK, j, fn,p
        
        call Read_Ribbon_info(rbi, "input_ribbonsall_v1")
        allocate(HT(rbi%Dim, rbi%Dim))
        allocate(H0T(rbi%Dim, rbi%Dim))
        allocate(vT(rbi%Dim, rbi%Dim, 2))
        allocate(vP(rbi%Dim, rbi%Dim, 2))
        allocate(vP2(rbi%Dim,rbi%Dim, 2))
        allocate(EV(rbi%Dim))
        
        NumK = 4801
        dk = 2.d0*Pi/dble(NumK)
        vP = 0
        vP2 = 0        

        call First_Available_File_Pointer(fn)
        open(fn, file=trim("test.txt"))
!!        open(1568,file=trim('vP2.dat'))  
!!        open(1551,file=trim('vP.dat'))
        open(6655,file=trim('vT.dat'))
        call Write_Ribbon_info(rbi, fn)
        do i = 0, NumK, 1
        vP = 0
        vP2 = 0
        vT = 0
          k = dk*dble(i)
!        do j = 1, NumK, 1
!        k = 3.6/(NumK-1) * dble(j-1) + 1.2 
!        do j = 1, 41,1
!                k = 4.14 + 0.0075 * dble(j-1)
          call construct_ribbon_HT(rbi, k, HT)          
          call MZHEEV(HT, rbi%Dim, 'V', EV)
          !call construct_ribbon_H0T(rbi, k, H0T)
!!          call MZHEEV(H0T, rbi%Dim, 'V', EV)         
          write(fn, '(1000(e13.6e3,1x))')k, EV
!          call construct_ribbon_velocity(rbi, k, vT)
!          vT(:, :, 1) = MatMul(dconjg(transpose(HT)), MatMul(vT(:, :,1), HT))
!          vT(:, :, 2) = MatMul(dconjg(transpose(HT)), MatMul(vT(:, :,2), HT))
!         HT(:,:) = MatMul(dconjg(transpose(H0T)), MatMul((HT(:,:)-H0T(:,:)), H0T))
!             write(1551,*)k,real(vT(25,28,1)),aimag(vT(25,28,1))
!!          call construct_ribbon_vp(rbi, k, HT, EV, H0T, vT, vP,vP2)
!!             write(1568,*)k,real(vP2(24,26,1)),real(vP2(24,26,2))
!!             write(1551,*)k,real(vP(24,26,1)),real(vP(24,26,2))
!             write(6655,*)k,real(vT(24,24,1)),real(vT(24,25,1)),real(vT(25,24,1)),real(vT(25,25,1))
!        write(fn, '(1000(e13.6e3,1x))')k, vT(rbi%num, rbi%num + 3, 1),vT(rbi%num + 3, rbi%num, 1)
!        write(*, '(1000(e13.6e3,1x))')k, vT(rbi%num, rbi%num + 2, 1),vT(rbi%num + 2,rbi%num, 1)
        end do
        close(fn)
!!        close(1568)
        close(1551)
        close(6655)
        
        deallocate(vP2)
        deallocate(vP)
        deallocate(H0T)
        deallocate(HT)
        deallocate(EV)
        deallocate(vT)
    end subroutine test
end module GrapheneRibbons

program main
    use GrapheneRibbons
    implicit none
    
    type(RibbonInfo_data) :: rbi        
    
    call Read_Ribbon_info(rbi, "input_ribbonsall_v1")
    
    call test
    
end program main

