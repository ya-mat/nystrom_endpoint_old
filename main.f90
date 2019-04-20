! -*- coding: utf-8 -*-
program main
  use lp_lap
  use mod_force_raise
  implicit none

  real(8) :: rad
  real(8),allocatable :: x(:, :)
  real(8),allocatable :: xn(:, :)
  complex*16,allocatable :: lp1(:, :)
  complex*16,allocatable :: lp2(:, :)
  complex*16,allocatable :: u(:)
  complex*16,allocatable :: kai(:)
  integer,allocatable :: ipiv(:)
  integer :: n
  integer :: info
  real(8),parameter :: pi = 3.1415926535897932384626433832795028841971d0
  real(8),parameter :: arctwopi = 0.159154943091895d0
  complex*16,parameter :: iunit = dcmplx(0d0, 1d0)
  integer :: exterior

  real(8) :: xm(2)
  real(8) :: th
  real(8) :: dth
  real(8) :: nth
  real(8) :: hy
  integer :: i
  integer :: j
  integer :: i0
  integer :: i1
  integer :: j0
  integer :: j1

  real(8) :: y1(2)
  real(8) :: y2(2)
  complex*16 :: result
  real(8) :: sx(2)
  complex*16 :: beshy1, beshy2
  real(8) :: doty1
  real(8) :: doty2
  real(8) :: dotn
  real(8) :: rr1(2)
  real(8) :: rr2(2)
  real(8) :: r1 !rr1の長さ
  real(8) :: r2 !rr2の長さ
  real(8) :: dnrm2 !blas
  integer :: ierr, nz

  real(8) :: s, r
  real(8) :: RR(2)
  complex*16 :: besh(2)
  integer,parameter :: bunten = 4
  real(8) :: p0(bunten), w(bunten)

  integer :: nr
  real(8) :: hr
  integer,parameter :: b_num = 2 ! 2 or 4
  real(8) :: beta(b_num)
  integer,parameter :: col_num = 20 ! 4 or 12 or 20
  real(8) :: wg(col_num)
  real(8),allocatable :: wn(:)
  real(8),allocatable :: w1(:)
  real(8),allocatable :: w2(:)
  integer :: tmpi

  real(8) :: tmp
  real(8) :: tmp2
  real(8) :: k_1
  real(8) :: phi
  real(8) :: cyr(2)
  real(8) :: cyi(2)
  real(8) :: ec
  real(8),allocatable :: en(:)
  integer :: e_ite2
  real(8),allocatable :: facgam(:)
  integer :: smooth_star
  real(8) :: rad1
  real(8) :: rad2
  integer :: daen_flag
  integer :: boundary_flag

  result = dcmplx(0d0, 0d0)

  hy = 1.1d0

  write(*,*) 'integral log(x)', hy*(1d0-log(hy*0.5d0))*arctwopi

  ! gauss
  call givedata_p0_and_w(p0, w, bunten)

  xm = (/0.0d0, 0d0/)
  y1 = (/-0.5d0*hy, 0d0/)
  y2 = (/0.5d0*hy, 0d0/)

  result = dcmplx(0d0, 0d0)
  do i=1,bunten
     s = 0.5d0 + 0.5d0*p0(i)
     RR = y1 + (y2 - y1) * s
     RR = xm - RR

     r = dnrm2(2, RR(1), 1)
     result = result + w(i)*(-log(r)*arctwopi)*hy*0.5d0
  end do
  write(*,*) 'gauss', result

  ! the corrected trapezoidal rule
  if(b_num .eq. 2) then
     beta(1) = 0.5694444444444444d-1
     beta(2) = -0.7638888888888889d-2
  else if(b_num .eq. 4) then
     beta(1) = 0.6965636022927689d-1
     beta(2) = -0.1877177028218695d-1
     beta(3) = 0.3643353174603175d-2
     beta(4) = -0.3440531305114639d-3
  end if

  if(col_num .eq. 4) then
     wg(1) = 0.7518812338640025d0
     wg(2) = -0.6032109664493744d0
     wg(3) = 0.1073866830872157d1
     wg(4) = -0.7225370982867850d0
  else if(col_num .eq. 12) then
     wg(1) = 0.2051970990601252d1
     wg(2) = -0.7407035584542865d1
     wg(3) = 0.1219590847580216d2
     wg(4) = -0.1064623987147282d2
     wg(5) = 0.4799117710681772d1
     wg(6) = -0.8837770983721025d0

     wg(7) = 0.2915391987686506d1
     wg(8) = -0.8797979464048396d1
     wg(9) = 0.1365562914252423d2
     wg(10) = -0.1157975479644601d2
     wg(11) = 0.5130987287355766d1
     wg(12) = -0.9342187797694916d0
  else if(col_num .eq. 20) then
     wg(1) = 0.3256353919777872d1
     wg(2) = -0.2096116396850468d2
     wg(3) = 0.6872858265408605d2
     wg(4) = -0.1393153744796911d3
     wg(5) = 0.1874446431742073d3
     wg(6) = -0.1715855846429547d3
     wg(7) = 0.1061953812152787d3
     wg(8) = -0.4269031893958787d2
     wg(9) = 0.1009036069527147d2
     wg(10) = -0.1066655310499552d1

     wg(11) = 0.4576078100790908d1
     wg(12) = -0.2469045273524281d2
     wg(13) = 0.7648830198138171d2
     wg(14) = -0.1508194558089468d3
     wg(15) = 0.1996415730837827d3
     wg(16) = -0.1807965537141134d3
     wg(17) = 0.1110467735366555d3
     wg(18) = -0.4438764193424203d2
     wg(19) = 0.1044548196545488d2
     wg(20) = -0.1100328792904271d1
  end if

  result = dcmplx(0d0, 0d0)
  nr = 10
  hr = (hy*0.5d0)/(nr-1)

  do i = 1, nr - 2
     result = result - log(abs(hr*i))*arctwopi
  end do
  result = result - log(abs(hr*(nr-1)))*arctwopi*0.5d0

  result = result*hr*2d0

  !collection
  do i = 1, col_num/2
     result = result - 2d0*hr*(wg(i)+wg(i+col_num/2))*log(abs(hr*i))*arctwopi
  end do

  do i = 1, b_num
     result = result - 2d0*hr*beta(i)*(log(abs(hr*(nr-1-i))) - log(abs(hr*(nr-1+i))))*arctwopi
  end do

  write(*,*) 'the corrected trapezoidal rule', result

  result = dcmplx(0d0, 0d0)

  nr = 12
  allocate(wn(-nr+1-b_num:nr-1+b_num))
  write(*,*) 'size(wn)', size(wn)
  hr = (hy*0.5d0)/(nr-1)

  wn = hr
  wn(nr-1) = 0.5d0*hr
  wn(-nr+1) = 0.5d0*hr
  wn(0) = 0d0
  wn(1) = wn(1) + 0.2051970990601252d1*2d0*hr
  wn(2) = wn(2) - 0.7407035584542865d1*2d0*hr
  wn(3) = wn(3) + 0.1219590847580216d2*2d0*hr
  wn(4) = wn(4) - 0.1064623987147282d2*2d0*hr
  wn(5) = wn(5) + 0.4799117710681772d1*2d0*hr
  wn(6) = wn(6) - 0.8837770983721025d0*2d0*hr
  wn(-1) = wn(-1) + 0.2915391987686506d1*2d0*hr
  wn(-2) = wn(-2) - 0.8797979464048396d1*2d0*hr
  wn(-3) = wn(-3) + 0.1365562914252423d2*2d0*hr
  wn(-4) = wn(-4) - 0.1157975479644601d2*2d0*hr
  wn(-5) = wn(-5) + 0.5130987287355766d1*2d0*hr
  wn(-6) = wn(-6) - 0.9342187797694916d0*2d0*hr
  wn(nr-1-1) = wn(nr-1-1) + 0.5694444444444444d-1*hr
  wn(nr-1+1) = -0.5694444444444444d-1*hr
  wn(-nr+1+1) = wn(-nr+1+1) + 0.5694444444444444d-1*hr
  wn(-nr+1-1) = -0.5694444444444444d-1*hr
  wn(nr-1-2) = wn(nr-1-2) - 0.7638888888888889d-2*hr
  wn(nr-1+2) = 0.7638888888888889d-2*hr
  wn(-nr+1+2) = wn(-nr+1+2) - 0.7638888888888889d-2*hr
  wn(-nr+1-2) = 0.7638888888888889d-2*hr

  do i = -nr+1-b_num, nr-1+b_num
!     write(*,*) 'i, wn(i)', i, wn(i)
  end do

  do i = -nr+1-b_num, -1
     result = result - wn(i)*log(abs(hr*i))*arctwopi
  end do
  do i = 1, nr-1+b_num
     result = result - wn(i)*log(abs(hr*i))*arctwopi
  end do

  write(*,*) 'the corrected trapezoidal rule2', result

  deallocate(wn)

! Laplace ------------------------------------------------

!  write(6,*) '# n, rad'
!  open(1000, file='input', status='old')
!  read(1000,*) n, rad
!  close(1000)

  rad = 2d0
  n = 200

  allocate(x(2, n))
  allocate(xn(2, n))
  allocate(lp1(n, n))
  allocate(lp2(n, n))
  allocate(u(n))
  allocate(kai(n))
  allocate(ipiv(n))
  allocate(w1(n))
  allocate(wn(-col_num/2:col_num/2))

  do i = 1, n

     th = 2d0*pi*(dble(i)+0.5d0)/dble(n)

     x(1,i) = rad*cos(th)
     x(2,i) = rad*sin(th)

     xn(1,i) = rad*cos(th)
     xn(2,i) = rad*sin(th)

     w1(i) = sqrt(xn(1,i)**2 + xn(2,i)**2)

     xn(:,i) = xn(:,i)/w1(i)

     w1(i) = w1(i)*2d0*pi/dble(n)
  enddo

  wn(:) = w1(1)
  wn(0) = 0d0
  if(col_num .eq. 12) then
     wn(1) = wn(1)*(1d0 + 0.2051970990601252d1 + 0.2915391987686506d1)
     wn(2) = wn(2)*(1d0 - 0.7407035584542865d1 - 0.8797979464048396d1)
     wn(3) = wn(3)*(1d0 + 0.1219590847580216d2 + 0.1365562914252423d2)
     wn(4) = wn(4)*(1d0 - 0.1064623987147282d2 - 0.1157975479644601d2)
     wn(5) = wn(5)*(1d0 + 0.4799117710681772d1 + 0.5130987287355766d1)
     wn(6) = wn(6)*(1d0 - 0.8837770983721025d0 - 0.9342187797694916d0)
     wn(-1) = wn(1)
     wn(-2) = wn(2)
     wn(-3) = wn(3)
     wn(-4) = wn(4)
     wn(-5) = wn(5)
     wn(-6) = wn(6)
  else if(col_num .eq. 20) then
     wn(1) = wn(1)*(1d0 + 0.3256353919777872d1 + 0.4576078100790908d1)
     wn(2) = wn(2)*(1d0 - 0.2096116396850468d2 - 0.2469045273524281d2)
     wn(3) = wn(3)*(1d0 + 0.6872858265408605d2 + 0.7648830198138171d2)
     wn(4) = wn(4)*(1d0 - 0.1393153744796911d3 - 0.1508194558089468d3)
     wn(5) = wn(5)*(1d0 + 0.1874446431742073d3 + 0.1996415730837827d3)
     wn(6) = wn(6)*(1d0 - 0.1715855846429547d3 - 0.1807965537141134d3)
     wn(7) = wn(7)*(1d0 + 0.1061953812152787d3 + 0.1110467735366555d3)
     wn(8) = wn(8)*(1d0 - 0.4269031893958787d2 - 0.4438764193424203d2)
     wn(9) = wn(9)*(1d0 + 0.1009036069527147d2 + 0.1044548196545488d2)
     wn(10) = wn(10)*(1d0 - 0.1066655310499552d1 - 0.1100328792904271d1)
     wn(-1) = wn(1)
     wn(-2) = wn(2)
     wn(-3) = wn(3)
     wn(-4) = wn(4)
     wn(-5) = wn(5)
     wn(-6) = wn(6)
     wn(-7) = wn(7)
     wn(-8) = wn(8)
     wn(-9) = wn(9)
     wn(-10) = wn(10)
  end if

  do i = 1, n
     !x**3*y - x*y**3
     u(i) = x(1, i)**3*x(2, i) - x(1, i)*x(2, i)**3

     !(3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
     kai(i) = (3*x(1, i)**2*x(2, i) - x(2, i)**3)*xn(1, i) + (x(1, i)**3 - 3*x(1, i)*x(2, i)**2)*xn(2, i)
  end do

  exterior = 0

  do j = 1, n
     do i = 1, n
        r= sqrt(dot_product(x(:,i)-x(:,j),x(:,i)-x(:,j)))
        if(abs(j-i) .gt. n/2) then
           if(j > i) then
              tmpi = i + n
           elseif(i > j) then
              tmpi = i - n
           end if
        else
           tmpi = i
        end if

        if(r.le.1d-10) then
           lp1(i, j) = 0d0
           lp2(i, j) = 0.5d0
        elseif(abs(j - tmpi) .le. col_num/2) then
           lp1(i, j) = -wn(j-tmpi)*log(r)*arctwopi
           lp2(i, j) = -wn(j-tmpi)*dot_product(x(:,i)-x(:,j),xn(:,j))*arctwopi/r**2
        else
           lp1(i, j) = -w1(j)*log(r)*arctwopi
           lp2(i, j) = -w1(j)*dot_product(x(:,i)-x(:,j),xn(:,j))*arctwopi/r**2
        endif
     end do
  end do

  u = matmul(lp2, u)
  deallocate(lp2)

  call ZGESV(N, 1, lp1, n, ipiv, u, n, info)
  if(info.ne.0) then
     write(6,*) 'dgesv, info=', info
     call force_raise()
  endif
  deallocate(lp1)
  deallocate(ipiv)

  write(*,*) 'Laplace number of dof', n
  write(*,*) 'Laplace rad', rad
  write(*,*) 'Laplace relative error', sqrt(dot_product(u - kai, u - kai)/dot_product(kai, kai))

  deallocate(u)
  deallocate(kai)
  deallocate(x)
  deallocate(xn)
  deallocate(w1)
  deallocate(wn)

! Helmholtz ------------------------------------------------

!  write(6,*) '# n, rad'
!  open(1000, file='input', status='old')
!  read(1000,*) n, rad
!  close(1000)

  rad1 = 2d0
  rad2 = 3d0
  n = 2000
  k_1 = 3d0
  phi = 1d0
  boundary_flag = 0 !0==circle, 1==daen, 2==smooth_star
!  daen_flag = 1
!  smooth_star = 0

  allocate(x(2, n))
  allocate(xn(2, n))
  allocate(lp1(n, n))
  allocate(lp2(n, n))
  allocate(u(n))
  allocate(kai(n))
  allocate(ipiv(n))
  allocate(w1(n))
  allocate(wn(-col_num/2:col_num/2))

!  select case(smooth_star)
  select case(boundary_flag)
  case(0)
     ! circle
     do i = 1, n
        th = 2d0*pi*(dble(i)+0.5d0)/dble(n)

        x(1,i) = rad1*cos(th)
        x(2,i) = rad1*sin(th)

        xn(1,i) = rad1*cos(th)
        xn(2,i) = rad1*sin(th)

        w1(i) = sqrt(xn(1,i)**2 + xn(2,i)**2)

        xn(:,i) = xn(:,i)/w1(i)

        w1(i) = w1(i)*2d0*pi/dble(n)
     enddo
  case(1)
     ! daen circle
     dth = 2d0*pi/dble(n)
     do i = 1, n
        th = 2d0*pi*(dble(i)+0.5d0)/dble(n)

        x(1,i) = rad1*cos(th)
        x(2,i) = rad2*sin(th)

        xn(1,i) = rad2*cos(th)
        xn(2,i) = rad1*sin(th)

        xn(:,i) = xn(:,i)/sqrt(xn(1,i)**2 + xn(2,i)**2)

        w1(i) = 2d0*pi/dble(n)
     enddo
  case(2)
     ! smooth_star
     do i = 1, n
        th = 2d0*pi*(dble(i)+0.5d0)/dble(n)

        x(1,i) = (1d0 + 0.3*cos(5d0*th))*cos(th)/2.6d0
        x(2,i) = (1d0 + 0.3*cos(5d0*th))*sin(th)/2.6d0

        xn(1,i) = (cos(th) + 0.3d0*(-5d0*sin(5d0*th)*sin(th) + cos(5d0*th)*cos(th)))/2.6d0 ! dy/dth
        xn(2,i) = -(-sin(th) + 0.3d0*(-5d0*sin(5d0*th)*cos(th) - cos(5d0*th)*sin(th)))/2.6d0 ! -dx/dth

        xn(:,i) = xn(:,i)/sqrt(xn(1,i)**2 + xn(2,i)**2)

        w1(i) = 2d0*pi/dble(n)
     end do
  end select

  wn(:) = w1(1)
  wn(0) = 0d0
  if(col_num .eq. 12) then
     wn(1) = wn(1)*(1d0 + 0.2051970990601252d1 + 0.2915391987686506d1)
     wn(2) = wn(2)*(1d0 - 0.7407035584542865d1 - 0.8797979464048396d1)
     wn(3) = wn(3)*(1d0 + 0.1219590847580216d2 + 0.1365562914252423d2)
     wn(4) = wn(4)*(1d0 - 0.1064623987147282d2 - 0.1157975479644601d2)
     wn(5) = wn(5)*(1d0 + 0.4799117710681772d1 + 0.5130987287355766d1)
     wn(6) = wn(6)*(1d0 - 0.8837770983721025d0 - 0.9342187797694916d0)
     wn(-1) = wn(1)
     wn(-2) = wn(2)
     wn(-3) = wn(3)
     wn(-4) = wn(4)
     wn(-5) = wn(5)
     wn(-6) = wn(6)
  else if(col_num .eq. 20) then
     wn(1) = wn(1)*(1d0 + 0.3256353919777872d1 + 0.4576078100790908d1)
     wn(2) = wn(2)*(1d0 - 0.2096116396850468d2 - 0.2469045273524281d2)
     wn(3) = wn(3)*(1d0 + 0.6872858265408605d2 + 0.7648830198138171d2)
     wn(4) = wn(4)*(1d0 - 0.1393153744796911d3 - 0.1508194558089468d3)
     wn(5) = wn(5)*(1d0 + 0.1874446431742073d3 + 0.1996415730837827d3)
     wn(6) = wn(6)*(1d0 - 0.1715855846429547d3 - 0.1807965537141134d3)
     wn(7) = wn(7)*(1d0 + 0.1061953812152787d3 + 0.1110467735366555d3)
     wn(8) = wn(8)*(1d0 - 0.4269031893958787d2 - 0.4438764193424203d2)
     wn(9) = wn(9)*(1d0 + 0.1009036069527147d2 + 0.1044548196545488d2)
     wn(10) = wn(10)*(1d0 - 0.1066655310499552d1 - 0.1100328792904271d1)
     wn(-1) = wn(1)
     wn(-2) = wn(2)
     wn(-3) = wn(3)
     wn(-4) = wn(4)
     wn(-5) = wn(5)
     wn(-6) = wn(6)
     wn(-7) = wn(7)
     wn(-8) = wn(8)
     wn(-9) = wn(9)
     wn(-10) = wn(10)
  end if

  tmp = k_1*sin(phi)
  tmp2 = k_1*cos(phi)

  do i = 1, n
     u(i) = exp(iunit*(tmp*x(1, i)+tmp2*x(2, i)))
     kai(i) = iunit*(tmp*xn(1, i)+tmp2*xn(2, i))*u(i)

!     !x**3*y - x*y**3
!     u(i) = x(1, i)**3*x(2, i) - x(1, i)*x(2, i)**3
! 
!     !(3x**2*y - y**3)nx + (x**3 - 3xy**2)ny
!     kai(i) = (3*x(1, i)**2*x(2, i) - x(2, i)**3)*xn(1, i) + (x(1, i)**3 - 3*x(1, i)*x(2, i)**2)*xn(2, i)
  end do

  exterior = 1

  do j = 1, n
     do i = 1, n
        r= sqrt(dot_product(x(:,i)-x(:,j),x(:,i)-x(:,j)))
        if(abs(j-i) .gt. n/2) then
           if(j > i) then
              tmpi = i + n
           elseif(i > j) then
              tmpi = i - n
           end if
        else
           tmpi = i
        end if

        if(r.le.1d-10) then
           lp1(i, j) = 0d0
!           lp2(i, j) = 0.5d0
           lp2(i, j) = -0.5d0
        else
           !  call zbesh(zr,zi,fnu,kode,m,n,cyr,cyi,NZ,ierr)
           call zbesh(k_1*r, 0d0, 0d0, 1, 1, 2, CYR, CYI, NZ, IERR)
           if(ierr.ne.0) then
              write(*,*) 'nz, ierr', nz, ierr
              write(*,*) 'cyr', cyr
              write(*,*) 'cyi', cyi
              stop 'zbesh error, in predirect_helmholtz'
           end if

           besh = dcmplx(cyr, cyi)
           select case(boundary_flag)
           case(0)
              besh = (iunit*0.25d0)*besh
           case(1)
              nth = 2d0*pi*(dble(j)+0.5d0)/dble(n)
              besh = (iunit*0.25d0)*besh*sqrt((rad2*cos(nth))**2 + (rad1*sin(nth))**2)
           case(2)
              nth = 2d0*pi*(dble(j)+0.5d0)/dble(n)
              besh = (iunit*0.25d0)*besh*sqrt(((-sin(nth) + 0.3d0*(-5d0*sin(5d0*nth)*cos(nth) - cos(5d0*nth)*sin(nth)))/2.6d0)**2 + ((cos(nth) + 0.3d0*(-5d0*sin(5d0*nth)*sin(nth) + cos(5d0*nth)*cos(nth)))/2.6d0)**2)
              !        xn(1,i) = (cos(th) + 0.3d0*(-5d0*sin(5d0*th)*sin(th) + cos(5d0*th)*cos(th)))/2.6d0 ! dy/dth
              !        xn(2,i) = (-sin(th) + 0.3d0*(-5d0*sin(5d0*th)*cos(th) - cos(5d0*th)*sin(th)))/2.6d0 ! dx/dth

           end select

           if(abs(j - tmpi) .le. col_num/2) then
              lp1(i, j) = -wn(j-tmpi)*besh(1)
              lp2(i, j) = -wn(j-tmpi)*k_1*besh(2)*(dot_product(x(:,i)-x(:,j),xn(:,j))/r)
           else
              lp1(i, j) = -w1(j)*besh(1)
              lp2(i, j) = -w1(j)*k_1*besh(2)*(dot_product(x(:,i)-x(:,j),xn(:,j))/r)
           end if
        endif
     end do
  end do

  u = matmul(lp2, u)
  deallocate(lp2)

  call ZGESV(N, 1, lp1, n, ipiv, u, n, info)
  if(info.ne.0) then
     write(6,*) 'dgesv, info=', info
     call force_raise()
  endif
  deallocate(lp1)
  deallocate(ipiv)

  write(*,*) 'Helmholtz number of dof', n
  write(*,*) 'Helmholtz rad, rad2', rad, rad2
  write(*,*) 'Helmholtz relative error', sqrt(dot_product(u - kai, u - kai)/dot_product(kai, kai))

  deallocate(u)
  deallocate(kai)
  deallocate(x)
  deallocate(xn)
  deallocate(wn)

  ! singularity check--------------------------------

!  r = 0.0000001d0
!  ec = 2d0*sqrt(pi)
!  k_1 = 1d0
!  e_ite2 = 3
!  allocate(en(0:e_ite2))
!  allocate(facgam(0:e_ite2))
! 
!  write(*,*) ''
!  do j = 1, 10
!     k_1 = 1d0*j
!     write(*,*) 'now k_1', k_1
!     write(*,*) 'laplace', -arctwopi*log(r)
! 
!     call zbesh(k_1*r, 0d0, 0d0, 1, 1, 1, CYR, CYI, NZ, IERR)
!     if(ierr.ne.0) then
!        write(*,*) 'nz, ierr', nz, ierr
!        write(*,*) 'cyr', cyr
!        write(*,*) 'cyi', cyi
!        stop 'zbesh error, in predirect_helmholtz'
!     end if
!     besh(1) = dcmplx(cyr(1), cyi(1))
!     write(*,*) 'helmholtz', (iunit*0.25d0)*besh(1)
! 
!     !call DEXINT (X, N, KODE, M, TOL, EN, NZ, IERR)
!     call DEXINT(0.25d0*r**2*ec**2, 1, 1, e_ite2+1, 1d-12, en, NZ, IERR)
!     if(ierr .ne. 0) then
!        write(*,*) 'nz, ierr', nz, ierr
!        stop '# force raise !!, dexint error'
!     end if
!     do i = 0, e_ite2
!        facgam(i) = k_1**(2*i)/(gamma(dble(i+1))*ec**(2*i))
!     end do
!     write(*,*) 'ewald', 0.50d0*arctwopi*sum(facgam*en)
!     write(*,*) ''
!  end do

end program main
