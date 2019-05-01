!-------------------------------------------
module lp_lap
  implicit none
contains
  !----------------------------------------------
  complex(c_double_complex) function slp_laplace(xm,y1,y2,hy,n) bind(c)
    use iso_c_binding
    implicit none

    real(c_double),intent(in) :: xm(2),y1(2),y2(2),hy,n(2)

    real(c_double) :: t(2), r1, r2
    real(c_double) :: rr1(2), rr2(2)
    real(c_double) :: uk1,uk2,uk3,uk4,uk5
    real(c_double) :: dnrm2 !blas
    real(c_double),parameter :: arctwopi = 0.159154943091895d0

    real(c_double) :: ym(2)

    ym = (y1 + y2)*0.5d0

    if(abs(xm(1)-ym(1))+abs(xm(2)-ym(2)) .lt. 1d-8) then
       slp_laplace = hy*(1d0-log(hy*0.5d0))*arctwopi
    else
       rr1 = xm - y1
       rr2 = xm - y2
       r1 = dnrm2(2, rr1(1), 1)
       r2 = dnrm2(2, rr2(1), 1)
       t(1) = -n(2)
       t(2) = n(1)

       uk1 = dot_product(rr1,n)
       uk2 = dot_product(rr1,t)
       uk3 = dot_product(rr2,n)
       uk4 = dot_product(rr2,t)
       uk5 = atan2(uk3,uk4)-atan2(uk1,uk2)
       slp_laplace = dcmplx((uk4*log(r2)-uk2*log(r1)+hy-uk1*uk5)*arctwopi, 0d0)
    end if

  end function slp_laplace
  !----------------------------------------
  complex(c_double_complex) function dlp_laplace(xm,y1,y2,hy,n,exterior) bind(c)
    use iso_c_binding
    implicit none

    real(c_double),intent(in) :: xm(2),y1(2),y2(2),hy,n(2)
    integer(c_int),optional,intent(in) :: exterior

    real(c_double) :: t(2),r1,r2
    real(c_double) :: rr1(2), rr2(2)
    real(c_double) :: uk1,uk2,uk3,uk4,uk5
    real(c_double) :: dnrm2 !blas
    real(c_double),parameter :: arctwopi = 0.159154943091895d0

    real(c_double) :: ym(2)

    ym = (y1 + y2)*0.5d0

    if(abs(xm(1)-ym(1))+abs(xm(2)-ym(2)) .lt. 1d-8) then
       if(present(exterior)) then
          if(exterior == 0) dlp_laplace = dcmplx(0.5d0, 0d0)
          if(exterior == 1) dlp_laplace = dcmplx(-0.5d0, 0d0)
       else
          stop 'if(present(exterior)) is false in dlp_laplace'
       end if
    else
       rr1 = xm - y1
       rr2 = xm - y2
       r1 = dnrm2(2, rr1(1), 1)
       r2 = dnrm2(2, rr2(1), 1)
       t(1) = -n(2)
       t(2) = n(1)

       uk1 = dot_product(rr1, n)
       uk2 = dot_product(rr1, t)
       uk3 = dot_product(rr2, n)
       uk4 = dot_product(rr2, t)
       uk5 = atan2(uk3,uk4)-atan2(uk1,uk2)
       dlp_laplace = dcmplx(uk5*arctwopi, 0d0)
    end if

    return
  end function dlp_laplace
  !----------------------------------------
  subroutine givedata_p0_and_w(p0, w, bunten)
    implicit none
    integer,intent(in) :: bunten
    real(8) :: p0(bunten), w(bunten)

    if(bunten .eq. 2) then
       !gauss point
       p0(1) = 0.577350269189625764509148780502d0
       p0(2) = -0.577350269189625764509148780502d0

       !weight
       w(1) = 1.0
       w(2) = 1.0
    elseif(bunten .eq. 4) then
       p0(1) = 0.861136311594052575223946488893d0
       p0(2) = 0.339981043584856264802665759103d0
       p0(3) = -0.339981043584856264802665759103d0
       p0(4) = -0.861136311594052575223946488893d0

       w(1) = 0.347854845137453857373063949222d0
       w(2) = 0.65214515486254614262693605078d0
       w(3) = 0.65214515486254614262693605078d0
       w(4) = 0.347854845137453857373063949222d0
    elseif(bunten .eq. 6) then
       p0(1) = -0.932469514203152027812301554494d0
       p0(2) = -0.66120938646626451366139959502d0
       p0(3) = -0.238619186083196908630501721681d0
       p0(4) = 0.238619186083196908630501721681d0
       p0(5) = 0.66120938646626451366139959502d0
       p0(6) = 0.932469514203152027812301554494d0

       w(1) = 0.171324492379170345040296142173d0
       w(2) = 0.36076157304813860756983351384d0
       w(3) = 0.467913934572691047389870344d0
       w(4) = 0.467913934572691047389870344d0
       w(5) = 0.36076157304813860756983351384d0
       w(6) = 0.171324492379170345040296142173d0
    elseif(bunten .eq. 8) then
       p0(1) = -0.9602898564975362316836d0
       p0(2) = -0.7966664774136267395916d0
       p0(3) = -0.5255324099163289858177d0
       p0(4) = -0.1834346424956498049395d0
       p0(5) = 0.1834346424956498049395d0
       p0(6) = 0.5255324099163289858177d0
       p0(7) = 0.7966664774136267395916d0
       p0(8) = 0.9602898564975362316836d0

       w(1) = 0.1012285362903762591525d0
       w(2) = 0.2223810344533744705444d0
       w(3) = 0.313706645877887287338d0
       w(4) = 0.3626837833783619829652d0
       w(5) = 0.3626837833783619829652d0
       w(6) = 0.313706645877887287338d0
       w(7) = 0.222381034453374470544d0
       w(8) = 0.1012285362903762591525d0
    elseif(bunten .eq. 32) then
       p0(1) = -0.997263861849481563545d0
       p0(2) = -0.9856115115452683354002d0
       p0(3) = -0.9647622555875064307738d0
       p0(4) = -0.9349060759377396891709d0
       p0(5) = -0.8963211557660521239653d0
       p0(6) = -0.8493676137325699701337d0
       p0(7) = -0.794483795967942406963d0
       p0(8) = -0.7321821187402896803874d0
       p0(9) = -0.6630442669302152009751d0
       p0(10) = -0.5877157572407623290408d0
       p0(11) = -0.5068999089322293900238d0
       p0(12) = -0.421351276130635345364d0
       p0(13) = -0.3318686022821276497799d0
       p0(14) = -0.2392873622521370745446d0
       p0(15) = -0.1444719615827964934852d0
       p0(16) = -0.0483076656877383162348d0
       p0(17) = 0.048307665687738316235d0
       p0(18) = 0.1444719615827964934852d0
       p0(19) = 0.2392873622521370745446d0
       p0(20) = 0.33186860228212764978d0
       p0(21) = 0.4213512761306353453641d0
       p0(22) = 0.5068999089322293900238d0
       p0(23) = 0.5877157572407623290408d0
       p0(24) = 0.6630442669302152009751d0
       p0(25) = 0.7321821187402896803874d0
       p0(26) = 0.7944837959679424069631d0
       p0(27) = 0.8493676137325699701337d0
       p0(28) = 0.8963211557660521239653d0
       p0(29) = 0.9349060759377396891709d0
       p0(30) = 0.9647622555875064307738d0
       p0(31) = 0.9856115115452683354002d0
       p0(32) = 0.997263861849481563545d0

       w(1) = 0.0070186100094700966004d0
       w(2) = 0.0162743947309056706052d0
       w(3) = 0.0253920653092620594558d0
       w(4) = 0.0342738629130214331027d0
       w(5) = 0.0428358980222266806569d0
       w(6) = 0.050998059262376176196d0
       w(7) = 0.0586840934785355471453d0
       w(8) = 0.065822222776361846838d0
       w(9) = 0.072345794108848506225d0
       w(10) = 0.0781938957870703064717d0
       w(11) = 0.0833119242269467552222d0
       w(12) = 0.087652093004403811143d0
       w(13) = 0.091173878695763884713d0
       w(14) = 0.09384439908080456563918d0
       w(15) = 0.0956387200792748594191d0
       w(16) = 0.0965400885147278005668d0
       w(17) = 0.0965400885147278005668d0
       w(18) = 0.0956387200792748594191d0
       w(19) = 0.0938443990808045656392d0
       w(20) = 0.091173878695763884713d0
       w(21) = 0.0876520930044038111428d0
       w(22) = 0.083311924226946755222d0
       w(23) = 0.078193895787070306472d0
       w(24) = 0.072345794108848506225d0
       w(25) = 0.065822222776361846838d0
       w(26) = 0.0586840934785355471453d0
       w(27) = 0.0509980592623761761962d0
       w(28) = 0.0428358980222266806569d0
       w(29) = 0.0342738629130214331027d0
       w(30) = 0.0253920653092620594558d0
       w(31) = 0.0162743947309056706052d0
       w(32) = 0.0070186100094700966004d0
    elseif(bunten .eq. 100) then
       p0(1) = -0.9997137267734412336782d0
       p0(2) = -0.9984919506395958184d0
       p0(3) = -0.996295134733125149186d0
       p0(4) = -0.993124937037443459652d0
       p0(5) = -0.9889843952429917480044d0
       p0(6) = -0.9838775407060570154961d0
       p0(7) = -0.9778093584869182885538d0
       p0(8) = -0.9707857757637063319309d0
       p0(9) = -0.9628136542558155272937d0
       p0(10) = -0.9539007829254917428493d0
       p0(11) = -0.9440558701362559779628d0
       p0(12) = -0.9332885350430795459243d0
       p0(13) = -0.921609298145333952667d0
       p0(14) = -0.9090295709825296904671d0
       p0(15) = -0.8955616449707269866985d0
       p0(16) = -0.8812186793850184155733d0
       p0(17) = -0.8660146884971646234107d0
       p0(18) = -0.8499645278795912842934d0
       p0(19) = -0.8330838798884008235429d0
       p0(20) = -0.815389238339176254394d0
       p0(21) = -0.7968978923903144763896d0
       p0(22) = -0.777627909649495475628d0
       p0(23) = -0.757598118519707176036d0
       p0(24) = -0.7368280898020207055124d0
       p0(25) = -0.71533811757305644646d0
       p0(26) = -0.6931491993558019659487d0
       p0(27) = -0.670283015603141015803d0
       p0(28) = -0.6467619085141292798326d0
       p0(29) = -0.622608860203707771604d0
       p0(30) = -0.5978474702471787212648d0
       p0(31) = -0.5725019326213811913169d0
       p0(32) = -0.546597012065094167468d0
       p0(33) = -0.520158019881763056647d0
       p0(34) = -0.493210789208190933569d0
       p0(35) = -0.465781649773358042249d0
       p0(36) = -0.437897402172031513109d0
       p0(37) = -0.4095852916783015425289d0
       p0(38) = -0.3808729816246299567634d0
       p0(39) = -0.3517885263724217209723d0
       p0(40) = -0.3223603439005291517225d0
       p0(41) = -0.2926171880384719647376d0
       p0(42) = -0.2625881203715034791689d0
       p0(43) = -0.2323024818449739696495d0
       p0(44) = -0.2017898640957359972361d0
       p0(45) = -0.1710800805386032748875d0
       p0(46) = -0.1402031372361139732075d0
       p0(47) = -0.1091892035800611150034d0
       p0(48) = -0.0780685828134366366948d0
       p0(49) = -0.046871682421591631615d0
       p0(50) = -0.015628984421543082872d0
       p0(51) = 0.0156289844215430828722d0
       p0(52) = 0.0468716824215916316149d0
       p0(53) = 0.0780685828134366366948d0
       p0(54) = 0.1091892035800611150034d0
       p0(55) = 0.140203137236113973208d0
       p0(56) = 0.1710800805386032748875d0
       p0(57) = 0.2017898640957359972361d0
       p0(58) = 0.23230248184497396965d0
       p0(59) = 0.262588120371503479169d0
       p0(60) = 0.292617188038471964738d0
       p0(61) = 0.3223603439005291517225d0
       p0(62) = 0.351788526372421720972d0
       p0(63) = 0.3808729816246299567634d0
       p0(64) = 0.4095852916783015425289d0
       p0(65) = 0.437897402172031513109d0
       p0(66) = 0.4657816497733580422492d0
       p0(67) = 0.4932107892081909335693d0
       p0(68) = 0.5201580198817630566468d0
       p0(69) = 0.546597012065094167468d0
       p0(70) = 0.5725019326213811913169d0
       p0(71) = 0.597847470247178721265d0
       p0(72) = 0.6226088602037077716042d0
       p0(73) = 0.6467619085141292798326d0
       p0(74) = 0.6702830156031410158026d0
       p0(75) = 0.6931491993558019659487d0
       p0(76) = 0.71533811757305644646d0
       p0(77) = 0.7368280898020207055124d0
       p0(78) = 0.7575981185197071760357d0
       p0(79) = 0.7776279096494954756276d0
       p0(80) = 0.7968978923903144763896d0
       p0(81) = 0.815389238339176254394d0
       p0(82) = 0.8330838798884008235429d0
       p0(83) = 0.8499645278795912842934d0
       p0(84) = 0.8660146884971646234107d0
       p0(85) = 0.8812186793850184155733d0
       p0(86) = 0.895561644970726986699d0
       p0(87) = 0.9090295709825296904671d0
       p0(88) = 0.921609298145333952667d0
       p0(89) = 0.9332885350430795459243d0
       p0(90) = 0.9440558701362559779628d0
       p0(91) = 0.953900782925491742849d0
       p0(92) = 0.9628136542558155272937d0
       p0(93) = 0.9707857757637063319309d0
       p0(94) = 0.9778093584869182885538d0
       p0(95) = 0.9838775407060570154961d0
       p0(96) = 0.9889843952429917480044d0
       p0(97) = 0.993124937037443459652d0
       p0(98) = 0.9962951347331251491861d0
       p0(99) = 0.9984919506395958184002d0
       p0(100) = 0.9997137267734412336782d0

       w(1) = 7.346344905056717304d-4
       w(2) = 0.00170939265351810524d0
       w(3) = 0.0026839253715534824194d0
       w(4) = 0.0036559612013263751823d0
       w(5) = 0.0046244500634221193511d0
       w(6) = 0.005588428003865515157d0
       w(7) = 0.006546948450845322764d0
       w(8) = 0.007499073255464711579d0
       w(9) = 0.008443871469668971403d0
       w(10) = 0.009380419653694457951418d0
       w(11) = 0.0103078025748689695858d0
       w(12) = 0.011225114023185977117d0
       w(13) = 0.0121314576629794974077d0
       w(14) = 0.0130259478929715422856d0
       w(15) = 0.013907710703718772688d0
       w(16) = 0.0147758845274413017689d0
       w(17) = 0.0156296210775460027239d0
       w(18) = 0.0164680861761452126431d0
       w(19) = 0.01729046056832358243934d0
       w(20) = 0.0180959407221281166644d0
       w(21) = 0.0188837396133749045529d0
       w(22) = 0.0196530874944353058654d0
       w(23) = 0.0204032326462094327668d0
       w(24) = 0.021133442112527641543d0
       w(25) = 0.021843002416247386314d0
       w(26) = 0.0225312202563362727018d0
       w(27) = 0.0231974231852541216225d0
       w(28) = 0.0238409602659682059626d0
       w(29) = 0.02446120270795705272d0
       w(30) = 0.02505754448157958970376d0
       w(31) = 0.025629402910208116076d0
       w(32) = 0.026176219239545676342d0
       w(33) = 0.02669745918357096266d0
       w(34) = 0.0271926134465768801365d0
       w(35) = 0.0276611982207923882942d0
       w(36) = 0.0281027556591011733176d0
       w(37) = 0.0285168543223950979909d0
       w(38) = 0.0289030896011252031349d0
       w(39) = 0.0292610841106382766201d0
       w(40) = 0.0295904880599126425118d0
       w(41) = 0.0298909795933328309168d0
       w(42) = 0.03016226510516914491907d0
       w(43) = 0.03040407952645482001651d0
       w(44) = 0.0306161865839804484965d0
       w(45) = 0.0307983790311525904277d0
       w(46) = 0.030950478850490988234d0
       w(47) = 0.031072337427566516588d0
       w(48) = 0.0311638356962099067838d0
       w(49) = 0.0312248842548493577324d0
       w(50) = 0.0312554234538633569476d0
       w(51) = 0.0312554234538633569476d0
       w(52) = 0.0312248842548493577324d0
       w(53) = 0.031163835696209906784d0
       w(54) = 0.031072337427566516588d0
       w(55) = 0.0309504788504909882341d0
       w(56) = 0.030798379031152590428d0
       w(57) = 0.0306161865839804484965d0
       w(58) = 0.0304040795264548200165d0
       w(59) = 0.030162265105169144919d0
       w(60) = 0.02989097959333283091684d0
       w(61) = 0.029590488059912642512d0
       w(62) = 0.02926108411063827662d0
       w(63) = 0.02890308960112520313488d0
       w(64) = 0.0285168543223950979909d0
       w(65) = 0.0281027556591011733176d0
       w(66) = 0.027661198220792388294d0
       w(67) = 0.0271926134465768801365d0
       w(68) = 0.02669745918357096266d0
       w(69) = 0.026176219239545676342d0
       w(70) = 0.025629402910208116076d0
       w(71) = 0.025057544481579589704d0
       w(72) = 0.02446120270795705272d0
       w(73) = 0.0238409602659682059626d0
       w(74) = 0.023197423185254121622d0
       w(75) = 0.0225312202563362727018d0
       w(76) = 0.021843002416247386314d0
       w(77) = 0.02113344211252764154267d0
       w(78) = 0.020403232646209432767d0
       w(79) = 0.0196530874944353058654d0
       w(80) = 0.0188837396133749045529d0
       w(81) = 0.018095940722128116664d0
       w(82) = 0.0172904605683235824393d0
       w(83) = 0.016468086176145212643d0
       w(84) = 0.0156296210775460027239d0
       w(85) = 0.0147758845274413017689d0
       w(86) = 0.013907710703718772688d0
       w(87) = 0.013025947892971542286d0
       w(88) = 0.0121314576629794974077d0
       w(89) = 0.0112251140231859771172d0
       w(90) = 0.010307802574868969586d0
       w(91) = 0.009380419653694457951d0
       w(92) = 0.008443871469668971402621d0
       w(93) = 0.00749907325546471157883d0
       w(94) = 0.0065469484508453227642d0
       w(95) = 0.0055884280038655151572d0
       w(96) = 0.0046244500634221193511d0
       w(97) = 0.0036559612013263751823d0
       w(98) = 0.0026839253715534824194d0
       w(99) = 0.00170939265351810524d0
       w(100) = 7.3463449050567173d-4
    end if

  end subroutine givedata_p0_and_w
  !--------------------------------------------
  function make_sol_mie_series_nystrom(k_1, x, rad, slp_or_dlp) result(res)
    implicit none
    !order of hankel and bessel function is equal to one

    complex*16,intent(in) :: k_1
    real(8),intent(in) :: x(:,:)
    real(8),intent(in) :: rad
    integer,intent(in) :: slp_or_dlp

    complex*16,allocatable :: res(:)
    integer :: len
    real(8) :: xm(2)
    integer,parameter :: n = 3
    real(8) :: cyr(n), cyi(n), zr, zi
    real(8) :: fnu
    integer :: kode, nz, ierr, m
    integer :: i
    integer :: ii
    complex*16 :: hankel(3)
    complex*16 :: bessel(3)
    complex*16 :: am
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    len = size(x(1,:))

    allocate(res(len))
    res = dcmplx(0d0, 0d0)

    zr = dble(k_1*rad)
    zi = dimag(k_1*rad) ! k_1 is real
    fnu = 0
    kode = 1
    call ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
    bessel(:) = dcmplx(cyr(:), cyi(:))
    m = 1
    call ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
    hankel(:) = dcmplx(cyr(:), cyi(:))

    do i = 1, len
       if(slp_or_dlp .eq. 1 .or. slp_or_dlp .eq. 3) then
          res(i) = -bessel(2)/hankel(2)*0.5d0*k_1*(hankel(1) - hankel(3))*exp(iunit*1d0*atan2(x(2, i), x(1, i)))&
               & + 0.5d0*k_1*(bessel(1) - bessel(3))*exp(iunit*1d0*atan2(x(2, i), x(1, i)))
       else if(slp_or_dlp .eq. 2 .or. slp_or_dlp .eq. 4) then
          am = -(bessel(1) - bessel(3))/(hankel(1) - hankel(3))
          res(i) = (am*hankel(2) + bessel(2))*exp(iunit*1d0*atan2(x(2, i), x(1, i)))
       end if
    end do

  end function make_sol_mie_series_nystrom
  !---------------------------------------
  function d_diff_kernel_hel(xx, yy, xn, yn, k_1) result(res)
    implicit none

    real(8),intent(in) :: xx(2)
    real(8),intent(in) :: yy(2)
    real(8),intent(in) :: xn(2)
    real(8),intent(in) :: yn(2)
    complex*16,intent(in) :: k_1

    complex*16 :: res
    real(8) :: r
    real(8) :: cyr(3)
    real(8) :: cyi(3)
    integer :: ierr
    integer :: nz
    complex*16 :: besh(3)
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)

    r = sqrt(dot_product(xx - yy, xx - yy))
    call zbesh(dble(k_1*r), 0d0, 0d0, 1, 1, 3, CYR, CYI, NZ, IERR)
    if(ierr.ne.0) then
       write(*,*) 'nz, ierr', nz, ierr
       write(*,*) 'cyr', cyr
       write(*,*) 'cyi', cyi
       stop 'zbesh error'
    end if

    besh = dcmplx(cyr, cyi)
    besh = (iunit*0.25d0)*besh

!    block
!      integer :: i
!      do i = 1, 3
!         write(*,*) 'i, besh', i, besh(i)
!      end do
!    end block

    res = k_1*besh(2)*dot_product(xn, yn)/r + (dot_product(xx-yy,yn)/r)*(dot_product(xx-yy,xn)/r)*(-1d0/r*k_1*besh(2) + 0.5d0*k_1**2*(besh(1) - besh(3)))

  end function d_diff_kernel_hel
  !----------------------------------------
  function lp_hel(xm, ym, y1, y2, hy, nx, ny, k, slp_or_dlp, msw_flag, exterior) result(result)
!    use my_slatec_func
!    use several_func
    implicit none

    real(8),intent(in) :: xm(2)
    real(8),intent(in) :: ym(2)
    real(8),intent(in) :: y1(2)
    real(8),intent(in) :: y2(2)
    real(8),intent(in) :: hy
    real(8),intent(in) :: nx(2)
    real(8),intent(in) :: ny(2)
    complex*16,intent(in) :: k
    integer,intent(in) :: slp_or_dlp
    integer,intent(in) :: msw_flag
    integer,intent(in) :: exterior

    integer,parameter :: bunten = 4
    real(8),parameter :: arctwopi = 0.159154943091895d0
    complex*16,parameter :: iunit = dcmplx(0d0, 1d0)
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
    real(8) :: zr, zi
    real(8) :: cyr, cyi
    real(8) :: fnu
    integer :: kode
    integer :: mmm, nnn
    integer :: ierr, nz

    real(8) :: s, r
    real(8) :: laplace
    real(8) :: RR(2)
    complex*16 :: besh
    real(8) :: p0(bunten), w(bunten)
    integer :: i, j

    result = dcmplx(0d0, 0d0)

    ! near
    if(abs(xm(1)-ym(1))+abs(xm(2)-ym(2)) .lt. 10d-8) then
       select case(slp_or_dlp)
       case(1, 4)
          result = hy*(1d0-log(hy*0.5d0))*arctwopi
       case(2, 3)
          if(exterior == 0) then
             if(slp_or_dlp .eq. 2) then
                result = dcmplx(0.5d0, 0d0)
             else !if(slp_or_dlp .eq. 3)
                result = dcmplx(-0.5d0, 0d0)
             end if
          elseif(exterior == 1) then
             if(slp_or_dlp .eq. 2) then
                result = dcmplx(-0.5d0, 0d0)
             else !if(slp_or_dlp .eq. 3)
                result = dcmplx(0.5d0, 0d0)
             end if
          end if
          return
       end select
    else
       select case(slp_or_dlp)
       case(1, 4)
          result = slp_laplace(xm,y1,y2,hy,ny)
       case(2, 3)
          ! result = dcmplx(0d0, 0d0)
       end select
    end if

    ! far
    call givedata_p0_and_w(p0, w, bunten)

    do i=1,bunten
       s = 0.5d0 + 0.5d0*p0(i)
       RR = y1 + (y2 - y1) * s
       RR = xm - RR
       r = dnrm2(2, RR(1), 1)

       zr = dble(k*r)
       zi = dimag(k*r)
       kode = 1
       if(msw_flag .eq. 0) then
          mmm = 1
       elseif(msw_flag .eq. 1) then
          mmm = 2
       end if

       select case(slp_or_dlp)
       case(1, 4)
          fnu = 0d0
       case(2, 3)
          fnu = 1d0
       end select

!       if(zi .eq. 0d0) then
!          call my_zbesh(zr,zi,fnu,kode,mmm,1,cyr,cyi,NZ,ierr)
!       else if(zi .ne. 0d0) then
          call zbesh(zr,zi,fnu,kode,mmm,1,cyr,cyi,NZ,ierr)
          if(ierr.ne.0) then
             !$omp critical
             write(*,*) 'nz, ierr', nz, ierr
             write(*,*) 'zr, zi', zr, zi
             write(*,*) 'cyr', cyr
             write(*,*) 'cyi', cyi
             stop 'zbesh error, in predirect_helmholtz'
             !$omp end critical
          end if
!       end if

       besh = dcmplx(cyr, cyi)
       if(msw_flag .eq. 0) then
          besh = (iunit*0.25d0)*besh
       elseif(msw_flag .eq. 1) then
          besh = -(iunit*0.25d0)*besh
       end if

       select case(slp_or_dlp)
       case(1, 4)
          laplace = -log(r)*arctwopi
          result = result + w(i)*(besh-laplace)*hy*0.5d0
       case(2)
          laplace = 0d0
          result = result + w(i)*(k*besh*(dot_product(RR, ny)/r)-laplace)*hy*0.5d0
       case(3)
          laplace = 0d0
          result = result - w(i)*(k*besh*(dot_product(RR, nx)/r)-laplace)*hy*0.5d0
       end select
    end do

    !d_dlp
    select case(slp_or_dlp)
    case(1, 2, 3)
       return
    case(4)
       dotn = dot_product(nx, ny)
       result = dotn*k**2*result

       rr1 = xm - y1
       rr2 = xm - y2
       r1 = dnrm2(2, rr1(1), 1)
       r2 = dnrm2(2, rr2(1), 1)
       sx(1) = -nx(2)
       sx(2) = nx(1)
       doty1 = dot_product(rr1, sx)
       doty2 = dot_product(rr2, sx)

       if(msw_flag .eq. 0) then
          mmm = 1
       elseif(msw_flag .eq. 1) then
          mmm = 2
       end if

       zr = dble(k*r1)
       zi = dimag(k*r1)
!       if(zi .eq. 0d0) then
!          call my_zbesh(zr, zi, 1d0, 1, mmm, 1, cyr, cyi, nz, ierr)
!       elseif(zi .ne. 0d0) then
          call zbesh(zr, zi, 1d0, 1, mmm, 1, cyr, cyi, nz, ierr)
!       end if
       if(ierr .ne. 0) stop 'zbesh error, in subroutine d_dlp_helmloltz of predirect_helmholtz'
       beshy1 = dcmplx(cyr, cyi)

       zr = dble(k*r2)
       zi = dimag(k*r2)
!       if(zi .eq. 0d0) then
!          call my_zbesh(zr, zi, 1d0, 1, mmm, 1, cyr, cyi, nz, ierr)
!       elseif(zi .ne. 0d0) then
          call zbesh(zr, zi, 1d0, 1, mmm, 1, cyr, cyi, nz, ierr)
!       end if
       if(ierr .ne. 0) stop 'zbesh error, in subroutine d_dlp_helmloltz of predirect_helmholtz'
       beshy2 = dcmplx(cyr, cyi)

       if(msw_flag .eq. 0) then
          result = result - 0.25d0*iunit*k*(beshy1*doty1/r1 - beshy2*doty2/r2)
       elseif(msw_flag .eq. 1) then
          result = result + 0.25d0*iunit*k*(beshy1*doty1/r1 - beshy2*doty2/r2)
       end if
    end select

  end function lp_hel
!--------------------------------------------
  subroutine rank_check_svd(a)
!    use mod_globals
    implicit none

    complex*16, intent(in) :: a(:,:)

    !    .. Scalar Arguments ..
    !    CHARACTER          JOBU, JOBVT
    INTEGER :: INFO, LDA, LDU, LDVT, LWORK, M, N
    !     .. Array Arguments ..
    DOUBLE PRECISION,allocatable :: RWORK(:), S(:)
    COMPLEX*16,allocatable :: dummy(:,:), U(:,:), VT(:,:),  WORK(:)

    integer :: i

    m = size(a(:,1))
    n = size(a(1,:))

    lda = m
    ldu = 1
    ldvt = 1
    lwork = max(1, 2*min(m, n) + max(m, n)) + 1

    allocate(dummy(m, n))
    allocate(s(min(m, n)))
    allocate(u(1,1))
    allocate(vt(1,1))
    allocate(work(lwork))
    allocate(rwork(5*min(m, n)))

    dummy = a

    call ZGESVD('N', 'N', M, N, dummy, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, INFO)

    do i = 1, min(m, n)
       write(*,*) 'singular value, i', s(i), i
    end do
    write(*,*) 'condition number', s(1)/s(min(m, n))

  end subroutine rank_check_svd
  !--------------------------------------------------------------------
end module lp_lap
