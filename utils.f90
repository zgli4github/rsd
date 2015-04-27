!include 'integral.f90'

module const
!====================================================
! Some constants: physical and cosmological
!====================================================
  implicit none
  real(kind=4),parameter :: pi = 3.14159265358979323846
  real(kind=4),parameter :: pih = 1.57079632679489661923
  real(kind=4),parameter :: twopi = 6.28318530717958647692
  real(kind=4),parameter :: r2d = 57.2957795130823208767981
  real(kind=4),parameter :: d2r = 0.0174532925199432957692 

  real(kind=4),parameter :: Mpc2cm = 3.08425188e24 
  real(kind=4),parameter :: Year2sec = 3.155815e7
  real(kind=4),parameter :: Jansky = 1e-30                   ! [W cm^-2 Hz^-1]
  real(kind=4),parameter :: KeV2K = 1.1605e7 
  real(kind=4),parameter :: SolarMass = 1.989e33             ! [g]
  real(kind=4),parameter :: CriticalDensityhm2=1.8804448e-29 ! [h^2 g cm^-3], (h is Hubble parameter)

  real(kind=4),parameter :: LightSpeed = 2.99792458e5        ! [km s^-1]
  real(kind=4),parameter :: Alpha_inv = 137.03599976         ! (= 1/alpha)
  real(kind=4),parameter :: ThomsonCS = 6.65245854e-25       ! [cm^2]
  real(kind=4),parameter :: NewtonG = 6.67384e-8             ! [cm^3 g^-1 s^-2]
  real(kind=4),parameter :: Planck = 6.58211889e-16          ! [eV*s], (= hbar or reduced Planck constant)
  real(kind=4),parameter :: Planck_cmgs = 1.0546e-27         ! [cm^2 g s^-1], (= hbar or reduced Planck constant)
  real(kind=4),parameter :: Boltzmann = 8.617342e-5          ! [eV K^-1]
  real(kind=4),parameter :: ElectronCharge = 1.602176565e-19 ! [C]
  real(kind=4),parameter :: ProtonMass = 1.672621777e-24     ! [g]
  real(kind=4),parameter :: ElectronMass = 9.109389e-28      ! [g]
  real(kind=4),parameter :: NeutronMass = 1.674927351e-24    ! [g]
  real(kind=4),parameter :: MuonMass = 0.188353147e-24       ! [g]
  real(kind=4),parameter :: BohrRadius = 0.52917721092e-8    ! [cm]
  real(kind=4),parameter :: ProtonRadius = 0.875e-13         ! [cm]

  real(kind=4),parameter :: MassofEnergy = 1.782661845e-33   ! [g/eV], (= eV/(g*c^2), Energy to Mass.)
  real(kind=4),parameter :: EnergyofMass = 0.5609588845e33   ! [eV/g], (= g*c^2/eV, Mass to energy.)
  real(kind=4),parameter :: eV2GHz = 2.4179893e5             ! Frequence to energy.
  real(kind=4),parameter :: GHz2eV = 4.1356676e-6            ! Energy ro frequency.
    
  real(kind=4),parameter :: badvalue=-1.63750E+30


  !-----------------------------
  ! Double precision counterpart
  !-----------------------------
  real(kind=8),parameter :: pid = 3.14159265358979323846d0
  real(kind=8),parameter :: pihd =  1.57079632679489661923d0
  real(kind=8),parameter :: twopid = 6.28318530717958647692d0
  real(kind=8),parameter :: r2dd = 57.2957795130823208767981d0 
  real(kind=8),parameter :: d2rd = 0.0174532925199432957692d0 

  real(kind=8),parameter :: Mpc2cmd = 3.0842518868d24 
  real(kind=8),parameter :: Year2secd = 3.155815d7
  real(kind=8),parameter :: Janskyd = 1d-30
  real(kind=8),parameter :: KeV2Kd = 1.1605d0
  real(kind=8),parameter :: SolarMassd = 1.989d33
  real(kind=8),parameter :: CriticalDensityhm2d = 1.8804448d-29 

  real(kind=8),parameter :: LightSpeedd = 2.99792458d5  
  real(kind=8),parameter :: Alpha_invd = 137.03599976d0 
  real(kind=8),parameter :: ThomsonCSd = 6.65245854d-25
  real(kind=8),parameter :: NewtonGd = 6.67384d-8 
  real(kind=8),parameter :: Planckd = 6.58211889d-16
  real(kind=4),parameter :: Planckd_cmgs = 1.0546d-27
  real(kind=8),parameter :: Boltzmannd = 8.617342d-5
  real(kind=8),parameter :: ProtonMassd = 1.672621777d-24
  real(kind=8),parameter :: ElectronMassd = 9.10938291d-28
  real(kind=8),parameter :: NeutronMassd = 1.674927351d-24
  real(kind=8),parameter :: MuonMassd = 0.188353147d-24 
  real(kind=8),parameter :: BohrRadiusd = 0.52917721092d-8
  real(kind=8),parameter :: ProtonRadiusd = 0.875d-13
  real(kind=8),parameter :: ElectronCharged = 1.602176565d-19
  real(kind=8),parameter :: MassofEnergyd = 1.782661845d-33
  real(kind=8),parameter :: EnergyofMassd = 0.5609588845d33
  real(kind=8),parameter :: eV2GHzd = 2.4179893d5
  real(kind=8),parameter :: GHz2eVd = 4.1356676d-6
  real(kind=8),parameter :: badvalued=-1.63750d30

end module const


module coordtransform
  use const
  implicit none
  contains
  !======================================================================
  ! Coordinate tranformations (all in degrees)
  !======================================================================
  ! SDSS <==> J2000: "sdss2j2k(lamda,eta,alpha,delta)"
  !                  "j2k2sdss(alpha,delta,lamda,eta)"
  !----------------------------------------------------------------------
  ! Equator <==> Galactic: "equator2galactic(ra,dec,l,b)"
  !                        "galactic2equator(l,b,ra,dec)"
  !----------------------------------------------------------------------
  ! Euclide <==> Shperical: "Euclide2Shpere(x,y,z,r,ra,dec)"
  !                         "Sphere2Euclide(r,ra,dec,x,y,z)"
  !----------------------------------------------------------------------
  ! Unit vector <==> position on sphere: "angle2vectro(ra,dec,v)"
  !                                      "vector2angle(v,ra,dec)"
  ! ==== v is a 3-vector.
  !----------------------------------------------------------------------
  ! Rotation: "rotation_zxz(vin,vout,alpha,beta,gama)"
  !      1) rotate vector by alpha about z-axis;
  !      2)          then by beta about x-axis;
  !      3)          then by gama about z-axis.
  ! Including 3 subroutines: "rotation_x(vin,vout,alpha)"
  !                          "rotation_y(vin,vout,alpha)"
  !                          "rotation_z(vin,vout,alpha)"
  ! ==== vin and vout are 3-vectors.
  !======================================================================

  subroutine b19502j2k(ra,dec,raj2k,decj2k)
    implicit none
    real*8,intent(in) :: ra,dec
    real*8,intent(out) :: raj2k,decj2k
    real*8 :: a = 0.640265d0
    real*8 :: b = 0.278369d0

    raj2k = ra + a + b*sin(ra*d2rd)*tan(dec*d2rd)
    decj2k = dec + b*cos(ra*d2rd)
    return
  end subroutine

  subroutine sdss2j2k(lambda,eta,alpha,delta)
    implicit none
    real*8,intent(in) :: lambda,eta
    real*8,intent(out) :: alpha,delta
    real*8 :: calpha,salpha,cdelta,sdelta,alpha0

    sdelta = cos(lambda*d2rd)*sin((eta+32.5d0)*d2rd)
    delta = r2dd*asin(sdelta)
    cdelta = cos(delta*d2rd)
    if((delta-90.d0).ge.-1.d-8) then
      delta = 90.d0
      alpha = 0d0
    else if((delta+90.d0).lt.1.d-8) then
      delta = -90.d0
      alpha = 0d0
    else
      salpha = cos(lambda*d2rd)*cos((eta+32.5d0)*d2rd)/cdelta
      if(salpha.ge.1d0) salpha = 1d0    ! avoid read-out error
      if(salpha.le.-1d0) salpha = -1d0  ! avoid read-out error
      alpha = r2dd*asin(salpha)
      calpha = -sin(lambda*d2rd)/cdelta
      if(salpha.ge.0d0) then
        if(calpha.ge.0d0) then
          alpha = alpha
        else
          alpha = 180d0 - alpha
        end if
      else
        if(calpha.ge.0d0) then
          alpha = alpha + 360d0
        else
          alpha = 180d0 - alpha
        end if
      end if
      alpha = alpha + 95d0
      if(alpha.ge.360d0) alpha = alpha - 360d0
    end if
    return
  end subroutine sdss2j2k

  subroutine j2k2sdss(alpha,delta,lamda,eta)
    implicit none
    real*8,intent(in) :: alpha,delta
    real*8,intent(out) :: lamda,eta
    real*8 :: slam,clam,seta325,ceta325,a,b,cdelta,sdelta,tgeta325

    cdelta = cos(delta*d2rd)
    sdelta = sin(delta*d2rd)
    slam = - cos((alpha - 95d0)*d2rd)*cos(delta*d2rd)
    if((1d0-slam).lt.1d-8) then
      lamda = 90d0 
      eta = 0d0
    else if((1d0+slam).lt.1d-8) then
      lamda = -90d0
      eta = 0d0
    else if((1d0+sdelta).lt.1d-8) then
      eta = 57.5d0
      lamda = 180d0
    else if((1d0-sdelta).lt.1d-8) then
      eta = 57.5d0
      lamda = 0d0
    else
      tgeta325 = tan(delta*d2rd)/sin((alpha-95d0)*d2rd)
      a = tan((32.5d0 - 90d0)*d2rd)
      if(tgeta325.gt.a) then
        eta = atan(tgeta325)*r2dd - 32.5d0
      else 
        eta = atan(tgeta325)*r2dd + 180d0 - 32.5d0
      end if
      clam = sin(delta*d2rd)/sin((eta+32.5d0)*d2rd)
      if(clam.ge.0d0) then
        lamda = asin(slam)*r2dd
      else
        if(slam.ge.0d0) then
          lamda = 180d0 - asin(slam)*r2dd
        else
          lamda = -180d0 - asin(slam)*r2dd
        end if
      end if 
    end if
  end subroutine j2k2sdss
          
  subroutine equator2galactic(ra,dec,l,b)
     implicit none
     real*8 :: ra,dec,l,b,rar,decr,lr,br
     real*8 :: sinra,cosra,sindec,cosdec,sinl,cosl,sinb,cosb
     real*8,parameter :: ragp=192.859508d0
     real*8,parameter :: decgp=27.128336d0
     real*8,parameter :: lcp = 122.932d0
     real*8 :: ragpr,decgpr,lcpr
     real*8 :: ramgp,lcpml
  
     ragpr=ragp*d2rd
     decgpr=decgp*d2rd
     lcpr=lcp*d2rd
     rar=ra*d2rd
     decr=dec*d2rd
     lr=l*d2rd
     br=b*d2rd
     ramgp=rar-ragpr
     lcpml=lcpr-lr
     
     sinb=sin(decgpr)*sin(decr)+cos(decgpr)*cos(decr)*cos(ramgp)
     if(sinb.gt.1d0) sinb=1d0
     if(sinb.lt.-1d0) sinb=-1d0
     br=asin(sinb)
     b=br*r2dd
     cosb=cos(br)
     if(cosb.eq.0d0) then
        write(*,*) "Pole"
        l=0d0
     else
        sinl=cos(decr)*sin(ramgp)/cosb
        cosl=(cos(decgpr)*sin(decr)-sin(decgpr)*cos(decr)*cos(ramgp))/cosb
        if(sinl.gt.1d0) sinl=1d0
        if(sinl.lt.-1d0) sinl=-1d0
        if(sinl.ge.0d0) then
          if(cosl.ge.0d0) then
            l=lcp-asin(sinl)*r2dd
          else
            l=lcp-180d0+asin(sinl)*r2dd
          end if
        else 
          if(cosl.ge.0d0) then
            l=lcp-360d0-asin(sinl)*r2dd
          else
            l=lcp-180d0+asin(sinl)*r2dd
          end if
        end if
     end if
     if(l.ge.360d0) l=l-360d0
     if(l.lt.0d0) l=l+360d0
     return
  end subroutine equator2galactic
  
  subroutine galactic2equator(l,b,ra,dec)
     implicit none
     real*8 :: ra,dec,l,b,rar,decr,lr,br
     real*8 :: sinra,cosra,sindec,cosdec,sinl,cosl,sinb,cosb
     real*8,parameter :: ragp=192.859508d0
     real*8,parameter :: decgp=27.128336d0
     real*8,parameter :: lcp=122.932d0
     real*8 :: ragpr,decgpr,lcpr
     real*8 :: ramgp,lcpml
  
     ragpr=ragp*d2rd
     decgpr=decgp*d2rd
     lcpr=lcp*d2rd
     rar=ra*d2rd
     decr=dec*d2rd
     lr=l*d2rd
     br=b*d2rd
     ramgp=rar-ragpr
     lcpml=lcpr-lr
     
     sindec=sin(decgpr)*sin(br)+cos(decgpr)*cos(br)*cos(lcpml)
     decr=asin(sindec)
     dec=decr*r2dd
     cosdec=cos(decr)
     if(cosdec.eq.0.) then
        write(*,*) "Pole"
     else
        sinra=cos(br)*sin(lcpml)/cosdec
        cosra=(cos(decgpr)*sin(br)-sin(decgpr)*cos(br)*cos(lcpml))/cosdec
        if(sinra.gt.1d0) sinra=1d0
        if(sinra.lt.-1d0) sinra=-1d0
        if(sinra.ge.0d0) then
          if(cosra.ge.0d0) then
            ra=ragp+asin(sinra)*r2dd
          else
            ra=ragp+180d0-asin(sinra)*r2dd
          end if
        else
          if(cosra.ge.0d0) then
            ra=ragp+360d0+asin(sinra)*r2dd
          else 
            ra=ragp+180d0-asin(sinra)*r2dd
          end if
        end if
     end if
     if(ra.ge.360d0) ra=ra-360d0
     if(ra.lt.0d0) ra=ra+360d0
     return
  end subroutine galactic2equator
 
  subroutine Sphere2Euclide(r,ra,dec,x,y,z)
    implicit none
    real*8 :: ra,dec,r
    real*8 :: x,y,z
    x=r*cos(dec*d2rd)*cos(ra*d2rd)
    y=r*cos(dec*d2rd)*sin(ra*d2rd)
    z=r*sin(dec*d2rd)
    return
  end subroutine

  subroutine Euclide2Sphere(x,y,z,r,ra,dec)
    implicit none
    real*8 :: x,y,z,a,b
    real*8 :: r,ra,dec
    real*8 :: cosra,sinra,cosdec
    r=sqrt(x**2+y**2+z**2)
    if(r.eq.0) then
      write(*,*) "Warning: r is 0!"
      return
    end if
    a=z/r
    if(a.gt.1d0) a=1d0
    if(a.lt.-1d0) a=-1d0
    dec = asin(a)
    if(abs(a).eq.1d0) then
      ra=0
      dec=dec*r2dd
      return
    end if
    cosdec = cos(dec)
    cosra = x/r/cosdec
    sinra = y/r/cosdec
    if(cosra.gt.1d0) cosra=1d0
    if(cosra.lt.-1d0) cosra=-1d0
    if(sinra.gt.1d0) sinra=1d0
    if(sinra.lt.-1d0) sinra=-1d0
    ra = (1d0-sign(1d0,sinra))*pid+sign(1d0,sinra)*acos(cosra)
    ra = ra*r2dd
    dec = dec*r2dd
    return
  end subroutine

  subroutine angle2vector(ra,dec,vout)
    implicit none
    real*8 :: ra,dec,vout(3)
    vout(1)=cos(dec*d2rd)*cos(ra*d2rd)
    vout(2)=cos(dec*d2rd)*sin(ra*d2rd)
    vout(3)=sin(dec*d2rd)
    return
  end subroutine

  subroutine vector2angle(vin,ra,dec)
    implicit none
    real*8 :: vin(3),ra,dec,r,x,y,z,a,b
    real*8 :: cosra,sinra,cosdec
    r=sqrt(sum(vin**2))
    if(r.gt.1d0) then
      write(*,*) "Warning: amplitude of input vector is larger than 1."
      write(*,*) "It will be set 1.",r
      x=vin(1)/r
      y=vin(2)/r
      z=vin(3)/r
      r=1d0
    end if
    if(r.eq.0d0) then
      write(*,*) "Warning: r is 0!"
      return
    end if
    if(z.gt.1d0) z=1d0
    if(z.lt.-1d0) z=-1d0
    dec = asin(z)
    if(abs(z).eq.1d0) then
      ra=0
      dec=dec*r2dd
      return
    end if
    cosdec = cos(dec)
    cosra = x/r/cosdec
    sinra = y/r/cosdec
    if(cosra.gt.1d0) cosra=1d0
    if(cosra.lt.-1d0) cosra=-1d0
    if(sinra.gt.1d0) sinra=1d0
    if(sinra.lt.-1d0) sinra=-1d0
    ra = (1d0-sign(1d0,sinra))*pid+sign(1d0,sinra)*acos(cosra)
    ra = ra*r2dd
    dec = dec*r2dd
    return
  end subroutine

  subroutine rotation_zxz(vin,vout,alpha,beta,gama)
    implicit none
    real*8 :: vin(3),vout(3),alpha,beta,gama,vt(3),vtt(3)
    call rotation_z(vin,vt,alpha)
    call rotation_x(vt,vtt,beta)
    call rotation_z(vtt,vout,gama)
    return
  end subroutine

  subroutine rotation_x(vin,vout,theta)
    implicit none
    real*8 :: vin(3),vout(3),theta,r,costh,sinth
    r=theta*d2rd
    costh=cos(r)
    sinth=sin(r)
    vout(1)=vin(1)
    vout(2)=vin(2)*costh-vin(3)*sinth
    vout(3)=vin(2)*sinth+vin(3)*costh
    return
  end subroutine
  subroutine rotation_y(vin,vout,theta)
    implicit none
    real*8 :: vin(3),vout(3),theta,r,costh,sinth
    r=theta*d2rd
    costh=cos(r)
    sinth=sin(r)
    vout(2)=vin(2)
    vout(1)=vin(1)*costh+vin(3)*sinth
    vout(3)=vin(3)*costh-vin(1)*sinth
    return
  end subroutine
  subroutine rotation_z(vin,vout,theta)
    implicit none
    real*8 :: vin(3),vout(3),theta,r,costh,sinth
    r=theta*d2rd
    costh=cos(r)
    sinth=sin(r)
    vout(3)=vin(3)
    vout(1)=vin(1)*costh-vin(2)*sinth
    vout(2)=vin(1)*sinth+vin(2)*costh
    return
  end subroutine

end module coordtransform


module distances
  use const
  use integral
  implicit none
  type cosmoparams
     real*8 :: omega_m,omega_lambda,h
  end type
  type(cosmoparams) :: params4cdist
  contains
  !========================================================================= 
  ! -------- Distance calculators --------
  ! A. angle distance on a unit shpere: "angledist(theta,ra1,dec1,ra2,dec2)"
  ! B. Comoving distance from z=0 to z: "comovingdist(r,z,om,omL,h)"
  ! C. Line-of-sight and perpendicular seperation: "los_proj(s,rp,pii,ra1,dec1,r1,ra2,dec2,r2)"
  !     1. s -- output, = sqrt(r1^2+r2^2-r1*r2*cos(theta))
  !     2. rp -- output, perpendicular-to-los distance
  !     3. pii -- output, Los distance
  !=======================================================
   function angledistance(ra1,dec1,ra2,dec2) result(res)
     implicit none
     real*8 :: res,ra1,dec1,ra2,dec2,costheta
     costheta=cos(dec1*d2rd)*cos(dec2*d2rd)*cos((ra1-ra2)*d2rd)&
                +sin(dec1*d2rd)*sin(dec2*d2rd)
     res=acos(costheta)*r2dd
   end function

   ! Comoving distance in a flat lambda-CDM model
   ! Unit Mpc/h.
   function comovingdistance(om,omL,z) result(res)
     implicit none
     real*8 :: res,z,om,omL,epsr=1d-4,epsa=0d0,abserr,zmin=0d0
     integer :: n=2,key=6,nevals
     params4cdist%omega_m=om
     params4cdist%omega_lambda=omL
     call one_integral(n,func_cdist,zmin,z,res,epsr,epsa,key,nevals,abserr)
     res=res*LightSpeed/100d0
   end function
   function func_cdist(x) result(res)
     implicit none
     real*8,intent(in) :: x
     real*8 :: res,om,oml
     om=params4cdist%omega_m
     oml=params4cdist%omega_lambda
     res=1d0/(sqrt(oml+om*(1d0+x)**3))
   end function

   ! Line-of-sight and projected distance of two points given ra,dec,and r.
   ! All angles in degrees.
   subroutine los_proj(s,rp,pii,ra1,dec1,r1,ra2,dec2,r2)   
     implicit none
     real :: r1,ra1,dec1,z1,r2,ra2,dec2,z2
     real :: s,pii,rp
     real :: costheta,s2
    
     costheta=cos(dec1*d2r)*cos(dec2*d2r)*cos((ra1-ra2)*d2r)&
                +sin(dec1*d2r)*sin(dec2*d2r)  
     s2=r1**2+r2**2-2*r1*r2*costheta
     pii=abs(r1**2-r2**2)/sqrt(r1**2+r2**2+2*r1*r2*costheta)
     rp=sqrt(s2-pii**2)
     s=sqrt(s2)
     return
   end subroutine los_proj

   ! LOS projection. 
   subroutine LOS_projection(v1,v2,rp2,rlos2)
     implicit none
     real*8 :: v1(3),v2(3),vlos(3),vdiff(3),rp2,rlos2
     real*8 :: rplus2,rdiff2
     vlos=v2+v1
     vdiff=v2-v1
     rplus2=sum(vlos*vlos)
     if (rplus2.eq.0) then
        rlos2=sum(vdiff*vdiff)
        rp2=0d0
        return
     end if
     rlos2=(sum(vlos*vdiff))**2/rplus2
     rdiff2=sum(vdiff*vdiff)
     rp2=rdiff2-rlos2
   end subroutine

end module distances


module randomx
  implicit none
  contains
    !================================================================================
    ! Random number generator for "uniform","Normal" and "Poisson" nembers
    ! A. Set a large integer "iseed" (smaller than 2^30).
    ! B. Initialize : 'call randini(iseed)'.
    ! C. Uniform ==== "call randa(x)"
    !    Normal  ==== "Call randg(x)"
    !    Poisson ==== "call randp(alam,n)", alam is mean.
    !================================================================================
    ! Written by E. Bertschinger, 1992.
    !================================================================================
    subroutine randini(iseed)
      !c  randini is used to initialize the random number generators.  iseed is a
      !c  positive integer.  The basic random number generator is a
      !c  subtract-with-borrow lagged Fibonacci generator with base b=2**32-5 and
      !c  period b**r-b**s=10**414.
      !c  See Marsaglia and Zaman, Ann. Appl. Prob. 1, 462 (1991).
      !c  This generator is shuffled with an independent one to break up serial
      !c  correlations.
      integer iseed,r,s,c,nshuf,n,irptr,i,j,nmr,nms,inorm
      parameter (r=43,s=22,nshuf=157)
      real*8 randar(nshuf),x,xnorm
      double precision xt(r),xn,b,binv
      common /randnos/ xt,randar,irptr,n,c
      common /randnor/ xnorm,inorm
      save /randnos/,/randnor/
      parameter (b=4294967291.0d0,binv=1.0d0/b)

      xn=iseed
      xn=dmod(xn,b)
      if (xn.lt.0) xn=xn+b
      if (xn.lt.0.or.xn.ge.b) then
        write(*,*) 'Error with seed in randini'
        stop
      end if
      xt(1)=xn
      !c  Initialize xt using multiplicative congruential generator.
      do 10 i=2,r
        xn=dmod(16807.0d0*xn,b)
        xt(i)=xn
      10 continue
      n=r
      c=0
      !c  Warm up generator.
      do 30 j=1,5
        do 20 i=1,nshuf
        n=n+1
        if (n.gt.r) n=1
        nmr=n-r
        if (nmr.lt.1) nmr=nmr+r
        nms=n-s
        if (nms.lt.1) nms=nms+r
        xn=xt(nms)-xt(nmr)-c
        if (xn.ge.0) then
          c=0
        else
          c=1
          xn=xn+b
        end if
        xt(n)=xn
      !c  Fill shuffling table.
          randar(i)=xn*binv
        20  continue
      30  continue
      irptr=1
      !c  Initialize shuffling generator.
      call randin1(iseed)
      !c  Run a few times to shuffle.
      do 40 i=1,5*nshuf
        call randa(x)
      40  continue
      xnorm=0.0
      inorm=0
      return
    end subroutine
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine randa(x)
      integer r,s,c,nshuf,n,irptr,nmr,nms
      parameter (r=43,s=22,nshuf=157)
      real*8 randar(nshuf),x,x1
      double precision xt(r),xn,b,binv
      common /randnos/ xt,randar,irptr,n,c
      save /randnos/
      parameter (b=4294967291.0d0,binv=1.0d0/b)

      !c  Extract random number from shuffling table.
      10  x=randar(irptr)
      !c  Generate a new random number.
      n=n+1
      if (n.gt.r) n=1
      nmr=n-r
      if (nmr.lt.1) nmr=nmr+r
      nms=n-s
      if (nms.lt.1) nms=nms+r
      xn=xt(nms)-xt(nmr)-c
      if (xn.ge.0) then
        c=0
      else
        c=1
        xn=xn+b
      end if
      xt(n)=xn
      !c  Refill shuffling table.
      randar(irptr)=xn*binv
      !c  Use independent generator to get pointer into shuffling table.
      call randa1(x1)
      irptr=int(nshuf*x1)+1
      irptr=min(irptr,nshuf)
      return
    end subroutine
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine randg(x)
      !c  randg generates standard normal deviates.  It uses randa, which must be
      !c  initialized by a call to randini.
      real*8 x,xr,xnorm
      real*8,parameter :: twopi=6.28318531
      integer :: inorm
      !parameter (twopi=6.28318531)
      common /randnor/ xnorm,inorm
      save /randnor/

      if (inorm.eq.1) then
        x=xnorm
        inorm=0
      else
      10  call randa(x)
        if (x.le.0.0) go to 10
        xr=sqrt(-2.0*log(x))
        call randa(x)
        xnorm=xr*sin(twopi*x)
        x=xr*cos(twopi*x)
        inorm=1
      end if
      return
    end subroutine
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine randp(alam,np)
      !c  randp generates a sample np from a Poisson distribution with mean alam.
      !c  It uses randa, which must be initialized by a call to randini.
      real*8 alam,x,t,tt
      integer :: np
      if (alam.le.0.0) then
        np=-1
        return
      end if

      if (alam.gt.50.0) then
      !c  Use asymptotic Gaussian distribution for alam > 50.
        call randg(x)
        np=alam+sqrt(alam)*x
        return
      end if

      np=0
      t=1.0
      tt=exp(-alam)
      10  call randa(x)
        t=t*x
        if (t.lt.tt) go to 20
        np=np+1
        go to 10
      20 return
    end subroutine
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine randin1(iseed)
      !c  randini is used to initialize the random number generators.  iseed is a
      !c  positive integer less than 1.0e9.  The basic random number generator is
      !c  taken from Press et al., Numerical Recipes, p. 199 and is based on
      !c  Knuth s suggestion for a portable random number generator.  mseed is
      !c  any large number less than m=1.0e9.
      integer,parameter :: m=1000000000,mseed=161803398,nrand=157
      real,parameter :: rm=1.0/m
      !parameter (rm=1.0/m,)
      !c  The dimension 55 is special and should not be modified; see Knuth.
      integer ma(55),iseed,isee,mj,mk,i,ii,inext,inextp,irptr,k
      real*8 randar(nrand)
      common /randno1/ randar,irptr,ma,inext,inextp
      save /randno1/

      isee=mod(iseed,m)
      !c  Initialize ma(55).
      mj=mseed-isee
      mj=mod(mj,m)
      if (mj.lt.0) mj=mj+m
      ma(55)=mj
      mk=1
      !c  Now initialize the rest of the table, in a slightly random order,
      !c  with numbers that are not especially random.
        do 10 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if (mk.lt.0) mk=mk+m
        mj=ma(ii)
      10  continue
      !c  Randomize them by "warming up the generator."
        do 30 k=1,4
          do 20 i=1,55
          ma(i)=ma(i)-ma(1+mod(i+30,55))
          if (ma(i).lt.0) ma(i)=ma(i)+m
      20    continue
      30  continue
      inext=0
      inextp=31
      !c  Exercise generator before storing in shuffling table.
        do 40 i=1,nrand
        inext=inext+1
        if (inext.eq.56) inext=1
        inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if (mj.lt.0) mj=mj+m
        ma(inext)=mj
      40  continue
      !c  Now fill shuffling table.
        do 50 i=1,nrand
        inext=inext+1
        if (inext.eq.56) inext=1
        inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if (mj.lt.0) mj=mj+m
        ma(inext)=mj
        randar(i)=mj*rm
      50  continue
      irptr=1
      return
    end subroutine
    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine randa1(x)
      !c  randa generates uniform random numbers in the interval [0,1).
      !c  It must be initialized with a call to randin1.
      integer,parameter :: m=1000000000,mseed=161803398,nrand=157
      real,parameter :: rm=1.0/m
      !parameter (m=1000000000,mseed=161803398,rm=1.0/m,nrand=157)
      !c  The dimension 55 is special and should not be modified; see Knuth.
      integer ma(55),mj,inext,inextp,irptr
      real*8 randar(157),x
      common /randno1/ randar,irptr,ma,inext,inextp
      save /randno1/

      !c  Extract random number from shuffling table.
      x=randar(irptr)
      !c  Generate a new random number.
      inext=inext+1
      if (inext.eq.56) inext=1
      inextp=inextp+1
      if (inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if (mj.lt.0) mj=mj+m
      ma(inext)=mj
      !c  Shuffle.
      randar(irptr)=mj*rm
      irptr=int(nrand*x)+1
      irptr=min(irptr,nrand)
      return
    end subroutine

    !=========================================================================
    ! Random number generator for "uniform","Gaussian" and "Normal" nembers
    ! A. set a large integer "iseed".
    ! B. initialize "ran2": xx=ran2(iseed).
    ! C. x = ran2(iseed) or gauss_ran2(iseed,mean,sigma) or normal_ran2(iseed)
    !=========================================================================
    ! 'ran2(iseed)': Generate "uniform" random numbers within range (0,1)
    !    1. iseed -- input, seed
    !    2. need initialize, xx=ran2(iseed)
    ! 'gauss_ran2(iseed,mean,sigma)': "Gaussian" random numbers
    !    1. iseed --input, seed
    !    2. mean,sigma -- input
    ! 'normal_ran2(iseed)': "uniform" random numbers
    !    1. iseed -- input, seed
    !=========================================================================
    FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
         idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
    END function ran2

   function gauss_ran2(idum,mean,sigma)
     implicit none
     integer :: idum
     real*8 :: mean,sigma,gauss_ran2
     gauss_ran2 = normal_ran2(idum)*sigma + mean
   end function

   function normal_ran2(idum)
     implicit none
     integer :: idum
     real*8 :: normal_ran2
     real*8 :: v1,v2,rsq
     do while(.true.)
       v1 = ran2(idum)*2d0 - 1d0
       v2 = ran2(idum)*2d0 - 1d0
       rsq = v1**2 + v2**2
       if(0.lt.rsq.and.rsq.lt.1d0) exit
     end do
     normal_ran2 = v1*sqrt(-2d0*log(rsq)/rsq)
   end function 

end module


module hunt
  implicit none
  contains
  !=======================================================================
  ! 'HUNTD(R)(xx,x,jlo)': find jlo such that xx(jlo) < x <= xx(jlo+1), 
  !    1. xx(:) -- input increasing array. (one-indexed)
  !    2. x -- input.
  !    3. jlo  -- input/output. 
  !               input is taken as the initial guess.
  !               Output 0 or size(xx) if x is out of range. 
  !=======================================================================
  ! 'HUNTJ(ind,ip,ipout)': find whether a given integer ip is in ind.
  !    1. ind(:) -- input array with increasing order. (one-indexed)
  !    2. ip -- input
  !    3. ipout -- output, the sequence number of ip in index.
  !                        -1 if ip is not in ind. 
  ! Note :: "xx","ind" should be one-indexed.
  !=======================================================================
    subroutine HUNTD(xx,x,jlo)
      implicit none
      real(8), dimension(:), intent(in) :: xx
      real(8), intent(in) :: x
      integer, intent(inout) :: jlo
      integer :: n,inc,jhi,jm
      logical :: hunting

      n = size(xx)
      if (jlo<1) jlo= 1
      if (jlo>n) jlo= n
      inc= 1
      hunting= .true.
      if (xx(jlo)<=x) then
         do while (hunting)
            jhi= jlo+inc
            if (jhi>n) then
               jhi= n
               hunting= .false.
            else
               if (xx(jhi)<=x) then
                  jlo= jhi
                  inc= inc+inc
               else
                  hunting= .false.
               endif
            endif
         end do
      else
         jhi= jlo
         do while (hunting)
            jlo= jhi-inc
            if (jlo<1) then
               jlo= 0
               hunting= .false.
            else
               if (xx(jlo)>x) then
                  jhi= jlo
                  inc= inc+inc
               else
                  hunting= .false.
               endif
            endif
         end do
      endif
      do while ((jhi-jlo)>1)
         jm= (jhi+jlo)/2
         if (xx(jm)<x) then
            jlo= jm
         else
            jhi= jm
         endif
      end do
    end subroutine HUNTD
    subroutine HUNTR(xx,x,jlo)
      implicit none
      real, dimension(:), intent(in) :: xx
      real, intent(in) :: x
      integer, intent(inout) :: jlo
      integer :: n,inc,jhi,jm
      logical :: hunting
  
      n = size(xx)
      if (jlo<1) jlo= 1
      if (jlo>n) jlo= n
      inc= 1
      hunting= .true.
      if (xx(jlo)<=x) then
         do while (hunting)
            jhi= jlo+inc
            if (jhi>n) then
               jhi= n
               hunting= .false.
            else
               if (xx(jhi)<=x) then
                  jlo= jhi
                  inc= inc+inc
               else
                  hunting= .false.
               endif
            endif
         end do
      else
         jhi= jlo
         do while (hunting)
            jlo= jhi-inc
            if (jlo<1) then
               jlo= 0
               hunting= .false.
            else
               if (xx(jlo)>x) then
                  jhi= jlo
                  inc= inc+inc
               else
                  hunting= .false.
               endif
            endif
         end do
      endif
      do while ((jhi-jlo)>1)
         jm= (jhi+jlo)/2
         if (xx(jm)<x) then
            jlo= jm
         else
            jhi= jm
         endif
      end do
    end subroutine HUNTR

    subroutine HUNTJ(ind,ip,ipout)
      implicit none
      integer,dimension(:),intent(in) :: ind
      integer,intent(in) :: ip
      integer,intent(out) :: ipout
      integer :: i,j,k,indi,indj,indk,n

      n=size(ind)
      i=1
      j=n
      ipout = -1 ! initialize ipout.
      indi = ind(i)
      indj = ind(j)
      do while(.true.)
 999     if (j.eq.(i+1)) goto 888
         if (indi.lt.ip.and.ip.lt.indj) then
            k=(i+j)/2
            indk=ind(k)
            if (ip.eq.indk) then
               ipout = k
               return
            else if(ip.lt.indk) then
               j=k
               indj=indk
            else
               i=k
               indi=indk
            end if
            goto  999
         else
888         if (indi.eq.ip) then
               ipout = i
            else if (indj.eq.ip) then
               ipout = j
            else
               ipout = -1
            end if
            return
         end if
      end do
    end subroutine HUNTJ
end module hunt


module interpolation
  implicit none
  contains
  !================================================================
  ! -------- Cubic spline interpolation ---------
  ! A. Call "spline" to calculate the second derivative y2. (Only once)
  ! B. Call "splint" to get interpolation value at given point.
  ! 'spline(x,y,n,yp1,ypn,yw)': 
  !    1. x(1:n),y(1:n) -- input, give function: yi=f(xi), x1<x2<...<xn
  !    2. yp1,ypn -- input, the first derivative of the interpolation at points 1 and n.
  !    3. If yp1 and/or ypn >= 1e30, the routine will set the corresponding boundary for 
  !       a natural spline, with zero second derivative on point 1 and/or n.
  ! 'splint(xa,ya,y2a,n,x,y)':
  !    1. xa(1:n),ya(1:n) -- input, give function y=f(x), x1<x2<...<xn
  !    2. ya2(1:n) -- second derivatives calculated by "spline".
  !    3. x -- input, value.  
  !    4. y -- output, returned f(x).
  !================================================================
  subroutine spline(x,y,n,yp1,ypn,y2)
    implicit none
    integer :: n
    real :: yp1,ypn,x(n),y(n),y2(n)
    integer,parameter :: nmax = 800 ! The largest anticipated value of n.
    integer :: i,k
    real :: p,qn,sig,un,u(nmax)
    
    if(yp1.gt..99e30) then
      y2(1) = 0.   ! The lower boundary condition is set 
	u(1) = 0.    ! either to be 'natural'.
    else
      y2(1) = -0.5 ! Or else to have a specified first derivative.
	u(1) = 3./(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1)) -yp1)
    end if
    do i=2,n-1
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
	p = sig*y2(i-1)+2.
      y2(i) = (sig-1.)/p
      u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
	        /(x(i+1)-x(i-1)) - sig*u(i-1))/p
    end do
    if(ypn.gt..99e30) then
      qn=0.       ! The upper boundary condition is set
	un=0.       ! either to be 'natural'.
    else          
      qn=0.5      ! Or else to have a specified first derivative.
	un= 3./(x(n)-x(n-1))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1  ! This the backsubstitution loopof the tridiagonal algorithm.
      y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    return
  end subroutine spline

  subroutine splint(xa,ya,y2a,n,x,y)
    implicit none
    integer :: n
    real :: x,y,xa(n),ya(n),y2a(n)
    integer :: k,khi,klo
    real :: a,b,h

    klo = 1
    khi = n
    101 if(khi-klo.gt.1) then
       k=(khi+klo)/2
	 if(xa(k).gt.x) then
	   khi=k
	 else
	   klo=k
	 end if
	 goto 101
    end if              ! klo and khi now bracket the input value x
    h = xa(khi)-xa(klo) 
    if(h.eq.0) pause 'bad xa input in splint' ! xa must be distinct.
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h**2/6.
    return
  end subroutine splint

  !----- for double precision -----
  subroutine splined(x,y,n,yp1,ypn,y2)
    implicit none
    integer :: n
    real*8 :: yp1,ypn,x(n),y(n),y2(n)
    integer,parameter :: nmax = 800 ! The largest anticipated value of n.
    integer :: i,k
    real*8 :: p,qn,sig,un,u(nmax)
    
    if(yp1.gt..99e30) then
      y2(1) = 0d0   ! The lower boundary condition is set 
	u(1) = 0d0    ! either to be 'natural'.
    else
      y2(1) = -0.5d0 ! Or else to have a specified first derivative.
	u(1) = 3d0/(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1)) -yp1)
    end if
    do i=2,n-1
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
	p = sig*y2(i-1)+2d0
      y2(i) = (sig-1d0)/p
      u(i) = (6d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
	        /(x(i+1)-x(i-1)) - sig*u(i-1))/p
    end do
    if(ypn.gt..99d30) then
      qn=0d0       ! The upper boundary condition is set
	un=0d0       ! either to be 'natural'.
    else          
      qn=0.5d0      ! Or else to have a specified first derivative.
	un= 3d0/(x(n)-x(n-1))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1d0)
    do k=n-1,1,-1  ! This the backsubstitution loopof the tridiagonal algorithm.
      y2(k) = y2(k)*y2(k+1)+u(k)
    end do
    return
  end subroutine splined
  subroutine splintd(xa,ya,y2a,n,x,y)
    implicit none
    integer :: n
    real*8 :: x,y,xa(n),ya(n),y2a(n)
    integer :: k,khi,klo
    real*8 :: a,b,h

    klo = 1
    khi = n
    101 if(khi-klo.gt.1) then
       k=(khi+klo)/2
	 if(xa(k).gt.x) then
	   khi=k
	 else
	   klo=k
	 end if
	 goto 101
    end if              ! klo and khi now bracket the input value x
    h = xa(khi)-xa(klo) 
    if(h.eq.0) pause 'bad xa input in splint' ! xa must be distinct.
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*h**2/6d0
    return
  end subroutine splintd


end module interpolation



module special_function
  use integral
  implicit none
  real*8 :: param_IncompleteGammaFunction
  contains
  !=======================================================
  ! Legendre Polynomial P(L,x),use iteration formula:
  !   P(L,x) = (2*L+1)*x*P(L-1,x) - L/(L+1)*P(L-2,x)
  !   and  P(0,x) = 1, P(1,x) = x
  ! Return are values of P(L,x) between Lmin and Lmax.
  !=======================================================
  subroutine LegendrePolynomial(x,Lmin,Lmax,res)
    implicit none
    real,intent(in) :: x
    integer,intent(in) :: Lmin,Lmax
    real,intent(out) :: res(*)
    real :: f(Lmax+1)
    integer :: L,nL

    nL = Lmax - Lmin + 1
    f(0) = 1.
    f(1) = x
    do L=2,Lmax
      f(L) = (2*L+1)*x*f(L-1)/(L+1) - L*f(L-2)/(L+1)
    end do

    res(1:nL) = f(Lmin:Lmax)
    return
  end subroutine LegendrePolynomial
  ! ### Legendre Polynomial for double precision.
  subroutine LegendrePolynomiald(x,Lmin,Lmax,res)
    implicit none
    real*8,intent(in) :: x
    integer,intent(in) :: Lmin,Lmax
    real*8,intent(out) :: res(*)
    real*8 :: f(Lmax+1)
    integer :: L,nL

    nL = Lmax - Lmin + 1
    f(0) = 1d0
    f(1) = x
    do L=2,Lmax
      f(L) = (2*L+1)*x*f(L-1)/(L+1) - L*f(L-2)/(L+1)
    end do

    res(1:nL) = f(Lmin:Lmax)
    return
  end subroutine LegendrePolynomiald

  FUNCTION gammln(xx)
    REAL gammln,xx
    INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof/76.18009172947146d0,-86.50532032941677d0,&
    &24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
    &-.5395239384953d-5/

    stp = 2.5066282746310005d0
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do  j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
    end do
    gammln=tmp+log(stp*ser/x)
    return
  END function

  ! Approximation of Integral of [t^(a-1) * e^(-t)] over (b1,b2);
  ! Here we cut down the integral at t = large = (a - 1)*4 - log(small) if a>1 
  !                                              (a - 1)*2 - log(small) if a<1
  ! It will be Incomplete Gamma Function with parameters a and b1 when b2->infinity.
  subroutine IncompleteGammaFunction(a,b1,b2,res,small)
    implicit none
    real*8,intent(in) :: a,b1,b2  ! a is the power of t, b is the lower integral limit.
    real*8,intent(out) :: res
    real*8,intent(in),optional :: small  ! The 
    real*8 :: large = 30d0  !  e^(-30) = 1e-13, small enough; e^(-40) = 4.24*10^-18
    real*8 :: epsr,epsa,abserr,lower,upper,ss
    integer :: n,neval,key

    if(present(small)) then
      ss = log(small)
    else
      ss = -1d0*large
    end if
    if(a.gt.1) then
      large = (a - 1.)*4.0 - ss
    else
      large = (a - 1.)*2.0 - ss
    end if
     
    epsr = 1d-4  ! relative error.
    epsa = 1d-4  ! absolute error.
    n = 2 ! Adaptive Gauss-Kronrod integral
    key = 4 ! (20-41) point Gauss-Kronrod rule use in adaptive dgage.

    param_IncompleteGammaFunction = a  ! Transform a to Integrand_of_GammaFunction.

    if(b1.ge.large) then
      res = 0
      write(*,*) "b1=",b1,"large=",large
    else if(b2.ge.large) then
      lower = b1
      upper = large
      call one_integral(n,Integrand_of_GammaFunction,lower,upper,res,epsr,epsa,key,neval,abserr)
    else 
      lower = b1
      upper = b2
      call one_integral(n,Integrand_of_GammaFunction,lower,upper,res,epsr,epsa,key,neval,abserr)
    end if
  end subroutine 

  function Integrand_of_GammaFunction(x) result(res)
    implicit none
    real*8,intent(in) :: x
    real*8 :: a,res
    a = param_IncompleteGammaFunction ! Accept a.
    res = x**(a-1d0)*exp(-x)
  end function 

  !-----  sine and cosine integral -----
  !-----  In NR
  !-----  Relative error (EPSR) and maximum iterations (Nits) are allowed to input.
  !-----  If absent, eps=6d-8 and maxit=100
  subroutine cisi(x,ci,si,epsr,nits)
    implicit none
    real*8 :: ci,si,x
    real*8,intent(in),optional :: epsr
    integer,intent(in),optional :: nits
    real*8,parameter :: euler = 0.57721556d0 ! Euler's constant
    real*8,parameter :: piby2 = 1.5707963    ! pi/2
    real*8,parameter :: fpmin = 1d-30   ! near smallest representable floating-point number
    real*8,parameter :: tmin = 2d0      ! the dividing line between using the series and continued fraction
    integer :: i,k
    real*8 :: a,err,fact,signa,suma,sumc,sums,t,term,absc,eps,maxit
    double complex :: h,b,c,d,del
    logical :: odd

    absc(h) = dabs(dble(h))+dabs(aimag(h))  ! statement function

    if(present(nits))  then
      maxit = nits
    else
      maxit = 100
    end if
    if(present(epsr)) then
      eps = epsr
    else
      eps = 6d-8
    end if

    t = dabs(x)
    if(t.eq.0) then
      si = 0
      ci = -1d0/fpmin
      return
    end if
    if(t.gt.tmin) then
      b = cmplx(1d0,t)
      c = 1d0/fpmin
      d = 1d0/b
      h = d
      do i=2,maxit
        a = -(i-1)**2
        b = b+2
        d = 1d0/(a*d+b)
        c = b+a/c
        del = c*d
        h = h*del
        if(absc(del-1d0).lt.eps) goto 101
      end do
      pause 'cf failed in cisi ...'
      101 continue
      h = cmplx(dcos(t),-dsin(t))*h
      ci = -dble(h)
      si = piby2 + aimag(h)
    else
      if(t.lt.sqrt(fpmin))then
        sumc = 0d0
        sums = t
      else
        suma = 0
        sums = 0
        sumc = 0
        signa = 1
        fact = 1
        odd = .true.
        do k=1,maxit
          fact = fact*t/k
          term = fact/k
          suma = suma+signa*term
          err = term/dabs(suma)
          if(odd) then
            signa = -signa
            sums = suma
            suma = sumc
          else
            sumc = suma
            suma = sums
          end if
          if(err.lt.eps) goto 102
          odd = .not.odd
        end do
        pause 'maxits exceeded in cisi .....'
      end if
      102 si = sums
      ci = sumc+dlog(t)+euler
    end if
    if(x.lt.0) si = -si
    return
  end subroutine cisi
end module special_function


module matrix
  implicit none
  contains

  subroutine inverse(a,c,n)
    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input:
    !   a(n,n) -- array of coefficients for matrix A
    !   n      -- dimension
    ! output:
    !   c(n,n) -- inverse matrix of A
    ! comments:
    !   The original matrix a(n,n) will be destroyed 
    !   during the calculation
    !===========================================================
    implicit none 
    integer n
    double precision a(n,n), c(n,n)
    double precision L(n,n), U(n,n), b(n), d(n), x(n)
    double precision coeff
    integer i, j, k
    
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 allows such operations on matrices
    L=0.0
    U=0.0
    b=0.0
    
    ! step 1: forward elimination
    do k=1, n-1
       do i=k+1,n
          coeff=a(i,k)/a(k,k)
          L(i,k) = coeff
          do j=k+1,n
             a(i,j) = a(i,j)-coeff*a(k,j)
          end do
       end do
    end do
    
    ! Step 2: prepare L and U matrices 
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
      do i=1,j
        U(i,j) = a(i,j)
      end do
    end do
  
    ! Step 3: compute columns of the inverse matrix C
    do k=1,n
      b(k)=1.0
      d(1) = b(1)
      ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,n
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
      ! Step 3b: Solve Ux=d using the back substitution
      x(n)=d(n)/U(n,n)
      do i = n-1,1,-1
        x(i) = d(i)
        do j=n,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
      ! Step 3c: fill the solutions x(n) into column k of C
      do i=1,n
        c(i,k) = x(i)
      end do
      b(k)=0.0
    end do
  end subroutine inverse
 
end module matrix


module mapprojection
  use const
  implicit none
  contains

  subroutine Lambert_equalarea_azimuthal_Nouthpole(theta,phi,x,y)
    implicit none
    real*8 :: theta,phi,x,y
    x=2d0*sin(theta*d2rd/2d0)*cos(phi*d2rd)
    y=2d0*sin(theta*d2rd/2d0)*sin(phi*d2rd)
    return
  end subroutine

  subroutine Lambert_equalarea_azimuthal_Nouthpole_inverse(x,y,theta,phi,iout)
    implicit none
    real*8 :: theta,phi,x,y,r,sinth,sinphi,cosphi,ad
    integer :: iout
    iout=0
    r=sqrt(x**2+y**2)
    sinth=r/2d0
    if(r.gt.2d0) then
      write(*,*) "Warning: input (x,y) is not the projection of any point on the sphere."
      iout=1
    else if(r.eq.0d0) then
      theta=0d0
      phi=0d0
      write(*,*) "Nouth Pole: theta=0 and phi=0."
    else if(r.eq.2d0) then
      theta=180d0
      phi=0d0
      write(*,*) "South Pole: theta=180 and phi=0."
    else
      theta=asin(r/2d0)*2d0*r2dd
      sinphi=y/r
      cosphi=x/r
      phi=(1d0-sign(1d0,sinphi))*pid+sign(1d0,sinphi)*acos(cosphi)
      phi=phi*r2dd
    end if
    return
  end subroutine

end module mapprojection
