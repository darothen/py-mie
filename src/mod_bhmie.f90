!This is Craig Bohren's bhmie code, modified by Dong-Chul Kim to
!calculate the assymetry parameter GSCA. Some additional minor 
!modifications by lz4ax for interfacing with bhmie_driver.

module bhmie_module

  use kinds_module

  private
  public :: bhmie_driver

  !implicit NONE
  !using implicit types as I'm too lazy to declare all variables
  !in Bohren's routine


  !number of angles between 0 and 90 degrees; matrix elements are 
  !calculated at 2*nang-1 angles, including 0, 90, and 180 degrees
  integer, parameter :: nang = 100
  real,    parameter :: ref_air = 1.0 ! air

CONTAINS

!=======================================================================
 subroutine bhmie_driver( radius, refr_ind, wavelength, &
                          Qext, Qsca, Qgsa )
!=======================================================================
!
! Purpose: bhmie driver. Calls bhmie routine to obtain the extinction
!          efficiency Qext, scattering efficiency Qsca and assymetry
!          parameter Qgsa for a particle with radius "radius" composed
!          of material with refractive index refre + j*refim
!
!-----------------------------------------------------------------------

   !use bhmie_mod, only: bhmie   
   implicit NONE

!----------------------------- arguments -------------------------------

   real (kind=rkind), intent(in) :: &
                      radius,   & ! particle radius [um]       
                      wavelength  ! in micron
   complex (kind=ckind), intent(in) :: &
                      refr_ind  ! particle complex refractive index

   real, intent(out) :: &
                      Qext,   & ! extinction efficiency
                      Qsca,   & ! scattering efficiency
                      Qgsa      ! assymetry parameter

!------------------------- local workspace -----------------------------

   !refractive index of the medium

   !integer :: iw, ir
   real    :: x, dum
   complex :: ref_ndx   

!-----------------------------------------------------------------------

  ! refractive index at specific wavelength
  ref_ndx = refr_ind / ref_air

  ! size parameter less than size
  dum = 2. * pi * ref_air / wavelength

  ! size parameter 
  x = dum * radius

  !get Mie coefficients for a particle with size radius(ir) and 
  call bhmie( x, ref_ndx, Qext, Qsca, Qgsa ) 

end subroutine bhmie_driver

! ===
!  Subroutine BHMIE calculates amplitude scattering matrix
!	elements and efficiencies for extinction, total scattering
!	and backscattering for a given size parameter and
!	relative refractive index
! ===

subroutine bhmie (x, refrel, qext, qsca, GSCA)
!subroutine bhmie (x,refrel,nang,s1,s2,qext,qsca,qback,GSCA)

  real, dimension(100) :: amu, theta, pi, tau, pi0, pi1
  complex :: d(3000), y, refrel, xi, xi0, xi1, an, bn, s1(200), s2(200), p, t
  double precision :: psi0, psi1, psi, dn, dx

!++dckim
  complex :: AN1,BN1
  GSCA=0.E0
!--dckim

!lz4ax: keep the gfortran compiler happy
  AN  = 0.0
  BN  = 0.0
  AN1 = 0.0
  BN1 = 0.0
!lz4ax:end
  
  dx = x
  y = x*refrel
  
  ! series terminated after nstop terms
  
  xstop = x + 4.0*x**0.3333 + 2.0
  nstop = xstop
  ymod  = cabs(y)
  nmx   = amax1(xstop, ymod) + 15
  dang  = 1.570796327/float(nang - 1)
  do 555 j=1,nang
  	theta(j) = (float(j) - 1.0)*dang
  555   amu(j)   = cos(theta(j))
  
  ! logarithmic derivative d(j) calculated by downward
  ! recurrence beginning with initial value 0.0 + i*0.0
  ! at j = nmx
  
  d(nmx) = cmplx (0.0, 0.0)
  nn = nmx - 1
  do 120 n=1,nn
  	rn = nmx - n + 1
  120	d(nmx-n) = (rn/y) - (1.0/(d(nmx-n+1) + rn/y))
  do 666 j=1,nang
  	pi0(j) = 0.0
  666	pi1(j) = 1.0
  nn = 2*nang - 1
  do 777 j=1,nn
  	s1(j) = cmplx (0.0, 0.0)
  777	s2(j) = cmplx (0.0, 0.0)

  ! Riccati-Bessel functions with real argument x
  ! calculated by upward recurrence
  
  psi0 = dcos(dx)
  psi1 = dsin(dx)
  chi0 = -sin(x)
  chi1 = cos(x)
  apsi0 = psi0
  apsi1 = psi1
  xi0  = cmplx (apsi0, -chi0)
  xi1  = cmplx (apsi1, -chi1)
  qsca = 0.0
  n    = 1
  
  200 dn = n
  rn = n
  fn = (2.0*dn + 1.0)/(rn*(rn + 1.0))
  psi = (2.0*dn - 1.0)*psi1/dx - psi0
  apsi= psi
  chi = (2.0*rn - 1.0)*chi1/x - chi0
  xi  = cmplx ( apsi, -chi)
!
!++dckim
!*** Store previous values of AN and BN for use
!    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
              AN1=AN
              BN1=BN
          ENDIF
!--dckim


  an  = (d(n)/refrel + rn/x)*apsi - apsi1
  an  = an/((d(n)/refrel + rn/x)*xi - xi1)
  bn  = (refrel*d(n) + rn/x)*apsi - apsi1
  bn  = bn/((refrel*d(n) + rn/x)*xi - xi1)
  qsca = qsca + (2.0*rn + 1.0)*(cabs(an)*cabs(an) + cabs(bn)*cabs(bn))

!++dckim
  !QSCA=QSCA+(2.*EN+1.)*(ABS(AN)**2+ABS(BN)**2)
!  GSCA=GSCA+((2.*rn+1.)/(rn*(rn+1.)))*                &
!              (REAL(AN)*REAL(BN)+AIMAG(AN)*AIMAG(BN))
!          IF(N.GT.1)THEN
!              GSCA=GSCA+((rn-1.)*(rn+1.)/rn)*         &
!             (REAL(AN1)*REAL(AN)+AIMAG(AN1)*AIMAG(AN)+  &
!              REAL(BN1)*REAL(BN)+AIMAG(BN1)*AIMAG(BN))
!          ENDIF
  GSCA=GSCA+((2.*rn+1.)/(rn*(rn+1.)))*                &
              (REAL(AN)*REAL(BN)+IMAG(AN)*IMAG(BN))
          IF(N.GT.1)THEN
              GSCA=GSCA+((rn-1.)*(rn+1.)/rn)*         &
             (REAL(AN1)*REAL(AN)+IMAG(AN1)*IMAG(AN)+  &
              REAL(BN1)*REAL(BN)+IMAG(BN1)*IMAG(BN))
          ENDIF
!--dckim

  do 789 j=1,nang
  	jj = 2*nang - j
	pi(j) = pi1(j)
	tau(j) = rn*amu(j)*pi(j) - (rn + 1.0)*pi0(j)
	p = (-1.0)**(n-1)
	s1(j) = s1(j) + fn*(an*pi(j) + bn*tau(j))
	t = (-1.0)**n
	s2(j) = s2(j) + fn*(an*tau(j) + bn*pi(j))
	if ( j.eq.jj) go to 789
	s1(jj) = s1(jj) + fn*(an*pi(j)*p + bn*tau(j)*t)
	s2(jj) = s2(jj) + fn*(an*tau(j)*t + bn*pi(j)*p)
  789 continue
  psi0 = psi1
  psi1 = psi
  apsi1 = psi1
  chi0 = chi1
  chi1 = chi
  xi1 = cmplx (apsi1, -chi1)
  n = n + 1
  rn = n
  do 999 j=1,nang
  	pi1(j) = ((2.0*rn - 1.0)/(rn - 1.0))*amu(j)*pi(j)
	pi1(j) = pi1(j) - rn*pi0(j)/(rn - 1.0)
  999	pi0(j) = pi(j)
  if (n-1-nstop) 200,300,300
!+-dckim, uncomment in the original code
!  300 qsca = (2.0/(x*x))*qsca

!++dckim
  300 GSCA=2.*GSCA/QSCA
  qsca = (2.0/(x*x))*qsca
!--dckim
  qext = (4.0/(x*x))*real(s1(1))
  qback = (4.0/(x*x))*cabs(s1(2*nang - 1))*cabs(s1(2*nang - 1))

end subroutine bhmie

end module	
