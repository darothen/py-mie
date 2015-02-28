module dmiess_module
!----------------------------------------------------------------------------------------------------
!
! This module is modified from a program which reads a lookup-table of
! core-shell Mie theory output for aerosol/droplet systems. The lookup-
! table analysis has been removed, leaving the core-shell theory calculation
! as implemented for
!
!       Toon and Ackerman, Applied Optics, Vol. 20, Pg. 3657
!
!----------------------------------------------------------------------------------------------------

  use kinds_module

  implicit none

  private 
  public :: dmiess_driver

contains

!=======================================================================
  subroutine dmiess_driver(radius_shell_um,radius_core_um,      &
                           refi_shell,refi_core, rad_lambda_um, &
                           QEXTc, QSCATc, asym                     )
!=======================================================================
!
! Purpose: dmiess driver. Calls the dmiess helper routines in order to
!          calculate the extinction efficency Qext, scattering
!          efficiency Qsca, and asymmetry parameter for a particle
!          consisting of a core and shell with total radius 
!          `radius_shell_um` and different composition / refractive
!          indices for the core and shell
!
!-----------------------------------------------------------------------
! NANG: Number of scattering angles between 0 and 90 degrees,
!        inclusive
! NOTE: NAN=2*NANG-1 is the number of scattering angles between
!       0 and 180 degrees, inclusive.
!----------------------------------------------------------


!----------------------------- arguments -------------------------------

    real (kind=rkind), intent(in) :: &
        radius_shell_um,  &  ! total particule radius, in microns
        radius_core_um,   &  ! particle core radius, in microns
        rad_lambda_um        ! wavelength of incident radiation, in microns

    complex (kind=ckind), intent(in) :: &
        refi_core,        &  ! complex refractive index of core
        refi_shell           ! complex refractive index of shell

    real (kind=rkind), intent(out) :: &
        qscatc,           &  ! scattering efficiency
        qextc,            &  ! extinction efficiency
        asym                 ! asymmetry parameter

!-- F2PY VARIABLE BINDINGS
    ! f2py intent(in) :: radius_shell_um, radius_core_um, rad_lambda_um
    ! f2py intent(in) :: refi_core, refi_shell
    ! f2py intent(out) :: qscatc, qextc, asym

!----------------------------------------------------------
    real (kind=rkind) :: VLAMBc, RGcmin, RGcmax, RGc, SIGMAGc, SHELRc, SHELIc, &
                         RINc, CORERc, COREIc,                                 &
                         QBACKc, EXTc, SCATc, BACKc, GSCAc,                    &
                         S11N, absorp, absorpbc, ext, scat, ssa, x, pi
    integer (kind=ikind) :: NRGFLAGc, NANG
    real (kind=rkind), dimension(200) :: ANGLESc, S1R, S1C, S2R, S2C, &
                                         S11, S12, S33, S34, SPOL, SP
    real mbc
!----------------------------------------------------------

    nang=2 ! only one angle
    nrgflagc=0 ! size distribution
    rgc=radius_shell_um     ! radius in microns

    vlambc=rad_lambda_um  ! radiation wavelength in microns
    rgcmin=0.001
    rgcmax=5.0
    sigmagc=1.0

    shelrc=real(refi_shell)
    shelic=aimag(refi_shell)

    corerc=real(refi_core)
    coreic=aimag(refi_core)

    RINc=radius_core_um/radius_shell_um

    pi = 4d0 * atan(1.0d0)
    x=2*pi*rgc/vlambc

    ! if ( debug ) then
    !   print *, "Shell radius", radius_shell_um
    !   print *, "Core radius", radius_core_um
    !   print *, "core fraction", RINc
    !   print *, "Shell ref", refi_shell
    !   print *, "Core ref", refi_core
    !   print *, "wavelength", rad_lambda_um

    !   print *, "VLAMBc", VLAMBc
    !   print *, "NRGFLAGc", NRGFLAGc
    !   print *, "RGcmin", RGcmin
    !   print *, "RGcmax", Rgcmax
    !   print *, "RGc", RGc
    !   print *, "SIGMAGc", SIGMAGc
    !   print *, "SHELRc", SHELRc
    !   print *, "SHELIc", SHELIc
    !   print *, "RINc", RINc
    !   print *, "CORERc", CORERc
    !   print *, "COREIc", COREIc
    !   print *, "NANG", NANG
    !   print *, "QEXTc", QEXTc
    !   print *, "QSCATc", QSCATc
    ! end if
    CALL ACKMIEPARTICLE( VLAMBc,NRGFLAGc,RGcmin,RGcmax,                 &
                         RGc,SIGMAGc,SHELRc,                            &
                         SHELIc, RINc,CORERc,COREIc,NANG,QEXTc,QSCATc,  &
                         QBACKc, EXTc,SCATc,BACKc, GSCAc,               &
                         ANGLESc,S1R,S1C,S2R,S2C,S11N,S11,S12,          &
                         S33,S34,SPOL,SP)

    ext=qextc*pi*(rgc*1.0e-04)**2  ! extinction cross section, cm^2
    scat=qscatc*pi*(rgc*1.0e-04)**2 ! scattering cross section, cm^2
    absorp=ext-scat ! absorption cross section, cm^2

    if(abs(absorp)/ext.le.0.00001)then  ! correct for round-off errors
      absorp=0.0
      scat=ext
    endif
    
    absorpbc=absorp
    mbc=1.333333*pi*(radius_core_um*1.e-4)**3*1.8  ! mass of BC, g,  assume density of 1.8 g/cm^3
    ssa=scat/ext
    asym=gscac

    ! if ( debug ) then
    !   print *, "-------------------------------------"
    !   print *, "RESULTS:"
    !   print *, "   sca - ", qscatc
    !   print *, "   ext - ", qextc
    !   print *, "  asym - ", asym
    ! end if

  end subroutine dmiess_driver

!**********************************************
!     /*--------------------------------------------------------*/
!     /* The Toon-Ackerman SUBROUTINE DMIESS for calculating the*/
!     /* scatter off of a coated sphere of some sort.           */
!     /* Toon and Ackerman, Applied Optics, Vol. 20, Pg. 3657   */
!     /*--------------------------------------------------------*/


      SUBROUTINE DMIESS(  RO,      RFR,     RFI,     THETD,     JX,   &
                          QEXT,    QSCAT,   CTBRQS,  ELTRMX,    PIE,   &
                          TAU,     CSTHT,   SI2THT,  ACAP, QBS, IT,   &
                          LL,      R,       RE2,     TMAG2,     WVNO  )
!
! **********************************************************************
!    THIS SUBROUTINE COMPUTES MIE SCATTERING BY A STRATIFIED SPHERE,
!    I.E. A PARTICLE CONSISTING OF A SPHERICAL CORE SURROUNDED BY A
!    SPHERICAL SHELL.  THE BASIC CODE USED WAS THAT DESCRIBED IN THE
!    REPORT: SUBROUTINES FOR COMPUTING THE PARAMETERS OF THE
!    ELECTROMAGNETIC RADIATION SCATTERED BY A SPHERE J.V. DAVE,
!    I B M SCIENTIFIC CENTER, PALO ALTO , CALIFORNIA.
!    REPORT NO. 320 - 3236 .. MAY 1968 .
!
!    THE MODIFICATIONS FOR STRATIFIED SPHERES ARE DESCRIBED IN
!        TOON AND ACKERMAN, APPL. OPTICS, IN PRESS, 1981
!
!    THE PARAMETERS IN THE CALLING STATEMENT ARE DEFINED AS FOLLOWS :
!      RO IS THE OUTER (SHELL) RADIUS;
!      R  IS THE CORE RADIUS;
!      RFR, RFI  ARE THE REAL AND IMAGINARY PARTS OF THE SHELL INDEX
!          OF REFRACTION IN THE FORM (RFR - I* RFI);
!      RE2, TMAG2  ARE THE INDEX PARTS FOR THE CORE;
!          ( WE ASSUME SPACE HAS UNIT INDEX. )
!      THETD(J): ANGLE IN DEGREES BETWEEN THE DIRECTIONS OF THE INCIDENT
!          AND THE SCATTERED RADIATION.  THETD(J) IS< OR= 90.0
!          IF THETD(J) SHOULD HAPPEN TO BE GREATER THAN 90.0, ENTER WITH
!          SUPPLEMENTARY VALUE, SEE COMMENTS BELOW ON ELTRMX;
!      JX: TOTAL NUMBER OF THETD FOR WHICH THE COMPUTATIONS ARE
!          REQUIRED.  JX SHOULD NOT EXCEED IT UNLESS THE DIMENSIONS
!          STATEMENTS ARE APPROPRIATEDLY MODIFIED;
!
!      THE DEFINITIONS FOR THE FOLLOWING SYMBOLS CAN BE FOUND IN LIGHT
!          SCATTERING BY SMALL PARTICLES, H.C.VAN DE HULST, JOHN WILEY
!          SONS, INC., NEW YORK, 1957.
!      QEXT: EFFIECIENCY FACTOR FOR EXTINCTION,VAN DE HULST,P.14 127.
!      QSCAT: EFFIECINCY FACTOR FOR SCATTERING,V.D. HULST,P.14 127.
!      CTBRQS: AVERAGE(COSINE THETA) * QSCAT,VAN DE HULST,P.128
!      ELTRMX(I,J,K): ELEMENTS OF THE TRANSFORMATION MATRIX F,V.D.HULST
!          ,P.34,45 125. I=1: ELEMENT M SUB 2..I=2: ELEMENT M SUB 1..
!          I = 3: ELEMENT S SUB 21.. I = 4: ELEMENT D SUB 21..
!      ELTRMX(I,J,1) REPRESENTS THE ITH ELEMENT OF THE MATRIX FOR
!          THE ANGLE THETD(J).. ELTRMX(I,J,2) REPRESENTS THE ITH ELEMENT
!          OF THE MATRIX FOR THE ANGLE 180.0 - THETD(J) ..
!      QBS IS THE BACK SCATTER CROSS SECTION.
!
!      IT: IS THE DIMENSION OF THETD, ELTRMX, CSTHT, PIE, TAU, SI2THT,
!          IT MUST CORRESPOND EXACTLY TO THE SECOND DIMENSION OF ELTRMX.
!      LL: IS THE DIMENSION OF ACAP
!          IN THE ORIGINAL PROGRAM THE DIMENSION OF ACAP WAS 7000.
!          FOR CONSERVING SPACE THIS SHOULD BE NOT MUCH HIGHER THAN
!          THE VALUE, N=1.1*(NREAL**2 + NIMAG**2)**.5 * X + 1
!      WVNO: 2*PIE / WAVELENGTH
!
!    ALSO THE SUBROUTINE COMPUTES THE CAPITAL A FUNCTION BY MAKING USE O
!    DOWNWARD RECURRENCE RELATIONSHIP.
!
!      TA(1): REAL PART OF WFN(1).  TA(2): IMAGINARY PART OF WFN(1).
!      TA(3): REAL PART OF WFN(2).  TA(4): IMAGINARY PART OF WFN(2).
!      TB(1): REAL PART OF FNA.     TB(2): IMAGINARY PART OF FNA.
!      TC(1): REAL PART OF FNB.     TC(2): IMAGINARY PART OF FNB.
!      TD(1): REAL PART OF FNAP.    TD(2): IMAGINARY PART OF FNAP.
!      TE(1): REAL PART OF FNBP.    TE(2): IMAGINARY PART OF FNBP.
!      FNAP, FNBP  ARE THE PRECEDING VALUES OF FNA, FNB RESPECTIVELY.
! **********************************************************************

!     /*--------------------------------------------------------------*/
!     /* Initially, make all types undefined.                         */
!     /*--------------------------------------------------------------*/

!      IMPLICIT UNDEFINED(A-Z)

!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/

      INTEGER*4  JX, IT, LL

      REAL*8     RO, RFR, RFI, THETD(IT), QEXT, QSCAT, CTBRQS,   &
                 ELTRMX(4,IT,2), PIE(3,IT), TAU(3,IT), CSTHT(IT),   &
                 SI2THT(IT), QBS, R,  RE2, TMAG2, WVNO

      COMPLEX*16 ACAP(LL)

!     /*--------------------------------------------------------*/
!     /* Variables used in the calculations below.              */
!     /*--------------------------------------------------------*/

      INTEGER*4  IFLAG, J, K, M, N, NN, NMX1, NMX2

      REAL*8     T(5), TA(4), TB(2), TC(2), TD(2), TE(2), X,   &
                 RX, X1, Y1, X4, Y4, SINX1, SINX4, COSX1, COSX4,   &
                 EY1, E2Y1, EY4, EY1MY4, EY1PY4, AA, BB,   &
                 CC, DD, DENOM, REALP, AMAGP, QBSR, QBSI, RMM,   &
                 PIG, RXP4

      COMPLEX*16 FNAP,      FNBP,      W,   &
                 FNA,       FNB,       RF,           RRF,   &
                 RRFX,      WM1,       FN1,          FN2,   &
                 TC1,       TC2,       WFN(2),       Z(4),   &
                 K1,        K2,        K3,   &
                 RC,        U(8),      DH1,   &
                 DH2,       DH4,       P24H24,       P24H21,   &
                 PSTORE,    HSTORE,    DUMMY,        DUMSQ

!     /*--------------------------------------------------------*/
!     /* Define the common block.                               */
!     /*--------------------------------------------------------*/

      COMMON / WARRAY / W(3,9000)

!
!      EQUIVALENCE   (FNA,TB(1)),(FNB,TC(1)),(FNAP,TD(1)),(FNBP,TE(1))
!
!   IF THE CORE IS SMALL SCATTERING IS COMPUTED FOR THE SHELL ONLY
!

!     /*--------------------------------------------------------*/
!     /* Begin the Mie calculations.                            */
!     /*--------------------------------------------------------*/

      IFLAG = 1
      IF ( R/RO .LT. 1.0D-06 )   IFLAG = 2
      IF ( JX .LE. IT )   GO TO 20
         WRITE( 6,7 )
         WRITE( 6,6 )
         STOP 30
   20 RF =  CMPLX( RFR,  -RFI )
      RC =  CMPLX( RE2,-TMAG2 )
      X  =  RO * WVNO
      K1 =  RC * WVNO
      K2 =  RF * WVNO
      K3 =  CMPLX( WVNO, 0.0D0 )
      Z(1) =  K2 * RO
      Z(2) =  K3 * RO
      Z(3) =  K1 * R
      Z(4) =  K2 * R
      X1   =  DREAL( Z(1) )
      Y1   =  DIMAG( Z(1) )
      X4   =  DREAL( Z(4) )
      Y4   =  DIMAG( Z(4) )
      RRF  =  1.0D0 / RF
      RX   =  1.0D0 / X
      RRFX =  RRF * RX
      T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
      T(1) =  DSQRT( T(1) )
      NMX1 =  1.30D0* T(1)
!
      IF ( NMX1 .LE. LL-1 )   GO TO 21
         WRITE(6,8)
         STOP 32
   21 NMX2 = T(1) * 1.2  ! JCB
      IF ( NMX1 .GT.  150 )   GO TO 22  ! JCB
         NMX1 = 150     ! JCB
        NMX2 = 135     ! JCB

!    21 NMX2 = T(1) * 1.5  ! JCB
!       IF ( NMX1 .GT.  200 )   GO TO 22  ! JCB
!         NMX1 = 200     ! JCB
        NMX2 = 180    ! JCB
!
   22 ACAP( NMX1+1 )  =  ( 0.0D0,0.0D0 )
      IF ( IFLAG .EQ. 2 )   GO TO 26
         DO 29   N = 1,3
   29    W( N,NMX1+1 )  =  ( 0.0D0,0.0D0 )
   26 CONTINUE
      DO 23   N = 1,NMX1
         NN = NMX1 - N + 1
         ACAP(NN) = (NN+1)*RRFX - 1.0D0 / ((NN+1)*RRFX + ACAP(NN+1))
         IF ( IFLAG .EQ. 2 )   GO TO 23
            DO 31   M = 1,3
   31       W( M,NN ) = (NN+1) / Z(M+1)  -   &
                         1.0D0 / ((NN+1) / Z(M+1) + W( M,NN+1 ))
   23 CONTINUE
!
      DO 30   J = 1,JX
      IF ( THETD(J) .LT. 0.0D0 )  THETD(J) =  DABS( THETD(J) )
      IF ( THETD(J) .GT. 0.0D0 )  GO TO 24
      CSTHT(J)  = 1.0D0
      SI2THT(J) = 0.0D0
      GO TO 30
   24 IF ( THETD(J) .GE. 90.0D0 )  GO TO 25
      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0D0
      CSTHT(J)  =  DCOS( T(1) )
      SI2THT(J) =  1.0D0 - CSTHT(J)**2
      GO TO 30
   25 IF ( THETD(J) .GT. 90.0 )  GO TO 28
      CSTHT(J)  =  0.0D0
      SI2THT(J) =  1.0D0
      GO TO 30
   28 WRITE( 6,5 )  THETD(J)
      WRITE( 6,6 )
      STOP 34
   30 CONTINUE
!
      DO 35  J = 1,JX
      PIE(1,J) =  0.0D0
      PIE(2,J) =  1.0D0
      TAU(1,J) =  0.0D0
      TAU(2,J) =  CSTHT(J)
   35 CONTINUE
!
! INITIALIZATION OF HOMOGENEOUS SPHERE
!
      T(1)   =  DCOS(X)
      T(2)   =  DSIN(X)
      WM1    =  CMPLX( T(1),-T(2) )
      WFN(1) =  CMPLX( T(2), T(1) )
      TA(1)  =  T(2)
      TA(2)  =  T(1)
      WFN(2) =  RX * WFN(1) - WM1
      TA(3)  =  DREAL(WFN(2))
      TA(4)  =  DIMAG(WFN(2))
!
      IF ( IFLAG .EQ. 2 )   GO TO 560
      N = 1
!
! INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
!
      SINX1   =  DSIN( X1 )
      SINX4   =  DSIN( X4 )
      COSX1   =  DCOS( X1 )
      COSX4   =  DCOS( X4 )
      EY1     =  DEXP( Y1 )
      E2Y1    =  EY1 * EY1
      EY4     =  DEXP( Y4 )
      EY1MY4  =  DEXP( Y1 - Y4 )
      EY1PY4  =  EY1 * EY4
      EY1MY4  =  DEXP( Y1 - Y4 )
      AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
      BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
      CC  =  SINX1 * ( E2Y1 + 1.0D0 )
      DD  =  COSX1 * ( E2Y1 - 1.0D0 )
      DENOM   =  1.0D0  +  E2Y1 * (4.0D0*SINX1*SINX1 - 2.0D0 + E2Y1)
      REALP   =  ( AA * CC  +  BB * DD ) / DENOM
      AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
      DUMMY   =  CMPLX( REALP, AMAGP )
      AA  =  SINX4 * SINX4 - 0.5D0
      BB  =  COSX4 * SINX4
      P24H24  =  0.5D0 + CMPLX( AA,BB ) * EY4 * EY4
      AA  =  SINX1 * SINX4  -  COSX1 * COSX4
      BB  =  SINX1 * COSX4  +  COSX1 * SINX4
      CC  =  SINX1 * SINX4  +  COSX1 * COSX4
      DD  = -SINX1 * COSX4  +  COSX1 * SINX4
      P24H21  =  0.5D0 * CMPLX( AA,BB ) * EY1 * EY4  +   &
                 0.5D0 * CMPLX( CC,DD ) * EY1MY4
      DH4  =  Z(4) / (1.0D0 + (0.0D0,1.0D0) * Z(4))  -  1.0D0 / Z(4)
      DH1  =  Z(1) / (1.0D0 + (0.0D0,1.0D0) * Z(1))  -  1.0D0 / Z(1)
      DH2  =  Z(2) / (1.0D0 + (0.0D0,1.0D0) * Z(2))  -  1.0D0 / Z(2)
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
!
! NOTE:  THE DEFINITIONS OF U(I) IN THIS PROGRAM ARE NOT THE SAME AS
!        THE USUBI DEFINED IN THE ARTICLE BY TOON AND ACKERMAN.  THE
!        CORRESPONDING TERMS ARE:
!          USUB1 = U(1)                       USUB2 = U(5)
!          USUB3 = U(7)                       USUB4 = DUMSQ
!          USUB5 = U(2)                       USUB6 = U(3)
!          USUB7 = U(6)                       USUB8 = U(4)
!          RATIO OF SPHERICAL BESSEL FTN TO SPHERICAL HENKAL FTN = U(8)
!
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0D0,-1.0D0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
!
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /   &
                     ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /   &
                     ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
!
!  Explicit equivalences added by J. Francis
      TB(1) = DREAL(FNA)
      TB(2) = DIMAG(FNA)
      TC(1) = DREAL(FNB)
      TC(2) = DIMAG(FNB)
      GO TO 561
  560 TC1  =  ACAP(1) * RRF  +  RX
      TC2  =  ACAP(1) * RF   +  RX
      FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
      FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
      TB(1) = DREAL(FNA)
      TB(2) = DIMAG(FNA)
      TC(1) = DREAL(FNB)
      TC(2) = DIMAG(FNB)
!
  561 CONTINUE
      FNAP = FNA
      FNBP = FNB
      TD(1) = DREAL(FNAP)
      TD(2) = DIMAG(FNAP)
      TE(1) = DREAL(FNBP)
      TE(2) = DIMAG(FNBP)
      T(1) = 1.50D0
!
!    FROM HERE TO THE STATMENT NUMBER 90, ELTRMX(I,J,K) HAS
!    FOLLOWING MEANING:
!    ELTRMX(1,J,K): REAL PART OF THE FIRST COMPLEX AMPLITUDE.
!    ELTRMX(2,J,K): IMAGINARY PART OF THE FIRST COMPLEX AMPLITUDE.
!    ELTRMX(3,J,K): REAL PART OF THE SECOND COMPLEX AMPLITUDE.
!    ELTRMX(4,J,K): IMAGINARY PART OF THE SECOND COMPLEX AMPLITUDE.
!    K = 1 : FOR THETD(J) AND K = 2 : FOR 180.0 - THETD(J)
!    DEFINITION OF THE COMPLEX AMPLITUDE: VAN DE HULST,P.125.
!
      TB(1) = T(1) * TB(1)
      TB(2) = T(1) * TB(2)
      TC(1) = T(1) * TC(1)
      TC(2) = T(1) * TC(2)
      DO 60 J = 1,JX
          ELTRMX(1,J,1) = TB(1) * PIE(2,J) + TC(1) * TAU(2,J)
          ELTRMX(2,J,1) = TB(2) * PIE(2,J) + TC(2) * TAU(2,J)
          ELTRMX(3,J,1) = TC(1) * PIE(2,J) + TB(1) * TAU(2,J)
          ELTRMX(4,J,1) = TC(2) * PIE(2,J) + TB(2) * TAU(2,J)
          ELTRMX(1,J,2) = TB(1) * PIE(2,J) - TC(1) * TAU(2,J)
          ELTRMX(2,J,2) = TB(2) * PIE(2,J) - TC(2) * TAU(2,J)
          ELTRMX(3,J,2) = TC(1) * PIE(2,J) - TB(1) * TAU(2,J)
          ELTRMX(4,J,2) = TC(2) * PIE(2,J) - TB(2) * TAU(2,J)
   60 CONTINUE
!
      QEXT   = 2.0D0 * ( TB(1) + TC(1))
      QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75D0
      CTBRQS = 0.0D0
      QBSR   = -2.0D0*(TC(1) - TB(1))
      QBSI   = -2.0D0*(TC(2) - TB(2))
      RMM    = -1.0D0
      N = 2
   65 T(1) = 2*N - 1
      T(2) =   N - 1
      T(3) = 2*N + 1
      DO 70  J = 1,JX
          PIE(3,J) = ( T(1)*PIE(2,J)*CSTHT(J) - N*PIE(1,J) ) / T(2)
          TAU(3,J) = CSTHT(J) * ( PIE(3,J) - PIE(1,J) ) -   &
                                T(1)*SI2THT(J)*PIE(2,J) + TAU(1,J)
   70 CONTINUE
!
! HERE SET UP HOMOGENEOUS SPHERE
!
      WM1    =  WFN(1)
      WFN(1) =  WFN(2)
      TA(1)  =  DREAL(WFN(1))
      TA(2)  =  DIMAG(WFN(1))
      WFN(2) =  T(1) * RX * WFN(1)  -  WM1
      TA(3)  =  DREAL(WFN(2))
      TA(4)  =  DIMAG(WFN(2))
!
      IF ( IFLAG .EQ. 2 )   GO TO 1000
!
! HERE SET UP STRATIFIED SPHERE
!
      DH2  =  - N / Z(2)  +  1.0D0 / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0D0 / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0D0 / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
!
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0D0,-1.0D0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
!
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /   &
                     ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /   &
                     ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
      TB(1) = DREAL(FNA)
      TB(2) = DIMAG(FNA)
      TC(1) = DREAL(FNB)
      TC(2) = DIMAG(FNB)
!
 1000 CONTINUE
      TC1  =  ACAP(N) * RRF  +  N * RX
      TC2  =  ACAP(N) * RF   +  N * RX
      FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
      FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
      M    =  WVNO * R
      IF ( N .LT. M )   GO TO 1002
      IF ( IFLAG .EQ. 2 )   GO TO 1001
      IF ( (ABS(  ( FN1-FNA ) / FN1  )  .LT.  1.0D-09).AND.(ABS(  ( FN2-FNB ) / FN2  )  .LT. 1.0D-09)) IFLAG = 2
      IF ( IFLAG .EQ. 1 )   GO TO 1002
 1001 FNA  =  FN1
      FNB  =  FN2
      TB(1) = DREAL(FNA)
      TB(2) = DIMAG(FNA)
      TC(1) = DREAL(FNB)
      TC(2) = DIMAG(FNB)
!
 1002 CONTINUE
      T(5)  =  N
      T(4)  =  T(1) / ( T(5) * T(2) )
      T(2)  =  (  T(2) * ( T(5) + 1.0D0 )  ) / T(5)
!
      CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2)   &
                         +           TE(1) * TC(1)  +  TE(2) * TC(2) )   &
                         +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
      QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
      T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT   =  QSCAT  +  T(3) * T(4)
      RMM     =  -RMM
      QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
      QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
!
      T(2)    =  N * (N+1)
      T(1)    =  T(3) / T(2)
      K = (N/2)*2
      DO 80 J = 1,JX
       ELTRMX(1,J,1)=ELTRMX(1,J,1)+T(1)*(TB(1)*PIE(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1)=ELTRMX(2,J,1)+T(1)*(TB(2)*PIE(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1)=ELTRMX(3,J,1)+T(1)*(TC(1)*PIE(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1)=ELTRMX(4,J,1)+T(1)*(TC(2)*PIE(3,J)+TB(2)*TAU(3,J))
      IF ( K .EQ. N )  THEN
       ELTRMX(1,J,2)=ELTRMX(1,J,2)+T(1)*(-TB(1)*PIE(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(-TB(2)*PIE(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(-TC(1)*PIE(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(-TC(2)*PIE(3,J)+TB(2)*TAU(3,J))
      ELSE
       ELTRMX(1,J,2)=ELTRMX(1,J,2)+T(1)*(TB(1)*PIE(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2)=ELTRMX(2,J,2)+T(1)*(TB(2)*PIE(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2)=ELTRMX(3,J,2)+T(1)*(TC(1)*PIE(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2)=ELTRMX(4,J,2)+T(1)*(TC(2)*PIE(3,J)-TB(2)*TAU(3,J))
      END IF
   80 CONTINUE
!
      IF ( T(4) .LT. 1.0D-14 )   GO TO 100
      N = N + 1
      DO 90 J = 1,JX
         PIE(1,J)  =   PIE(2,J)
         PIE(2,J)  =   PIE(3,J)
         TAU(1,J)  =  TAU(2,J)
         TAU(2,J)  =  TAU(3,J)
   90 CONTINUE
      FNAP  =  FNA
      FNBP  =  FNB
      TD(1) = DREAL(FNAP)
      TD(2) = DIMAG(FNAP)
      TE(1) = DREAL(FNBP)
      TE(2) = DIMAG(FNBP)
      IF ( N .LE. NMX2 )   GO TO 65
         WRITE( 6,8 )
         STOP 36
  100 CONTINUE
!      DO 120 J = 1,JX
!      DO 120 K = 1,2
!         DO  115  I= 1,4
!         T(I)  =  ELTRMX(I,J,K)
!  115    CONTINUE
!         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
!         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
!         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
!         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
!  120 CONTINUE
      T(1)    =  2.0D0 * RX**2
      QEXT    =   QEXT * T(1)
      QSCAT   =  QSCAT * T(1)
      CTBRQS  =  2.0D0 * CTBRQS * T(1)
!
! QBS IS THE BACK SCATTER CROSS SECTION
!
      PIG   = DACOS(-1.0D0)
      RXP4  = RX*RX/(4.0D0*PIG)
      QBS   = RXP4*(QBSR**2 + QBSI**2)
!
    
!    5 FORMAT(10X,' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN',   &
!       '90.0 DEGREES. IT IS ', E15.4 )
!    6 FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
!    7 FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
!    8 FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST',   &
!      ' GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )



 5    FORMAT(' THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN 90.0 DEGREES. IT IS ', E15.4 )
 6    FORMAT('PLEASE READ COMMENTS.')
 7    FORMAT('THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
 8    FORMAT('THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST GET DETAILED OUTPUT AND MODIFY SUBROUTINE')

!
      RETURN
      END SUBROUTINE DMIESS
! /*****************************************************************/
! /* SUBROUTINE ACKMIEPARICLE                                      */
! /*****************************************************************/

! THIS PROGRAM COMPUTES SCATTERING PROPERTIES FOR DISTRIBUTIONS OF
! PARTICLES COMPOSED OF A CORE OF ONE MATERIAL AND A SHELL OF ANOTHER.

! /*---------------------------------------------------------------*/
! /* INPUTS:                                                       */
! /*---------------------------------------------------------------*/

!   VLAMBc: Wavelength of the radiation
!   NRGFLAGc: Flag to indicate a number density of volume radius
!   RGc: Geometric mean radius of the particle distribution
!   SIGMAGc: Geometric standard deviation of the distribution
!   SHELRc: Real part of the index of refraction for the shell
!   SHELIc: Imaginary part of the index of refraction for the shell
!   RINc: Inner core radius as a fraction of outer shell radius
!   CORERc: Real part of the index of refraction for the core
!   COREIc: Imaginary part of the index of refraction for the core
!   NANG: Number of scattering angles between 0 and 90 degrees,
!         inclusive

! /*---------------------------------------------------------------*/
! /* OUTPUTS:                                                      */
! /*---------------------------------------------------------------*/

!   QEXTc: Extinction efficiency of the particle
!   QSCAc: Scattering efficiency of the particle
!   QBACKc: Backscatter efficiency of the particle
!   EXTc: Extinction cross section of the particle
!   SCAc: Scattering cross section of the particle
!   BACKc: Backscatter cross section of the particle
!   GSCA: Asymmetry parameter of the particle's phase function
!   ANGLES(NAN): Scattering angles in degrees
!   S1R(NAN): Real part of the amplitude scattering matrix
!   S1C(NAN): Complex part of the amplitude scattering matrix
!   S2R(NAN): Real part of the amplitude scattering matrix
!   S2C(NAN): Complex part of the amplitude scattering matrix
!   S11N: Normalization coefficient of the scattering matrix
!   S11(NAN): S11 scattering coefficients
!   S12(NAN): S12 scattering coefficients
!   S33(NAN): S33 scattering coefficients
!   S34(NAN): S34 scattering coefficients
!   SPOL(NAN): Degree of polarization of unpolarized, incident light
!   SP(NAN): Phase function
!
! NOTE: NAN=2*NANG-1 is the number of scattering angles between
!       0 and 180 degrees, inclusive.
! /*---------------------------------------------------------------*/

      SUBROUTINE ACKMIEPARTICLE( VLAMBc,NRGFLAGc,RGcmin,RGcmax,   &
                   RGc,SIGMAGc,SHELRc,   &
                   SHELIc, RINc,CORERc,COREIc,NANG,QEXTc,QSCATc,   &
                   QBACKc, EXTc,SCATc,BACKc, GSCAc,   &
                ANGLESc,S1R,S1C,S2R,S2C,S11N,S11,S12,S33,S34,SPOL,SP)

!     /*--------------------------------------------------------*/
!     /* Parameter statements.                                  */
!     /*--------------------------------------------------------*/

      IMPLICIT REAL*8 (A-H, O-Z)
      IMPLICIT INTEGER (M-N)

      PARAMETER(MXNANG=501)
      !INTEGER, PARAMETER :: MAXNANG=501

!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/


!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/

      REAL*8 VLAMBc,RGcmin,RGcmax,RGc,SIGMAGc,SHELRc,SHELIc
      REAL*8 RINc,CORERc,COREIc
      INTEGER*4 NANG
      REAL*8 QEXTc,QSCATc,QBACKc,EXTc,SCATc,BACKc,GSCAc
      REAL*8 ANGLESc(*),S1R(*),S1C(*),S2R(*),S2C(*)
      REAL*8 S11N,S11(*),S12(*),S33(*),S34(*),SPOL(*),SP(*)

!     /*--------------------------------------------------------*/
!     /* Define the types of the common block.                  */
!     /*--------------------------------------------------------*/

      INTEGER*4 IPHASE

      REAL*8 ALAMB, RGmin, RGmax, RGV, SIGMAG, RGCFRAC, RFRS,   &
             RFIS, RFRC, RFIC

!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/

!      IMPLICIT REAL*8 (A-H, O-Z)

!     /*--------------------------------------------------------*/
!     /* Input common block for scattering calculations.        */
!     /*--------------------------------------------------------*/

      COMMON / PHASE  / IPHASE

      COMMON / INPUTS / ALAMB, RGmin, RGmax, RGV, SIGMAG,   &
                        RGCFRAC, RFRS,  RFIS,  RFRC,  RFIC

!     /*--------------------------------------------------------*/
!     /* Output common block for scattering calculations.       */
!     /*--------------------------------------------------------*/

      COMMON / OUTPUTS / QEXT, QSCAT, QBS, EXT, SCAT, BSCAT, ASY

!     /*--------------------------------------------------------*/
!     /* Arrays to hold the results of the scattering and       */
!     /* moment calculations.                                   */
!     /*--------------------------------------------------------*/

      DIMENSION COSPHI(2*MXNANG-1), SCTPHS(2*MXNANG-1)

!     /*--------------------------------------------------------*/
!     /* Copy the input parameters into the common block INPUTS */
!     /*--------------------------------------------------------*/

         NSCATH  = NANG

         IPHASE  = 0
         ALAMB   = VLAMBc
         RGmin   = RGcmin
         RGmax   = RGcmax
         RGV     = RGc
         SIGMAG  = SIGMAGc
         RGCFRAC = RINc
         RFRS    = SHELRc
         RFIS    = SHELIc
         RFRC    = CORERc
         RFIC    = COREIc

!     /*--------------------------------------------------------*/
!     /* Calculate the particle scattering properties for the   */
!     /* given wavelength, particle distribution and indices of */
!     /* refraction of inner and outer material.                */
!     /*--------------------------------------------------------*/

         CALL PFCNPARTICLE(NSCATH, COSPHI, SCTPHS,   &
             ANGLESc,S1R,S1C,S2R,S2C,S11N,S11,S12,S33,S34,SPOL,SP)

!     /*--------------------------------------------------------*/
!     /* If IPHASE = 1, then the full phase function was        */
!     /* calculated; now, go calculate its moments.             */
!     /*--------------------------------------------------------*/

!        IF (IPHASE .gt. 0) CALL DISMOM (NSCATA,COSPHI,SCTPHS,RMOMS)

!     /*--------------------------------------------------------*/
!     /* Copy the variables in the common block OUTPUTS to the  */
!     /* variable addresses passed into this routine.           */
!     /*--------------------------------------------------------*/

         QEXTc  = QEXT
         QSCATc = QSCAT
         QBACKc = QBS

         EXTc  = EXT
         SCATc = SCAT
         BACKc = BSCAT

         GSCAc  = ASY

!     /*--------------------------------------------------------*/
!     /* FORMAT statements.                                     */
!     /*--------------------------------------------------------*/

 107   FORMAT(I6, 'IS AN INVALID MEAN RADIUS FLAG')

!     /*--------------------------------------------------------*/
!     /* DONE with this subroutine so exit.                     */
!     /*--------------------------------------------------------*/

      END SUBROUTINE ACKMIEPARTICLE

! /*****************************************************************/
! /* SUBROUTINE PFCNPARTICLE                                       */
! /*****************************************************************/
!
! THIS SUBROUTINE COMPUTES THE PHASE FUNCTION SCTPHS(I) AT NSCATA
! ANGLES BETWEEN 0.0 AND 180.0 DEGREES SPECIFIED BY COSPHI(I) WHICH
! CONTAINS THE ANGLE COSINES. THESE VALUES ARE RETURNED TO FUNCTION
! PHASFN FOR AZIMUTHAL AVERAGING.
! INPUT DATA FOR THIS ROUTINE IS PASSED THROUGH COMMON /SIZDIS/
! AND CONSISTS OF THE FOLLOWING PARAMETERS
! NEWSD  = 1 IF SIZE DIS VALUES HAVE NOT PREVIOUSLY BEEN USED IN
!          THIS ROUTINE,  = 0 OTHERWISE.
! RGV    = GEOMETRIC MEAN RADIUS FOR THE VOLUME DISTRIBUTION OF THE
!          SPECIFIED PARTICLES
! SIGMAG = GEOMETRIC STANDARD DEVIATION
! RFR,I  = REAL AND IMAGINARY INDEX OF REFRACTION OF PARTICLES
! ALAMB  = WAVELENGTH AT WHICH CALCULATIONS ARE TO BE PERFORMED
!
! /*---------------------------------------------------------------*/

      SUBROUTINE PFCNPARTICLE( NSCATH, COSPHI, SCTPHS,   &
           ANGLESc,S1R,S1C,S2R,S2C,S11N,S11,S12,S33,S34,SPOL,SP )

!     /*--------------------------------------------------------*/
!     /* Parameter statements.                                  */
!     /*--------------------------------------------------------*/

      IMPLICIT REAL*8 (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER(MXNANG=501, MXNWORK=500000)

!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/

!      IMPLICIT REAL*8 (A-H, O-Z)

!     /*--------------------------------------------------------*/
!     /* Dimension statements.                                  */
!     /*--------------------------------------------------------*/

      REAL*8 ANGLESc(*),S1R(*),S1C(*),S2R(*),S2C(*)
      REAL*8 S11N,S11(*),S12(*),S33(*),S34(*),SPOL(*),SP(*)

!     /*--------------------------------------------------------*/
!     /* Define the types of the common block.                  */
!     /*--------------------------------------------------------*/

      INTEGER*4 IPHASE

      REAL*8 ALAMB, RGmin, RGmax, RGV, SIGMAG, RGCFRAC, RFRS,   &
             RFIS, RFRC, RFIC

!     /*--------------------------------------------------------*/
!     /* Set reals to 8 bytes, i.e., double precision.          */
!     /*--------------------------------------------------------*/

!      IMPLICIT REAL*8 (A-H, O-Z)

!     /*--------------------------------------------------------*/
!     /* Input common block for scattering calculations.        */
!     /*--------------------------------------------------------*/

      COMMON / PHASE  / IPHASE

      COMMON / INPUTS / ALAMB, RGmin, RGmax, RGV, SIGMAG,   &
                        RGCFRAC, RFRS,  RFIS,  RFRC,  RFIC

!     /*--------------------------------------------------------*/
!     /* Output common block for scattering calculations.       */
!     /*--------------------------------------------------------*/

      COMMON / OUTPUTS / QEXT, QSCAT, QBS, EXT, SCAT, BSCAT, ASY

!     /*--------------------------------------------------------*/
!     /* Arrays to perform the scattering calculations and to   */
!     /* hold the subsequent results.                           */
!     /*--------------------------------------------------------*/

      REAL*8  THETA(MXNANG), ELTRMX(4,MXNANG,2), PII(3,MXNANG),   &
              TAU(3,MXNANG), CSTHT(MXNANG),      SI2THT(MXNANG)

      REAL*8   ROUT, RFRO, RFIO, DQEXT, DQSCAT, CTBRQS, DQBS,   &
               RIN,  RFRI, RFII, WNUM

      COMPLEX*16  ACAP(MXNWORK)

      DIMENSION  COSPHI(2*MXNANG-1), SCTPHS(2*MXNANG-1)
      INTEGER idum

!     /*--------------------------------------------------------*/
!     /* Obvious variable initializations.                      */
!     /*--------------------------------------------------------*/

         PIE    = DACOS( -1.0D0 )

!     /*--------------------------------------------------------*/
!     /* Maximum number of scattering angles between 0 and 90   */
!     /* degrees, inclusive.                                    */
!     /*--------------------------------------------------------*/

         IT = MXNANG

   DO idum = 1, IT
           THETA(idum) = 0.0
         ENDDO

!     /*--------------------------------------------------------*/
!     /* Maximum number of scattering angles between 0 and 180  */
!     /* degrees, inclusive.                                    */
!     /*--------------------------------------------------------*/

         IT2 = 2 * IT - 1

!     /*--------------------------------------------------------*/
!     /* Dimension of the work array ACAP.                      */
!     /*--------------------------------------------------------*/

         LL = MXNWORK

!     /*--------------------------------------------------------*/
!     /* NSCATA is the actual user-requested number of          */
!     /* scattering angles between 0 and 90 degrees, inclusive. */
!     /*--------------------------------------------------------*/

         NSCATA = 2 * NSCATH - 1

!     /*--------------------------------------------------------*/
!     /* If the user did not request a phase function, then we  */
!     /* can set NSCATA and NSCATH to 0.                        */
!     /*--------------------------------------------------------*/

         IF ( IPHASE  .le.  0 )  then
              NSCATH  =  0
              NSCATA  =  0
         ENDIF

!     /*--------------------------------------------------------*/
!     /* Check to make sure that the user-requested number of   */
!     /* scattering angles does not excede the current maximum  */
!     /* limit.                                                 */
!     /*--------------------------------------------------------*/

         IF ( NSCATA .gt. IT2  .OR.  NSCATH .gt. IT)  then
              WRITE( 6,105 )  NSCATA, NSCATH, IT2, IT
              STOP 11
         ENDIF

!     /*--------------------------------------------------------*/
!     /* Subroutine SCATANGLES was added by EEC[0495] in order  */
!     /* to facilitate changing the scattering angle locations  */
!     /* output by the Ackerman and Toon Mie code.              */
!     /*--------------------------------------------------------*/

!         CALL SCATANGLES(NSCATH,THETA,COSPHI)

!     /*--------------------------------------------------------*/
!     /* COMPUTE SCATTERING PROPERTIES OF THE PARTICLE.         */
!     /*--------------------------------------------------------*/

!     /*--------------------------------------------------------*/
!     /* DMIESS expects a wavenumber.                           */
!     /*--------------------------------------------------------*/

      WNUM = (2.D0*PIE) / ALAMB

!     /*--------------------------------------------------------*/
!     /* DMIESS assignments of the indices of refraction of the */
!     /* core and shell materials.                              */
!     /*--------------------------------------------------------*/

      RFRO = RFRS
      RFIO = RFIS
      RFRI = RFRC
      RFII = RFIC

!     /*--------------------------------------------------------*/
!     /* DMIESS core and shell radii.                           */
!     /*--------------------------------------------------------*/

       ROUT = RGV
       RIN  = RGCFRAC * ROUT

!     /*--------------------------------------------------------*/
!     /* Scattering angles are symmetric about 90 degrees.      */
!     /*--------------------------------------------------------*/

      IF ( NSCATH  .eq.  0.0 )  THEN
           JX  =  1
      ELSE
           JX   = NSCATH
      ENDIF

!     /*--------------------------------------------------------*/
!     /* Compute the scattering properties for this particle.   */
!     /*--------------------------------------------------------*/

      CALL DMIESS(  ROUT,    RFRO,    RFIO,    THETA,   JX,   &
                    DQEXT,   DQSCAT,  CTBRQS,  ELTRMX,  PII,   &
                    TAU,     CSTHT,   SI2THT,  ACAP,    DQBS,  IT,   &
                    LL,      RIN,     RFRI,    RFII,    WNUM   )

!     /*--------------------------------------------------------*/
!     /* Compute total cross-sectional area of the particle.    */
!     /*--------------------------------------------------------*/

      X = PIE * RGV * RGV

!     /*--------------------------------------------------------*/
!     /* Assign the final extinction efficiency.                */
!     /*--------------------------------------------------------*/

      QEXT = DQEXT

!     /*--------------------------------------------------------*/
!     /* Compute total extinction cross-section due to particle.*/
!     /*--------------------------------------------------------*/

      EXT = DQEXT * X

!     /*--------------------------------------------------------*/
!     /* Assign the final scattering efficiency.                */
!     /*--------------------------------------------------------*/

      QSCAT = DQSCAT

!     /*--------------------------------------------------------*/
!     /* Compute total scattering cross-section due to particle.*/
!     /*--------------------------------------------------------*/

      SCAT = DQSCAT * X

!     /*--------------------------------------------------------*/
!     /* Assign the final backscatter efficiency.               */
!     /*--------------------------------------------------------*/

      QBS = DQBS

!     /*--------------------------------------------------------*/
!     /* Compute backscatter due to particle.                   */
!     /*--------------------------------------------------------*/

      BSCAT = DQBS * X

!     /*--------------------------------------------------------*/
!     /* Compute asymmetry parameter due to particle.           */
!     /*--------------------------------------------------------*/

      ASY = (CTBRQS * X) / SCAT

!     /*--------------------------------------------------------*/
!     /* If IPHASE is 1, compute the phase function.            */
!     /* S33 and S34 matrix elements are normalized by S11. S11 */
!     /* is normalized to 1.0 in the forward direction.  The    */
!     /* variable SPOL is the degree of polarization for        */
!     /* incident unpolarized light.                            */
!     /*--------------------------------------------------------*/

      IF ( IPHASE  .gt.  0 )  THEN

        DO 355 J=1,NSCATA

          IF (J .LE. JX)  THEN
             JJ = J
             NINDEX = 1
          ELSE
             JJ = NSCATA - J + 1
             NINDEX = 2
          ENDIF

          ANGLESc(J) = COSPHI(J)

          S1R(J) = ELTRMX(1,JJ,NINDEX)
          S1C(J) = ELTRMX(2,JJ,NINDEX)
          S2R(J) = ELTRMX(3,JJ,NINDEX)
          S2C(J) = ELTRMX(4,JJ,NINDEX)

          S11(J) = 0.5D0*(S1R(J)**2+S1C(J)**2+S2R(J)**2+S2C(J)**2)
          S12(J) = 0.5D0*(S2R(J)**2+S2C(J)**2-S1R(J)**2-S1C(J)**2)
          S33(J) = S2R(J)*S1R(J) + S2C(J)*S1C(J)
          S34(J) = S2R(J)*S1C(J) - S1R(J)*S2C(J)

          SPOL(J) = -S12(J) / S11(J)

          SP(J) = (4.D0*PIE)*(S11(J) / (SCAT*WNUM**2))

  355   CONTINUE

!        /*-----------------------------------------------------*/
!        /* DONE with the phase function so exit the IF.        */
!        /*-----------------------------------------------------*/

      ENDIF

!     /*--------------------------------------------------------*/
!     /* END of the computations so exit the routine.           */
!     /*--------------------------------------------------------*/

      RETURN

!     /*--------------------------------------------------------*/
!     /* FORMAT statements.                                     */
!     /*--------------------------------------------------------*/

 100  FORMAT(7X, I3)
 105  FORMAT(1X,'NUMBER OF ANGLES SPECIFIED =',2I6, /   &
                   10X,'EXCEEDS ARRAY DIMENSIONS =',2I6 )

 120  FORMAT(/10X,'INTEGRATED VOLUME',           T40,'=',1PE14.5,/   &
               15X,'PERCENT VOLUME IN CORE',      T40,'=',0PF10.5,/   &
               15X,'PERCENT VOLUME IN SHELL',     T40,'=',0PF10.5,/   &
               10X,'INTEGRATED SURFACE AREA',     T40,'=',1PE14.5,/   &
               10X,'INTEGRATED NUMBER DENSITY',   T40,'=',1PE14.5 )
 125  FORMAT(10X,'CORE RADIUS COMPUTED FROM :', /, 20X, 9A8, /  )

 150  FORMAT(1X,'* * * WARNING * * *', /   &
               10X,'PHASE FUNCTION CALCULATION MAY NOT HAVE CONVERGED'/   &
               10X,'VALUES OF S1 AT NSDI-1 AND NSDI ARE :', 2E14.6, /   &
               10X,'VALUE OF X AT NSDI =', E14.6)

!     /*--------------------------------------------------------*/
!     /* DONE with this subroutine so exit.                     */
!     /*--------------------------------------------------------*/

      END SUBROUTINE PFCNPARTICLE

! /*****************************************************************/
! /*****************************************************************/
      SUBROUTINE SCATANGLESCONSTSINE(NSCATH,THETA,COSPHI)

      IMPLICIT REAL*8 (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      INTEGER*4 NSCATH
      REAL*8    THETA(*), COSPHI(*)

!      IMPLICIT REAL*8 (A-H, O-Z)

      PIE = DACOS(-1.D0)

      NSCATA = 2*NSCATH - 1

      DMU = ( 1.0D0 ) / ( NSCATH - 1.0D0 )

      DO  20  I=1,NSCATH
         SINANG    = ( I-1.0D0 ) * DMU
         ANG       = DASIN(SINANG)
         THETA(I)  = 180.0D0 * ANG / PIE
         COSPHI(I) = DCOS(ANG)
         J         = NSCATA - I + 1
         COSPHI(J) = -COSPHI(I)
   20 CONTINUE

      RETURN

      END SUBROUTINE SCATANGLESCONSTSINE
!**********************************************


end module dmiess_module











