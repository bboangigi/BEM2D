!! ---------------------------------------------------------------------
Module SOLVE_M
!! ---------------------------------------------------------------------
     Use glb2DBEM
!! ---------------------------------------------------------------------
      Implicit None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
      SUBROUTINE SOLVE(NB, NT, MX, NP, NQ, NEQ, AKB, SML, XP, YP, XQ, YQ, VN)
      Implicit None
!! -----------Input variables from app2DBEM_OSAKA-----------------------
      Integer, Intent(in)          :: NB, NT, MX, NP, NQ, NEQ     ! from app2DBEM_OSAKA
      Double Precision, Intent(in) :: AKB, SML                    ! from app2DBEM_OSAKA
!! -----------Variables declaration in SOLVE subroutine-----------------
      Integer :: I, J, K, M

      !! Variables declaration in SDSUB subroutine --------------------
      Double Precision :: XPI, YPI

!! -----------Input variables&arrays from OFFSET_M-----------------------------
       Double precision, Allocatable, Intent(in) :: XP(:), YP(:), XQ(:), YQ(:)
       Double Precision, Allocatable, Intent(in) :: VN(:,:)
!       Allocate(XP(MX), YP(MX), XQ(NQ), YQ(NQ))
!       Allocate(VN(3,NP))

!! ----------Arrays declaration in SOLVE subroutine------------------

       Complex(8), Allocatable :: ZAA(:,:), ZBB(:,:), ZSA(:,:), ZSB(:,:), ZFI(:,:)
       Complex(8), Allocatable :: ZS(:), ZD(:)
       Double Precision, Allocatable :: SS(:), DD(:)

       Allocate(ZSA(MX,NP), ZSB(MX,NEQ), ZAA(NP,NP), ZBB(NP,NEQ), ZFI(NEQ,NP))
       Allocate(ZS(NP), ZD(NP), SS(NP), DD(NP))

!       Complex(8) :: ZAA, ZSA, Z0, ZI, ZFI, ZBB, ZSB, ZS, ZD
!       Dimension :: ZSA(MX,NP), ZSB(MX,NEQ), ZAA(NP,NP), ZBB(NP,NEQ), ZFI(NEQ,NP)
! -------------------------------------------------
      DO 10 I=1,NB
      DO 20 J=1,NB
      20 ZAA(I,J)=(0.0d0)
      DO 10 M=1,NEQ
      ZBB(I,M)=(0.0d0)
      10 CONTINUE

      DO 30 I=1,NT
      DO 40 J=1,NB
      40 ZSA(I,J)=(0.0d0)
      DO 50 M=1,NEQ
      50 ZSB(I,M)=(0.0d0)
      IF(I.LE.NB) ZSA(I,I)=DCMPLX(PI,0.0D0)
      30 CONTINUE
!! -------------------------------------------------
      DO I=1,NT

            XPI = XP(I)
            YPI = YP(I)

            Call SDSUB(MX, NP, NEQ, NQ, NB, XP, YP, XQ, YQ, VN, XPI, YPI, SS, DD)
            Call SDCAL(MX, NP, NEQ, NQ, NB, AKB, XPI, YPI, XQ, YQ, ZS, ZD)

            DO J=1,NB
                  ZSA(I,J)=ZSA(I,J)+DD(J)+ZD(J)
            End Do

            DO M=1,3
                  DO J=1,NB
                        ZSB(I,M)=ZSB(I,M)+(SS(J)+ZS(J))*VN(M,J)
                  End Do
            End do

            ZSB(I,4)=PI2*ZEXP(-AKB*(YP(I)-ZI*XP(I)))

      End Do

      DO I=1,NB
            DO J=1,NB
            DO K=1,NT
                  ZAA(I,J)=ZAA(I,J)+ZSA(K,I)*ZSA(K,J)
            End do
            end do

            DO M=1,NEQ
            DO K=1,NT
                  ZBB(I,M)=ZBB(I,M)+ZSA(K,I)*ZSB(K,M)
            End do
            End do
      End do

      CALL ZSWEEP(NP,NB,ZAA,ZBB,NEQ,SML)

      IF(CDABS(ZAA(1,1)).LT.SML) WRITE(6,600)
      600 FORMAT(//10X,'*** ERROR: ZSWEEP IN SUBROUTINE (SOLVE)', &
            ' WAS ABNORMALLY DONE.',/23X,'PLEASE CHECK!'///)

      DO M=1,NEQ
      DO I=1,NB
            ZFI(M,I)=ZBB(I,M)
      End do
      End do

      ! write(*,*)
      ! DO I=1,NB
      !     write(*,"(i5,99(1pe12.4))") i, (ZFI(M,I), m = 1, NEQ)
      ! End do
      ! read(*,*)

      RETURN
      END Subroutine
!
! !! ---------------------------------------------------------------------
!!! ---------------------------------------------------------------------
      SUBROUTINE SDSUB(MX, NP, NEQ, NQ, NB, XP, YP, XQ, YQ, VN, XPI, YPI, SS, DD)
      !     Kernel Function: Rankine Source

      ! Implicit None
      ! IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!! -----------Input variables from app2DBEM_OSAKA------------------------------
      Integer :: MX, NP, NEQ, NQ, NB
      Double Precision :: AKB
!! -----------Input variables in SDSUB-----------------------------------------
      Integer :: J, L
      Double Precision :: XPI, YPI
      Double Precision :: DX, DY, D, SL, CDEL, SDEL, COEF, XA, XB, YA, YB, SUBA, SUBB, WA1, WA2, WA3, ABSC, SWA, DWA
!! -----------Input variables&arrays from OFFSET_M-----------------------------
      Double Precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:)
      Double Precision, Allocatable :: VN(:,:)
!! ----------Arrays declaration in SDSUB subroutine----------------------------
      Double Precision, Allocatable :: SS(:), DD(:)

      DO 100 J=1,NB
            SWA=0.0D0
            DWA=0.0D0

            IF(DABS(YPI).LT.1.0D-8) GOTO 10
            DX=XQ(J+1)-XQ(J)
            DY=YQ(J+1)-YQ(J)
            D=DSQRT(DX*DX+DY*DY)
            CDEL=DX/D
            SDEL=DY/D
            XA=XPI-XQ(J  )
            XB=XPI-XQ(J+1)

            SL=-1.0D0
            DO 200 L=1,2
            SL=-SL
            YA=SL*YPI-YQ(J  )
            YB=SL*YPI-YQ(J+1)
            SUBA=XA*CDEL+YA*SDEL
            SUBB=XB*CDEL+YB*SDEL
            COEF=XA*SDEL-YA*CDEL
            ABSC=DABS(COEF)
            WA1=0.5D0*(SUBB*DLOG(XB*XB+YB*YB)-SUBA*DLOG(XA*XA+YA*YA))
            IF(ABSC.LT.1.0D-10) THEN
            WA2=0.0D0
            WA3=0.0D0
            ELSE
            WA2=ABSC*(DATAN(SUBB/ABSC)-DATAN(SUBA/ABSC))
            WA3=WA2/COEF
            ENDIF
            SWA=SWA-(WA1+WA2)*SL
            DWA=DWA+ WA3*SL
      200 CONTINUE

      10 SS(J)=SWA
            DD(J)=DWA
      100 CONTINUE
            RETURN
            END Subroutine
!!! ---------------------------------------------------------------------
!! ---------------------------------------------------------------------
      SUBROUTINE SDCAL(MX, NP, NEQ, NQ, NB, AKB, XPI, YPI, XQ, YQ, ZS, ZD)
      !     Kernel Function: Wave Term
!! -----------Input variables from app2DBEM_OSAKA------------------------------
      Integer :: MX, NP, NEQ, NQ, NB
      Double Precision :: AKB
!! -----------Input variables in SDSUB-----------------------------------------
      Integer :: J, L
      Double Precision :: XPI, YPI
      Double Precision, Allocatable :: XQ(:), YQ(:)
!! -----------Input variables in SDCAL-----------------------------------------
      Double Precision :: XX, YY, XE, YE, RFL1, RFT1, RFL2, RFT2, EC, ES, SUB
      Double Precision :: DX, DY, D, CDEL, SDEL, SGNX
      Complex(8) :: ZETA, ZFC1, ZFC2, ZFS1, ZFS2, ZSUB
!! -----------Input variables&arrays from OFFSET_M-----------------------------
      Double Precision, Allocatable :: XP(:), YP(:)
      Double Precision, Allocatable :: VN(:,:)
!! -----------Input variables in SDSUB-----------------------------------------
      Double Precision, Allocatable :: SS(:), DD(:)

!! ----------Arrays declaration in SDCAL subroutine----------------------------
      Complex(8), Allocatable :: ZS(:), ZD(:)
      ! DIMENSION ZS(NB),ZD(NB)

      DO 100 J=1,NB
            ZS(J) = DCMPLX(0.0D0, 0.d0)
            ZD(J) = DCMPLX(0.0D0, 0.d0)
      100 CONTINUE

      XX=XPI-XQ(1)
      YY=YPI+YQ(1)

      SGNX=DSIGN(1.0D0,XX)
      IF(DABS(XX).LT.1.0D-10) SGNX=0.0D0
      XE=-AKB*YY
      YE=-AKB*DABS(XX)
      ZETA=CMPLX(XE,YE)

            CALL EZE1Z(XE,YE,EC,ES)
            RFL1=0.5D0*DLOG(XX**2+YY**2)
            RFT1=DATAN2(YY,XX)
            ZFC1= EC-PI*ZEXP(ZETA)*ZI
            ZFS1=(ES-PI*ZEXP(ZETA))*SGNX

            DO 200 J=1,NB
            XX=XPI-XQ(J+1)
            YY=YPI+YQ(J+1)
            SGNX=DSIGN(1.0D0,XX)
            IF(DABS(XX).LT.1.0D-10) SGNX=0.0D0
            XE=-AKB*YY
            YE=-AKB*DABS(XX)
            ZETA=CMPLX(XE,YE)
            CALL EZE1Z(XE,YE,EC,ES)
            RFL2=0.5D0*DLOG(XX**2+YY**2)
            RFT2=DATAN2(YY,XX)
            ZFC2= EC-PI*EXP(ZETA)*ZI
            ZFS2=(ES-PI*EXP(ZETA))*SGNX

            DX=XQ(J+1)-XQ(J)
            DY=YQ(J+1)-YQ(J)
            D =DSQRT(DX*DX+DY*DY)
            CDEL=DX/D
            SDEL=DY/D
            SUB =SDEL*(RFL2-RFL1)+CDEL*(RFT2-RFT1)
            ZSUB=SDEL*(ZFC2-ZFC1)-CDEL*(ZFS2-ZFS1)
            ZS(J)=ZS(J)+2.0D0/AKB*(SUB+ZSUB)
            ZD(J)=ZD(J)-2.0D0*(ZFS2-ZFS1)
            RFL1=RFL2
            RFT1=RFT2
            ZFC1=ZFC2
            ZFS1=ZFS2
      200 CONTINUE
            RETURN
            END Subroutine
!! ---------------------------------------------------------------------
!! !! ---------------------------------------------------------------------
SUBROUTINE ZSWEEP(NDIM,N,ZA,ZB,NEQ,EPS)

     IMPLICIT DOUBLE PRECISION (A-H,O-Y)
     IMPLICIT COMPLEX*16 (Z)

     Integer :: I, J, K, N, NDIM, NEQ, IP

     DIMENSION ZA(NDIM,NDIM), ZB(NDIM,NEQ)

!       Implicit None
!       Complex(16), Allocatable :: ZA(:,:), ZB(:,:)
!       Double Precision, Intent(in) :: EPS
! !! ---------------------------------------------------------------------
!       ! Integer, parameter :: MX=105, NP=100, NQ=101, NEQ=4
!       ! Double Precision, Parameter :: SML = 1.0D-14
!
!       Allocate(ZA(NDIM,NDIM), ZB(NDIM,NEQ))
!
!       ! Allocate(ZSA(MX,NP), ZSB(MX,NEQ), ZAA(NP,NP), ZBB(NP,NEQ))
!       ! Allocate(ZS(NP), SS(NP), DD(NP))

      ! DIMENSION ZA(NDIM,NDIM), ZB(NDIM,NEQ)

      DO 5 K=1,N
      P=0.0D0
      DO 1 I=K,N
      IF(P.GE.CDABS(ZA(I,K))) GOTO 1
      P=CDABS(ZA(I,K))
      IP=I
    1 CONTINUE
      IF(P.LE.EPS) GOTO 6
      IF(IP.EQ.K)  GOTO 7
        DO 2 J=K,N
        ZW=ZA(K,J)
        ZA(K,J)=ZA(IP,J)
    2   ZA(IP,J)=ZW
          DO 20 J=1,NEQ
          ZW=ZB(K,J)
          ZB(K,J)=ZB(IP,J)
   20     ZB(IP,J)=ZW
    7 CONTINUE
      IF(K.EQ.N) GOTO 70
      DO 3 J=K+1,N
    3 ZA(K,J)=ZA(K,J)/ZA(K,K)
   70   DO 30 J=1,NEQ
   30   ZB(K,J)=ZB(K,J)/ZA(K,K)
      DO 5 I=1,N
      IF(I.EQ.K) GOTO 5
      IF(K.EQ.N) GOTO 40
        DO 4 J=K+1,N
    4   ZA(I,J)=ZA(I,J)-ZA(I,K)*ZA(K,J)
   40   CONTINUE
        DO 45 J=1,NEQ
   45   ZB(I,J)=ZB(I,J)-ZA(I,K)*ZB(K,J)
    5 CONTINUE
      ZA(1,1)=(1.0D0,0.0D0)
      RETURN
    6 ZA(1,1)=CMPLX(DABS(P),0.0D0)
      RETURN
      END Subroutine
!!---------------------------------------------------------------------
!! --------------------------------------------------------------------
SUBROUTINE EZE1Z(XX,YY,EC,ES)
      ! Implicit None
      IMPLICIT DOUBLE PRECISION (A-H,O-Y)
      IMPLICIT COMPLEX*16 (Z)
      DOUBLE PRECISION  NEW

      Integer :: J, N

      DATA GAMMA /0.5772156649015D0/

      X =XX
      Y =DABS(YY)
      R =DSQRT(X*X+Y*Y)
      C =DATAN2(Y,X)

      IF(R.GT.25.0D0)  GO TO 30
      IF(X.GT.0.0D0.AND.R.GT.8.0D0)  GO TO 20
      IF(X.LE.0.0D0.AND.Y.GT.10.0D0) GO TO 20

      ER=-GAMMA-DLOG(R)+R*DCOS(C)
      EI=-C+R*DSIN(C)
      SB=-R
        DO 100 N=2,100
	    FN=DFLOAT(N)
	    CN=C*FN
	    SB=-SB*R*(FN-1.0D0)/FN/FN
	    ER=ER-SB*DCOS(CN)
	    EI=EI-SB*DSIN(CN)
	    IF(N.EQ.100)  GO TO 1
	    IF(EI.EQ.0.0D0)  GO TO 10
        IF(DABS(SB/EI).LE.1.0D-8) GO TO 10
	      GO TO 100
   10   IF(DABS(SB/ER).LE.1.0D-8) GO TO 1
  100   CONTINUE
    1 CC=DEXP(X)*DCOS(Y)
      SS=DEXP(X)*DSIN(Y)
      EC=CC*ER-SS*EI
      ES=CC*EI+SS*ER
      IF(YY.LT.0.0D0) ES=-ES
      RETURN

   20 Z =CMPLX(X,Y)
      Z1=(1.0D0,0.0D0)
      ZSUB=(10.0D0,0.0D0)
      ZS  =Z+ZSUB/(Z1+ZSUB/Z)
        DO 200 J=1,9
	      ZSUB=CMPLX(DFLOAT(10-J),0.0D0)
        ZS  =Z+ZSUB/(Z1+ZSUB/ZS)
  200   CONTINUE
      ZSUB=Z1/ZS
      EC=DREAL(ZSUB)
      ES=DIMAG(ZSUB)
      IF(YY.LT.0.0D0) ES=-ES
      RETURN

   30 OLD=-1.0D0/R
      EXC=OLD*DCOS(C)
      EXS=OLD*DSIN(C)
        DO 300 N=2,100
	      NEW=-OLD/R*DFLOAT(N-1)
	      IF(EXS.EQ.0.0D0) GO TO 31
	      IF(DABS(NEW/EXS).LE.1.0D-8) GO TO 31
	      GO TO 32
   31   IF(EXC.EQ.0.0D0) GO TO 32
        IF(DABS(NEW/EXC).LE.1.0D-8) GO TO 33
   32   IF(DABS(OLD).LT.DABS(NEW))  GO TO 33
        OLD=NEW
	      EXC=EXC+OLD*DCOS(C*DFLOAT(N))
        EXS=EXS+OLD*DSIN(C*DFLOAT(N))
  300   CONTINUE
   33 EC=-EXC
      ES=EXS
      IF(DABS(PI-DABS(C)).LT.1.0D-10) ES=-PI*DEXP(X)
      IF(YY.LT.0.0D0) ES=-ES
      RETURN
      END Subroutine
!! --------------------------------------------------------------------
End Module
