!! ---------------------------------------------------------------------
Module FORCE_M
!! ---------------------------------------------------------------------

   Use glb2DBEM

!! ---------------------------------------------------------------------
      IMPLICIT None
      ! Implicit None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
SUBROUTINE FORCE(NB,AK,IPRINT)
      Implicit None
      Integer, Intent(in)          :: NB
      Integer, Intent(in)          :: AK
      Integer, Intent(in)          :: IPRINT
!! ---------------------------------------------------------------------
    ! Integer :: IOFF, IAD, I,J, II
    ! Double Precision :: DTH, Sigma, RSUB, AMD, A1, A3, AMB, TH, D, DS, DX, DY, C22, CMAS, KZZ, OG, SUM, S1, S2, S3, OBM, GM
    ! Integer, parameter :: MX=105, NP=100, NQ=101
    ! Double precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:), VN(:,:)
    ! Allocate(XP(MX), YP(MX), XQ(NQ), YQ(NQ))
    ! Allocate(VN(3,NP))
    !   IMPLICIT DOUBLE PRECISION (A-H,O-Y)
    !   IMPLICIT COMPLEX*16 (Z)
    !   PARAMETER (MX=105,NP=100,NQ=101)
    !   DIMENSION A(3,3),B(3,3),BE(3,3),EAMP(3),EPHA(3)
    !   COMMON /PAI/ PI,PI05,PI2
    !   COMMON /FILEOUT/ IRAO,IDIF,IRAD
    !   COMMON /ELM/ XP(MX),YP(MX),XQ(NQ),YQ(NQ)
    !   COMMON /VN2/ VN(3,NP)
    !   COMMON /FAI/ ZFI(4,NP)
    !   COMMON /FCE/ ZAB(3,3),ZEXF(3)

      Z0=(0.0D0,0.0D0)
      ZI=(0.0D0,1.0D0)
      DO 10 I=1,3
      DO 11 J=1,3
   11 ZAB(I,J)=Z0
      ZEXF( I)=Z0
   10 CONTINUE

      DO 100 K=1,NB
      DX=XQ(K+1)-XQ(K)
      DY=YQ(K+1)-YQ(K)
      D =DSQRT(DX*DX+DY*DY)
      DO 110 I=1,3
      DO 120 J=1,3
  120 ZAB(I,J)=ZAB(I,J)-ZFI(J,K)*VN(I,K)*D
      ZEXF(I )=ZEXF(I )+ZFI(4,K)*VN(I,K)*D
  110 CONTINUE
  100 CONTINUE

      DO 150 I=1,3
      DO 160 J=1,3
      A (I,J)= DREAL(ZAB(I,J))
      B (I,J)=-DIMAG(ZAB(I,J))
  160 CONTINUE
      EAMP(I)=CDABS(ZEXF(I))
      EPHA(I)=DATAN2(DIMAG(ZEXF(I)),DREAL(ZEXF(I)))*180.0D0/PI
  150 CONTINUE

      write(IDIF, 650) AK, (EAMP(I), I=1,3), (EPHA(I), I=1,3)
      write(IRAD, 650) AK

      do I = 1, 3
      write(IRAD, 650) ( A(I, J), J = 1, 3)
      enddo

      do I = 1, 3
      write(IRAD, 650) ( B(I, J), J = 1, 3)
      end do

      IF(IPRINT.EQ.0) RETURN
      WRITE(6,600) NB,AK

      DO 300 I=1,3
      C1=B (I,I)
      C2=BE(I,I)
      CHK=DABS(C1-C2)/DABS(C1+C2)*200.0D0
  300 WRITE(6,610) I,I,A(I,I),B(I,I)
      WRITE(6,615)
      DO 310 I=1,3
      DO 310 J=1,3
      IF(I.EQ.J) GOTO 310
      WRITE(6,610) I,J,A(I,J),B(I,J)
  310 CONTINUE
      WRITE(6,630)
      DO 320 I=1,3
      WRITE(6,640) I,ZEXF(I),EAMP(I),EPHA(I)
  320 CONTINUE

  600 FORMAT(//5X,'++++++++ ADDED-MASS & DAMPING COEFF. ( ', &
         'NB=',I3,', K*B/2=',F8.4,' )  +++++++',//10X, &
         'I  J',8X,'ADDED-MASS',6X,'DAMPING')
  610 FORMAT(8X,'(',I2,',',I2,')',3X,E13.4,3(2X,E13.4))
  615 FORMAT(' ')
  620 FORMAT(8X,'(',I2,',',I2,')',3X,E13.4,2(2X,E13.4))
  630 FORMAT(//5X,'+++++ WAVE EXCITING FORCE +++++', &
       //17X,'PRESSURE INTEGRAL',12X,'AMP',5X,'PHASE(DEG)')
  640 FORMAT(8X,I2,2E13.4,2X,2E13.4,3X,E11.4,2X,F9.3)
  650 FORMAT(99(1pe15.6))
      RETURN
      END Subroutine
End Module