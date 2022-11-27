!! ---------------------------------------------------------------------
Module OFFSET_M
!! ---------------------------------------------------------------------
    Use glb2DBEM
!! ---------------------------------------------------------------------
      Implicit None
!! ---------------------------------------------------------------------
Contains
!! ---------------------------------------------------------------------
!! ---------------------------------------------------------------------
SUBROUTINE OFFSET(NB, NT, NP, NQ, MX, H0, SIG1, SIG2, OGD, KZZB, IPRINT, XP, YP, XQ, YQ, VN)
      Implicit None
      Integer, Intent(in)          :: NB, IPRINT, NQ, NP, NT, MX
      Double Precision, Intent(in) :: H0, SIG1, SIG2, OGD, KZZB
!! ---------------------------------------------------------------------

      Integer :: IOFF, IAD, I, J, II
      Double Precision :: DTH, Sigma, RSUB, AMD, A1, A3, AMB, TH, D, DS, DX, DY, C22, CMAS, KZZ, OG, SUM, S1, S2, S3, OBM, GM

       Double precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:)
       Double Precision, Allocatable :: VN(:,:)
       Allocate(XP(MX), YP(MX), XQ(NQ), YQ(NQ))
       Allocate(VN(3,NP))

      open(newunit = IOFF, file = "Offset.dat", status = "replace")

      IAD=NT-NB
      DTH=PI/DFLOAT(NB)

      SIGMA=SIG1
      RSUB=(H0+1.0D0)**2+8.0D0*H0*(1.0D0-4.0D0*SIGMA/PI)
      AMD =0.25D0*(3.0D0*(H0+1.0D0)-DSQRT(RSUB))
      A1  =0.5D0*(H0-1.0D0)/AMD
      A3  =0.5D0*(H0+1.0D0)/AMD-1.0D0
      AMB =AMD/H0

!! ! -------------------------------------------------------------------
!!         PRINT *, '  SIGMA', '  IAD', '  DTH', '  PI', '  PI05', '  PI2',  '  OGD', '  KZZB'
!!         PRINT *, SIGMA, IAD, DTH, PI, PI05, PI2, OGD, KZZB
!! ! -------------------------------------------------------------------

      Do 100 J=1,NB/2+1
      TH=PI05-DTH*DFLOAT(J-1)
      ! 값 출력은 잘 됨, 계산 여부와 상관 없고 Array bound 문제인데..
      XQ(J)=AMB*( (1.0D0+A1)*SIN(TH)-A3*SIN(3.0D0*TH))
      YQ(J)=AMB*( (1.0D0-A1)*COS(TH)+A3*COS(3.0D0*TH))
  100 CONTINUE
      SIGMA=SIG2
      RSUB=(H0+1.0D0)**2+8.0D0*H0*(1.0D0-4.0D0*SIGMA/PI)
      AMD=0.25D0*(3.0D0*(H0+1.0D0)-DSQRT(RSUB))
      A1=0.5D0*(H0-1.0D0)/AMD
      A3=0.5D0*(H0+1.0D0)/AMD-1.0D0
      AMB=AMD/H0

      DO 105 J=NB/2+2,NB+1
      TH = PI05-DTH*DFLOAT(J-1)
      XQ(J)=AMB*((1.0D0+A1)*DSIN(TH)-A3*DSIN(3.0D0*TH))
      YQ(J)=AMB*((1.0D0-A1)*DCOS(TH)+A3*DCOS(3.0D0*TH))
  105 CONTINUE

      Do J = 1,NB+1
          write(IOFF, "(i5,99(1pe15.6))") J, XQ(J), YQ(J)
      End do
      close(IOFF)

       DO 110 I=1, NB
       XP(I)=(XQ(I+1)+XQ(I))/2.0D0
       YP(I)=(YQ(I+1)+YQ(I))/2.0D0
       DX=XQ(I+1)-XQ(I)
       DY=YQ(I+1)-YQ(I)
       D =DSQRT(DX*DX+DY*DY)
       VN(1,I)= DY/D
       VN(2,I)=-DX/D
       VN(3,I)=XP(I)*VN(2,I)-YP(I)*VN(1,I)
   110 CONTINUE

       IF(IAD.EQ.0) GOTO 130
       DS=(XQ(1)-XQ(NB+1))/DFLOAT(IAD+1)
       DO 120 I=1,IAD
       II=NB+I
       XP(II)=XQ(NB+1)+DS*DFLOAT(I)
       YP(II)=0.0D0
   120 CONTINUE

   130 CMAS=(SIG1+SIG2)/H0
       C22 =(XQ(1)-XQ(NB+1))/XQ(1)
       OG  =OGD/H0
       KZZ =KZZB
          SUM=0.0D0

!!   -------------------------------------------------------------------
!!            PRINT *, CMAS, C22, OG, KZZ
!!   -------------------------------------------------------------------

 	   DO 200 J=1,NB
          S1 =YQ(J+1)-YQ(J)
          S2 =XQ(J  )*(2.0D0*YQ(J  )+YQ(J+1))
          S3 =XQ(J+1)*(2.0D0*YQ(J+1)+YQ(J  ))
 	   SUM=SUM+S1*(S2+S3)
   200    CONTINUE
       OBM=SUM/6.0D0
       GM =(2.0D0/3.0D0-OBM)/CMAS+OG

       WRITE(6,600) CMAS,C22,OGD,KZZ,GM
       IF(IPRINT.EQ.0) RETURN
       WRITE(6,610)
       DO 300 J=1,NB+1
   300 WRITE(6,620) J,XQ(J),YQ(J),XP(J),YP(J)
   600 FORMAT( &
          15X,'NONDIMENSIONAL MASS------- S/(B/2)**2=',F8.5,/ &
         /15X,'HEAVE RESTORING FORCE COEFF--AW/(B/2)=',F8.5,/ &
         /15X,'CENTER OF GRAVITY----------------OG/D=',F8.5,/ &
         /15X,'GYRATIONAL RADIUS-----------KZZ/(B/2)=',F8.5,/ &
         /15X,'METACENTRIC HEIGHT-----------GM/(B/2)=',F8.5/)
   610 FORMAT(/15X,'***** CHECK OF ORDINATES *****' &
         /8X,'J',6X,'XQ',8X,'YQ',10X,'XP',8X,'YP')
   620 FORMAT(7X,I2,1X,2F10.5,2X,2F10.5)
       RETURN
      END Subroutine
!! ---------------------------------------------------------------------
End Module
!! ---------------------------------------------------------------------
