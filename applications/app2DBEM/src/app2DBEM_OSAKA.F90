!! -------------------------------------------------------------------------- !!
Program app2DBEM_OSAKA
!! -------------------------------------------------------------------------- !!
    Use pkg2DBEM_OSAKA
!! -------------------------------------------------------------------------- !!
      Implicit None
      Integer :: IDIF, IRAD, IRAO
!! ----------Initial Variables Argument--------------------------------------- !!
      Integer, Parameter :: IPRINT=1, NPRINT=0, NB = 80, NTPlus = 3
      Double Precision, Parameter :: H0 = 1.0d0, SIG1 = 0.8d0, SIG2 = 0.8d0, OGD = 0.05d0, KZZB = 0.35d0, SML = 1.0D-14
!! ----------Input / Output Variables Argument related to OFFSET_M-------------!!
      Integer, Parameter :: MX=105, NP=100, NEQ = 4, NQ=101, NT = 83  !--Input
      Double Precision, Allocatable :: XP(:), YP(:), XQ(:), YQ(:)     !--Output
      Double Precision, Allocatable :: VN(:,:)                        !--Output

       !! -----Input / Output Variables Argument related to SDSUB----------!!
      Double Precision, Allocatable :: SS(:), DD(:)                   !--Output
       !! -----Input / Output Variables Argument related to SDSUB----------!!
      Complex(8), Allocatable :: ZS(:), ZD(:)                   !--Output


!! ----------Input Variables Argument related to SOLVE_M, Output text files----!!
      Double Precision :: dk
      Double Precision :: akb
      ! Variables related to SDSUB and SDCAL

      Integer :: iK
      Integer, Parameter :: nK  = 300
      Double Precision, Parameter :: kStart = 0.01D0, kEnd = 5.0D0
!! -------------------------------------------------------------------------------

!! -------------------------------------------------------------------------------
!! -------------------------------------------------------------------------------
!! -------------------------------------------------------------------------------
!! -------------------------------------------------------------------------------
!! -------------------------------------------------------------------------------
!! -------------------------------------------------------------------------------


!! ----------OFFSET Module Calling (Input/output variables corrected!)----------!!
     Call OFFSET(NB, NT, NP, NQ, MX, H0, SIG1, SIG2, OGD, KZZB, IPRINT, XP, YP, XQ, YQ, VN)
!! -------------------------------------------------------------------------------


!! ----------Output Text file writing command-----------------------------------!!

       open(newunit=IRAO,file="MotionRAO.dat",status="replace")
       write(IRAO,'(a)') "# kB/2 |X2|/A |X3|/A |X4|" // &
         "(kA) Phs(X2) Phs(X3) Phs(X4)"

       open(newunit=IDIF,file="WaveExtForce.dat",status="replace")
       write(IRAO,'(a)') "# kB/2 |F2|/A |F3|/A |F4|" // &
         "(kA) Phs(F2) Phs(F3) Phs(F4)"

       open(newunit=IRAD,file="RadiationForce.dat",status="replace")
       write(IRAD,'(a)') "# nFreq "
       write(IRAD,'(a)') "# kB/2 "
       write(IRAD,'(a)') "# AddedMass(3, 3) "
       write(IRAD,'(a)') "# Damping(3, 3) "
       write(IRAD,'(i5)') nK

!! -------------------------------------------------------------------------------
       dk = (kEnd - kStart) / ( nK - 1.D0 )
       AKB = kStart

       CALL SOLVE(NB, NT, MX, NP, NQ, NEQ, AKB, SML, XP, YP, XQ, YQ, VN)

!       do iK = 1, nK
!           CALL SOLVE (NB,NT,AKB)
!!           CALL FORCE (NB,AKB,IPRINT)
!!           CALL MOTION(AKB,IPRINT)
!           AKB = AKB + dk
!       end do
! -------------------------------------------------------------------------------
Contains
!! -------------------------------------------------------------------------- !!

End Program app2DBEM_OSAKA
