MODULE rk54DormandMod

USE parameters
USE routines

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk54Dormand

! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! A family of embedded Runge-Kutta formulae, by J. R. Dormand and P. J. Prince, 
! Journal of Computational and Applied Mathematics, Vol 6, No. 1, 1980, pages 19 to 26.
! downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5d_1.pdf
! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.2D0,&
 k2=0.3D0,&
 k3=0.8D0,&
 k4=0.888888888888888888888888888888888888888888888888888888888889D0,&
 k5=1.0D0,&
 k6=1.0D0
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.2D0,&
 a2_0=7.5D-2,&
 a2_1=0.225D0,&
 a3_0=0.977777777777777777777777777777777777777777777777777777777778D0,&
 a3_1=-3.73333333333333333333333333333333333333333333333333333333333D0,&
 a3_2=3.55555555555555555555555555555555555555555555555555555555556D0,&
 a4_0=2.95259868922420362749580856576741350403901844231062338058223D0,&
 a4_1=-11.5957933241883859167809785093735711019661636945587562871513D0,&
 a4_2=9.82289285169943606157597927145252248132906569120560890108215D0,&
 a4_3=-0.290809327846364883401920438957475994513031550068587105624143D0,&
 a5_0=2.84627525252525252525252525252525252525252525252525252525253D0,&
 a5_1=-10.7575757575757575757575757575757575757575757575757575757576D0,&
 a5_2=8.90642271774347246045359252906422717743472460453592529064227D0,&
 a5_3=0.278409090909090909090909090909090909090909090909090909090909D0,&
 a5_4=-0.27353130360205831903945111492281303602058319039451114922813D0,&
 a6_0=9.11458333333333333333333333333333333333333333333333333333333D-2,&
 a6_1=0.0D0,&
 a6_2=0.449236298292902066486972147349505840071877807726864330637916D0,&
 a6_3=0.651041666666666666666666666666666666666666666666666666666667D0,&
 a6_4=-0.322376179245283018867924528301886792452830188679245283018868D0,&
 a6_5=0.130952380952380952380952380952380952380952380952380952380952D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=9.11458333333333333333333333333333333333333333333333333333333D-2,&
 b1=0.0D0,&
 b2=0.449236298292902066486972147349505840071877807726864330637916D0,&
 b3=0.651041666666666666666666666666666666666666666666666666666667D0,&
 b4=-0.322376179245283018867924528301886792452830188679245283018868D0,&
 b5=0.130952380952380952380952380952380952380952380952380952380952D0,&
 b6=0.0D0,&
 b_0=8.99131944444444444444444444444444444444444444444444444444444D-2,&
 b_1=0.0D0,&
 b_2=0.453489068583408206049715483677747828691224917640011979634621D0,&
 b_3=0.6140625D0,&
 b_4=-0.271512382075471698113207547169811320754716981132075471698113D0,&
 b_5=8.90476190476190476190476190476190476190476190476190476190476D-2,&
 b_6=2.5D-2
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !

CONTAINS

SUBROUTINE rk54DormandEachStep(NEQ,y0,yn,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
 INTEGER, INTENT(IN) :: NEQ ! Dimension of y(:)
 REAL(KIND=8), DIMENSION(NEQ), INTENT(IN) :: y0     ! y(t)
 REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn    ! y(t+h)
 REAL(KIND=8), INTENT(IN) :: h           ! initial step size
 REAL(KIND=8), INTENT(OUT) :: hNew       ! new step size
 ! PleaseTerminate=1 -> stop the program =0 -> seems good
 ! PleaseRerun=1 -> re-run this step     =0 -> seems good
 LOGICAL, INTENT(OUT) :: PleaseTerminate, PleaseRerun
 REAL(KIND=8), INTENT(IN) :: EPS ! the epsilon for error
 REAL(KIND=8), INTENT(IN) :: C_t, C_phi ! parameters needed by geo_eqnst
 ! -------------------------------------------------------------------!
 ! the following are intenal variables & parameters
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded method
 REAL(KIND=8) :: err, errMax
 ! ------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------ !
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6
 INTEGER :: i
  ! To initialise the logical variables
  PleaseTerminate = .False.
  PleaseRerun = .False.
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! use y0 to get dy0
  CALL geo_eqns(NEQ,y0,C_t,C_phi,dy0)

  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL geo_eqns(NEQ,y1,C_t,C_phi,dy1)

  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL geo_eqns(NEQ,y2,C_t,C_phi,dy2)

  y3=y0+h*(a3_0*dy0+a3_1*dy1+a3_2*dy2)
  ! use y3 to get dy3
  CALL geo_eqns(NEQ,y3,C_t,C_phi,dy3)

  y4=y0+h*(a4_0*dy0+a4_1*dy1+a4_2*dy2+a4_3*dy3)
  ! use y4 to get dy4
  CALL geo_eqns(NEQ,y4,C_t,C_phi,dy4)

  y5=y0+h*(a5_0*dy0+a5_1*dy1+a5_2*dy2+a5_3*dy3+a5_4*dy4)
  ! use y5 to get dy5
  CALL geo_eqns(NEQ,y5,C_t,C_phi,dy5)

  y6=y0+h*(a6_0*dy0+a6_1*dy1+a6_2*dy2+a6_3*dy3+a6_4*dy4+a6_5*dy5)
  ! use y6 to get dy6
  CALL geo_eqns(NEQ,y6,C_t,C_phi,dy6)

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6)
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! Find the max value of yerr
  errMax = MAXVAL(ABS(yerr))
  ! using the error to propose the next step
  err = ABS(errMax/EPS)
  ! Increase the time step or decrease it
  IF (err.GT.1.D0) THEN
    PleaseRerun = .True.  ! the error is too large, PleaseRerun this step
    ! ReduceAtMost is suggested to be 0.1D0 or 0.05D0
    hNew = MAX(0.8D0*err**(-1.D0/6.D0), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-1.D0/6.D0),MAX(0.8D0*err**(-1.D0/6.D0), ReduceAtMost)
  ELSE
    PleaseRerun = .False. ! the error is fine, keep this step and move on
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-1.D0/6.D0))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-1.D0/6.D0),MIN(IncreaseAtMost,0.8D0*err**(-1.D0/6.D0))
  END IF

  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! adjust the step (make sure it is bounded with MinStepSize & MaxStepSize)
  hNew = MAX(MIN(hNew,MaxStepSize),MinStepSize)
  IF (PleaseRerun) THEN
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
      PleaseTerminate = .True. ! stop the program as the time step cannot be reduced any further
    END IF
    RETURN
  END IF

  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! check if any value have went crazy (Nan or Inf)
  DO i = 1, NEQ
    ! if any value is Nan or Inf, reRun this step
    IF (.NOT.IEEE_IS_NORMAL(yn(i))) THEN
      PleaseRerun = .True.
      IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
        PleaseTerminate = .True. ! stop the program
      ELSE
        ! hNew is likely not computated due to the presence of Nan or Inf value in yn & yerr
        IF (.NOT.IEEE_IS_NORMAL(hNew)) THEN
          hNew = MAX(MinStepSize,MIN(h, h*0.5D0)) ! just try to reduce it by half
        ELSE
          hNew = MAX(MinStepSize,MIN(hNew, h)) ! just take the smaller one
        END IF
      END IF
      RETURN
    END IF
  END DO
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

 RETURN
END SUBROUTINE rk54DormandEachStep


! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! This subroutine choose the time step based on accuracy required and suggest
!   the next step
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
SUBROUTINE rk54Dormand(NEQ,X0,XN,Y0,C_t,C_phi,EPS,YN,IERR,h)
  ! -------------------------------------------------------------------!
  ! Solve the system of ordinary differential equations                !
  ! -------------------------------------------------------------------!
  INTEGER, INTENT(IN) :: NEQ ! dimension of Y0 and YN
  REAL(KIND=8), INTENT(IN) :: X0    ! the current time (t)
  REAL(KIND=8), INTENT(INOUT) :: XN ! the new time == X0 + time step (t+h)
  REAL(KIND=8), INTENT(OUT) :: h    ! the time step suggested for the next step
  REAL(KIND=8), DIMENSION(NEQ), INTENT(IN) :: Y0  ! y(t)
  REAL(KIND=8), DIMENSION(NEQ), INTENT(OUT) :: YN ! y(t+h)
  REAL(KIND=8), INTENT(IN) :: EPS ! the epsilon for error
  REAL(KIND=8), INTENT(IN) :: C_t, C_phi ! parameters needed by geo_eqns
  INTEGER, INTENT(OUT) :: IERR ! the error flag 1==trouble, 0==good
  ! -------------------------------------------------------------------!
  ! the following are intenal variables & parameters
  INTEGER :: ITER ! for iterations
  REAL(KIND=8) :: hOld, hNew ! saving the time step
  LOGICAL :: PleaseRerun, PleaseTerminate
  ! -------------------------------------------------------------------!

  ! Initialise error flag
  IERR = 0
  ! Initialise the time step
  h = XN - X0

  ! Iterate until the error is smaller than EPS
  DO ITER = 1, MaxNumberOfIteration
    ! ----------------------------------------------------------- !
    ! Save the time step adopted (to compare with the new suggested time step)
    hOld = h
    ! ----------------------------------------------------------- !
    ! ----------------------------------------------------------- !
    ! Try to move the system forward by h and calculate the suggested time step -> hNew
    !   YN is updated to y(t+h)
    !   hNew is updated to the suggested time step
    !   PleaseRerun & PleaseTerminate are updated
    CALL rk54DormandEachStep(NEQ,Y0,YN,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
    ! ----------------------------------------------------------- !
    ! ----------------------------------------------------------- !
    IF (.NOT.PleaseRerun) THEN
      ! This solver is happy with the time step chosen
      ! At this point, hOld==h is the time step used for y(t) -> y(t+)
      !                hNew is the time step suggested for the next step
      EXIT
    END IF
    ! ----------------------------------------------------------- !
    ! ----------------------------------------------------------- !
    ! Now the program decide to re-run this step
    !   update the time step to the suggested one
    h = hNew
    ! ----------------------------------------------------------- !
    ! ----------------------------------------------------------- !
    ! If this is already the last iteration, raise an error
    !    If the solver terminate for this reason, you could increase MaxNumberOfIteration
    IF (ITER==MaxNumberOfIteration) THEN
      ! PRINT *, "Max iteration is reached, but the solver still gives a large error"
      IERR = 1
      EXIT
    END IF
    ! ----------------------------------------------------------- !
    ! If Nan or Inf value appears and persist even though time step has been reduced
    IF (PleaseTerminate) THEN
      ! PRINT *, "Serious problem happened, please terminate the program"
      ! PRINT *, "Most likely Nan or Inf value appears but decreasing the time step does not help"
      IERR = 2
      EXIT
    END IF
    ! ----------------------------------------------------------- !
    ! If the solver is not happy with the error obtained
    !   but the suggested time step is the same as the input time step
    IF ((Abs(hNew-hOld)/hOld).LT.1D-13) THEN
      ! PRINT *, "The solver is asking for re-run, but the suggested time step remains the same"
      IERR = 3
      EXIT
    END IF
    ! ----------------------------------------------------------- !
    ! If the solver is not happy with the error obtained
    !   but the time step we used already reaches the Min time step allowed
    IF ((Abs(hOld-MinStepSize)/MinStepSize).LT.1D-13) THEN
      ! PRINT *, "The solver is asking for re-run, but the time step used already reaches min time step allowed"
      IERR = 4
      EXIT
    END IF
    ! ----------------------------------------------------------- !
    ! ----------------------------------------------------------- !
  END DO

  ! ----------------------------------------------------------- !
  ! The values are updated even if the program wants to terminate
  ! ----------------------------------------------------------- !
  XN = X0 + hOld ! hOld is used to forward y(t) to y(t+)
  h = hNew       ! Save the suggested time step to h

  ! ----------------------------------------------------------- !
  ! It is recommended to terminate the program if IERR > 0
  ! ----------------------------------------------------------- !
  RETURN
END SUBROUTINE rk54Dormand


END MODULE
