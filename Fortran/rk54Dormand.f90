MODULE rk54DormandMod

USE GlobalCommonMod
USE DyDtMod

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk54DormandEachStep

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

SUBROUTINE rk54DormandEachStep(t,y0,yn,h,hNew,PleaseRerun,PleaseTerminate)
 REAL(KIND=8), INTENT(IN) :: t ! current time
 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0          ! y(t)
 REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn  ! y(t+h)
 REAL(KIND=8), INTENT(IN) :: h           ! initial step size
 REAL(KIND=8), INTENT(OUT) :: hNew       ! new step size
 ! PleaseTerminate=1 -> stop the program =0 -> seems good
 ! PleaseRerun=1 -> re-run this step     =0 -> seems good
 LOGICAL, INTENT(OUT) :: PleaseTerminate, PleaseRerun
 ! -------------------------------------------------------------------!
 ! the following are intenal variables & parameters
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded method (may not in use)
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yMax ! the abs of max value of y
 REAL(KIND=8),  DIMENSION(SIZE(y0)) :: tolh   ! the expected error
 REAL(KIND=8) :: err, errMax
 ! ------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------ !
 REAL(KIND=8), DIMENSION(SIZE(y0)) ::y1,y2,y3,y4,y5,y6
 REAL(KIND=8), DIMENSION(SIZE(y0)) ::dy0,dy1,dy2,dy3,dy4,dy5,dy6
 INTEGER :: i
  ! To initialise the logical variables
  PleaseTerminate = .False.
  PleaseRerun = .False.
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! use y0 to get dy0
  CALL dev(t,y0,dy0,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL dev(t,y1,dy1,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL dev(t,y2,dy2,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y3=y0+h*(a3_0*dy0+a3_1*dy1+a3_2*dy2)
  ! use y3 to get dy3
  CALL dev(t,y3,dy3,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y4=y0+h*(a4_0*dy0+a4_1*dy1+a4_2*dy2+a4_3*dy3)
  ! use y4 to get dy4
  CALL dev(t,y4,dy4,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y5=y0+h*(a5_0*dy0+a5_1*dy1+a5_2*dy2+a5_3*dy3+a5_4*dy4)
  ! use y5 to get dy5
  CALL dev(t,y5,dy5,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  y6=y0+h*(a6_0*dy0+a6_1*dy1+a6_2*dy2+a6_3*dy3+a6_4*dy4+a6_5*dy5)
  ! use y6 to get dy6
  CALL dev(t,y6,dy6,PleaseTerminate)
  IF (PleaseTerminate) THEN ! bad dy
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program
      hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step
      PleaseRerun = .True.
      PleaseTerminate = .False.
      RETURN
  END IF

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6)
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! construct yMax for rtol
  DO i = 1, SIZE(y0)
    yMax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))
  END DO

  ! ------------------------------------------------------------------------ !
  ! now define the tol-array
  tolh = rtol*yMax + atol(1:SIZE(y0)) ! atol might be longer
  ! using the error to estimate the next step
  err = MAXVAL(ABS(yerr/tolh))
  IF (err.GT.1.D0) THEN
    PleaseRerun = .True.  ! the error is too large, PleaseRerun this step
    ! ReduceAtMost is suggested to be 0.1D0 or 0.05D0
    hNew = MAX(0.8D0*err**(-1.D0/6.D0), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-1.D0/6.D0),MAX(0.8D0*err**(-1.D0/6.D0), ReduceAtMost)
  ELSE
    ! PleaseRerun = .False. ! the error is fine, keep this step and move on
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-1.D0/6.D0))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-1.D0/6.D0),MIN(IncreaseAtMost,0.8D0*err**(-1.D0/6.D0))
  END IF
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! adjust the step (make sure it is bounded with MinStepSize & MaxStepSize)
  hNew = MAX(MIN(hNew,MaxStepSize),MinStepSize)

  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! check if any value have went crazy (Nan or Inf)
  DO i = 1, SIZE(y0)
    ! if any value is Nan or Inf, reRun this step
    IF (.NOT.IEEE_IS_NORMAL(yn(i))) THEN
      PleaseRerun = .True.
      IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
        PleaseTerminate = .True. ! stop the program
      ELSE
        ! hNew is likely not computated due to the presence of Nan or Inf value in yn & yerr
        hNew = MAX(MinStepSize,MIN(hNew, h*0.5D0)) ! just take the smaller one
      END IF
      RETURN
    END IF
  END DO
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !

 RETURN
END SUBROUTINE rk54DormandEachStep


END MODULE rk54DormandMod
