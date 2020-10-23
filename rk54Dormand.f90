MODULE rk54DormandMod

USE GlobalCommonMod
USE DyDtMod

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk54DormandEachStep
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


CONTAINS

SUBROUTINE rk54DormandEachStep(y0,yn,h,hnew,rerun,test)
 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0     ! y(t)
 REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn ! y(t+h)
 REAL(KIND=8), INTENT(IN) :: h ! initial step size
 REAL(KIND=8), INTENT(OUT) :: hnew       ! new step size
 LOGICAL, INTENT(OUT) :: test, rerun ! test=stop the program, rerun=re-run this step, reject
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: tolh ! the tolerance, determined by the problem
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ymax ! the max value among y0 and yn
 REAL(KIND=8) :: err
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6
 INTEGER :: i
  test = .False.
  rerun = .False. 
  ! use y0 to get dy0
  CALL  dev(y0,dy0,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL  dev(y1,dy1,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL  dev(y2,dy2,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y3=y0+h*(a3_0*dy0+a3_1*dy1+a3_2*dy2)
  ! use y3 to get dy3
  CALL  dev(y3,dy3,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y4=y0+h*(a4_0*dy0+a4_1*dy1+a4_2*dy2+a4_3*dy3)
  ! use y4 to get dy4
  CALL  dev(y4,dy4,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y5=y0+h*(a5_0*dy0+a5_1*dy1+a5_2*dy2+a5_3*dy3+a5_4*dy4)
  ! use y5 to get dy5
  CALL  dev(y5,dy5,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y6=y0+h*(a6_0*dy0+a6_1*dy1+a6_2*dy2+a6_3*dy3+a6_4*dy4+a6_5*dy5)
  ! use y6 to get dy6
  CALL  dev(y6,dy6,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6)
  ! Find the max value of y among this step
  DO i = 1, SIZE(y0)
    ymax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))
  END DO
  tolh = rtol*ymax + atol(1:SIZE(y0)) ! atol might be a longer array

  ! using the error to estimate the next step
  err = MAXVAL(ABS(yerr/tolh))
  IF (err.GT.1.D0) THEN
    rerun = .True.
    hnew = MAX(0.9D0*err**(-1.6666666666667D-01), 0.1D0)*h ! no less than factor of 0.1
    ! PRINT *, 'Decrease time step by', 0.9D0*err**(-1.6666666666667D-01),MAX(0.9D0*err**(-1.6666666666667D-01), 0.1D0)
  ELSE
    rerun = .False.
    hnew = MIN(5.D0, 0.9D0*err**(-1.6666666666667D-01))*h ! no more than factor of 5
    ! PRINT *, 'Increase time step by', 0.9D0*err**(-1.6666666666667D-01),MIN(5.D0,0.9D0*err**(-1.6666666666667D-01))
  END IF

  ! adjust the step
  hnew = MAX(MIN(hnew,MaxStepSize),MinStepSize)
  IF (rerun) THEN
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
      test = .True. ! stop the program
    ELSE
      hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half
    END IF
    RETURN
  END IF

  ! check if any value have went crazy (Nan or Inf)
  DO i = 1, SIZE(y0)
    CALL checkNanInf(yn(i), rerun)
    IF (rerun) THEN
      IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
        test = .True. ! stop the program
      ELSE
        hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half
      END IF
      RETURN
    END IF
  END DO

  RETURN
END SUBROUTINE

END MODULE
