MODULE rk65DormandMod

USE GlobalCommonMod
USE DyDtMod

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk65DormandEachStep
! P.J. Prince and J. R. Dormand, High order embedded Runge-Kutta formulae,
! Journal of Computational and Applied Mathematics . 7 (1981), pp. 67-75.
! downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK6/RKcoeff6e_1.pdf
! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.1D0,&
 k2=0.222222222222222222222222222222222222222222222222222222222222D0,&
 k3=0.428571428571428571428571428571428571428571428571428571428571D0,&
 k4=0.6D0,&
 k5=0.8D0,&
 k6=1.0D0,&
 k7=1.0D0
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.1D0,&
 a2_0=-2.46913580246913580246913580246913580246913580246913580246914D-2,&
 a2_1=0.246913580246913580246913580246913580246913580246913580246914D0,&
 a3_0=0.448250728862973760932944606413994169096209912536443148688047D0,&
 a3_1=-0.787172011661807580174927113702623906705539358600583090379009D0,&
 a3_2=0.767492711370262390670553935860058309037900874635568513119534D0,&
 a4_0=0.589636363636363636363636363636363636363636363636363636363636D0,&
 a4_1=-0.981818181818181818181818181818181818181818181818181818181818D0,&
 a4_2=0.712573426573426573426573426573426573426573426573426573426573D0,&
 a4_3=0.279608391608391608391608391608391608391608391608391608391608D0,&
 a5_0=-0.713589225589225589225589225589225589225589225589225589225589D0,&
 a5_1=1.30909090909090909090909090909090909090909090909090909090909D0,&
 a5_2=0.120128342245989304812834224598930481283422459893048128342246D0,&
 a5_3=-0.652013468013468013468013468013468013468013468013468013468013D0,&
 a5_4=0.736383442265795206971677559912854030501089324618736383442266D0,&
 a6_0=2.34048821548821548821548821548821548821548821548821548821549D0,&
 a6_1=-3.18181818181818181818181818181818181818181818181818181818182D0,&
 a6_2=-0.763123754073980318324209726924659051355883935069455431446382D0,&
 a6_3=4.48261211722750184288645827107365568904030442491980953519415D0,&
 a6_4=-2.84586056644880174291938997821350762527233115468409586056645D0,&
 a6_5=0.96770216962524654832347140039447731755424063116370808678501D0,&
 a7_0=1.74913946007696007696007696007696007696007696007696007696008D0,&
 a7_1=-2.39042207792207792207792207792207792207792207792207792207792D0,&
 a7_2=-0.396252573783682380967448840752008172822652460661510435266091D0,&
 a7_3=3.27285832934871396409857948319486781025242563704102165640627D0,&
 a7_4=-2.06351637877373171490818549642079053843759726112667289137877D0,&
 a7_5=0.828193241053817976894899971823048746125669202592279515356438D0,&
 a7_6=0.0D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=7.06018518518518518518518518518518518518518518518518518518519D-2,&
 b1=0.0D0,&
 b2=0.305849410770225249863258913032668688777286062353935657103078D0,&
 b3=0.115103824238439623055007670392285776901161516546131930747315D0,&
 b4=0.187227668845315904139433551198257080610021786492374727668845D0,&
 b5=0.254252958579881656804733727810650887573964497041420118343195D0,&
 b6=-3.30357142857142857142857142857142857142857142857142857142857D-2,&
 b7=0.1D0,&
 b_0=7.60185185185185185185185185185185185185185185185185185185185D-2,&
 b_1=0.0D0,&
 b_2=0.274041072050121823877479986077271145144448311869126348764358D0,&
 b_3=0.19205895244356782818321279859741398202936664475126013587552D0,&
 b_4=0.1075708061002178649237472766884531590413943355119825708061D0,&
 b_5=0.290310650887573964497041420118343195266272189349112426035503D0,&
 b_6=6.0D-2,&
 b_7=0.0D0


CONTAINS

SUBROUTINE rk65DormandEachStep(y0,yn,h,hnew,rerun,test)
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
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6,y7
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7
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

  y7=y0+h*(a7_0*dy0+a7_1*dy1+a7_2*dy2+a7_3*dy3+a7_4*dy4+a7_5*dy5 + &
        &  a7_6*dy6)
  ! use y7 to get dy7
  CALL  dev(y7,dy7,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6+b7*dy7)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6+b_7*dy7)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6+(b7-b_7)*dy7)
  ! Find the max value of y among this step
  DO i = 1, SIZE(y0)
    ymax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))
  END DO
  tolh = rtol*ymax + atol(1:SIZE(y0)) ! atol might be a longer array

  ! using the error to estimate the next step
  err = MAXVAL(ABS(yerr/tolh))
  IF (err.GT.1.D0) THEN
    rerun = .True.
    hnew = MAX(0.8D0*err**(-1.6666666666667D-01), 0.1D0)*h ! no less than factor of 0.1
    ! PRINT *, 'Decrease time step by', 0.8D0*err**(-1.6666666666667D-01),MAX(0.8D0*err**(-1.6666666666667D-01), 0.1D0)
  ELSE
    rerun = .False.
    hnew = MIN(5.D0, 0.8D0*err**(-1.6666666666667D-01))*h ! no more than factor of 5
    ! PRINT *, 'Increase time step by', 0.8D0*err**(-1.6666666666667D-01),MIN(5.D0,0.8D0*err**(-1.6666666666667D-01))
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
