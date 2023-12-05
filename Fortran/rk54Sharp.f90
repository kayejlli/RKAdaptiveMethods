MODULE rk54SharpMod

USE CommonMod
USE DyDtMod ! Schrodinger equation

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk54SharpEachStep

! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! Sharp and Smart
! Explicit Runge-Kutta Pairs with One More Derivative Evaluation than the Minimum, by P.W.Sharp and E.Smart,
! Siam Journal of Scientific Computing, Vol. 14, No. 2, pages. 338-348, March 1993.
! downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5q_2.pdf
! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.153577661431064572425828970331588132635253054101221640488656D0,&
 k2=0.230366492146596858638743455497382198952879581151832460732984D0,&
 k3=0.454431960049937578027465667915106117353308364544319600499376D0,&
 k4=0.718701700154559505409582689335394126738794435857805255023184D0,&
 k5=0.791443850267379679144385026737967914438502673796791443850267D0,&
 k6=1.0D0
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.153577661431064572425828970331588132635253054101221640488656D0,&
 a2_0=5.75916230366492146596858638743455497382198952879581151832461D-2,&
 a2_1=0.172774869109947643979057591623036649214659685863874345549738D0,&
 a3_0=0.218063693370253761406524706148364170760293697007780282918016D0,&
 a3_1=-0.635546527442726732584466997919464222481400785293824161025013D0,&
 a3_2=0.871914794122410549205407959686206169074415452830363478606373D0,&
 a4_0=0.207773367227052300096100839738861023302168651384515360452001D0,&
 a4_1=1.60573966518619747181414417492312948411342671677485480013666D-2,&
 a4_2=-0.137972353445870120144308475535691314309989509400978076059254D0,&
 a4_3=0.63284328972151535073964888338299312290548102670651942262907D0,&
 a5_0=-0.241781872872781056363648186242282314413782723071257619065702D0,&
 a5_1=0.601452512951897753046819443854772165840736915878165234760847D0,&
 a5_2=0.301839605320121249768538737485444795772055391318602983936203D0,&
 a5_3=-0.219112410844540646656370984072649112109552926041401535130127D0,&
 a5_4=0.349046015712682379349046015712682379349046015712682379349046D0,&
 a6_0=0.330780778802901409063313671222252973494175389716377752870426D0,&
 a6_1=-0.665181949835911939574530779632670997043302696924948240941869D0,&
 a6_2=0.597480546417610316451751888728111752198014367470977005982738D0,&
 a6_3=0.392928364282947066694335599209494837370266188957227930691965D0,&
 a6_4=-0.188444044293808824653492939697305924185999056873365171156569D0,&
 a6_5=0.532436304626261972018622560170117358166845807653730722553309D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=7.47878314168636749281910572233152878314168636749281910572233D-2,&
 b1=0.0D0,&
 b2=0.305777916293097858742966882065510845480765767095516938732436D0,&
 b3=0.185788303106221145772277785260134517271863090401110641764411D0,&
 b4=0.179016875444510802735976924753265260654432606206920778016392D0,&
 b5=0.18296240707263985115392068403110742209485500595485678376287D0,&
 b6=7.16666666666666666666666666666666666666666666666666666666667D-2,&
 b_0=7.89527012746155394400067395205399578568069323395132117920419D-2,&
 b_1=0.0D0,&
 b_2=0.28732576214199055665014282034091586607217929570050753670103D0,&
 b_3=0.217113378313197876664439154320550916910672426124205010713208D0,&
 b_4=0.125860508812105266888506018829812295547576418639351664185772D0,&
 b_5=0.220838834868425106861464537504898288840728453032288838006123D0,&
 b_6=6.99088145896656534954407294832826747720364741641337386018237D-2
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !

CONTAINS

SUBROUTINE rk54SharpEachStep(t,y0,yn,h,hNew,PleaseRerun,PleaseTerminate)
 REAL(KIND=8), INTENT(IN) :: t ! current time
 COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: y0          ! y(t)
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn  ! y(t+h)
 REAL(KIND=8), INTENT(IN) :: h           ! initial step size
 REAL(KIND=8), INTENT(OUT) :: hNew       ! new step size
 ! PleaseTerminate=1 -> stop the program =0 -> seems good
 ! PleaseRerun=1 -> re-run this step     =0 -> seems good
 LOGICAL, INTENT(OUT) :: PleaseTerminate, PleaseRerun
 ! -------------------------------------------------------------------!
 ! the following are intenal variables & parameters
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded method (may not in use)
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)) :: yMax ! the abs of max value of y
 REAL(KIND=8),  DIMENSION(SIZE(y0)) :: tolh   ! the expected error
 REAL(KIND=8) :: err, errMax
 ! ------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------ !
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)) ::y1,y2,y3,y4,y5,y6
 COMPLEX(KIND=8), DIMENSION(SIZE(y0)) ::dy0,dy1,dy2,dy3,dy4,dy5,dy6
 INTEGER :: i
  ! To initialise the logical variables
  PleaseTerminate = .False.
  PleaseRerun = .False.
  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! use y0 to get dy0
  CALL dydt(t,y0,dy0,PleaseRerun)

  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL dydt(t,y1,dy1,PleaseRerun)

  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL dydt(t,y2,dy2,PleaseRerun)

  y3=y0+h*(a3_0*dy0+a3_1*dy1+a3_2*dy2)
  ! use y3 to get dy3
  CALL dydt(t,y3,dy3,PleaseRerun)

  y4=y0+h*(a4_0*dy0+a4_1*dy1+a4_2*dy2+a4_3*dy3)
  ! use y4 to get dy4
  CALL dydt(t,y4,dy4,PleaseRerun)

  y5=y0+h*(a5_0*dy0+a5_1*dy1+a5_2*dy2+a5_3*dy3+a5_4*dy4)
  ! use y5 to get dy5
  CALL dydt(t,y5,dy5,PleaseRerun)

  y6=y0+h*(a6_0*dy0+a6_1*dy1+a6_2*dy2+a6_3*dy3+a6_4*dy4+a6_5*dy5)
  ! use y6 to get dy6
  CALL dydt(t,y6,dy6,PleaseRerun)

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
    PleaseRerun = .False. ! the error is fine, keep this step and move on
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
  IF (PleaseRerun) THEN
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min
      PleaseTerminate = .True. ! stop the program as the time step cannot be reduced any further
    END IF
    RETURN
  END IF

  ! ------------------------------------------------------------------------ !
  ! ------------------------------------------------------------------------ !
  ! check if any value have went crazy (Nan or Inf)
  DO i = 1, SIZE(y0)
    ! if any value is Nan or Inf, reRun this step
    IF (.NOT.IEEE_IS_NORMAL(ABS(yn(i)))) THEN
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
END SUBROUTINE rk54SharpEachStep


END MODULE rk54SharpMod
