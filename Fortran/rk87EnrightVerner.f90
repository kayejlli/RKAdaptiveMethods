MODULE rk87EnrightVernerMod

USE parameters
USE routines

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk87EnrightVerner

! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! The Relative Efficiency of Alternative Defect Control Schemes for High-Order Continuous Runge-Kutta Formulas
! W. H. Enright SIAM Journal on Numerical Analysis, Vol. 30, No. 5. (Oct., 1993), pp. 1419-1445.
! downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8c_1.pdf
! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=5.56D-2,&
 k2=0.102577772963604852686308492201039861351819757365684575389948D0,&
 k3=0.153866659445407279029462738301559792027729636048526863084922D0,&
 k4=0.3846D0,&
 k5=0.4615D0,&
 k6=0.1538D0,&
 k7=0.8571D0,&
 k8=0.950522279549898543264080746155892215372083974054663210568629D0,&
 k9=0.7222D0,&
 k10=0.9375D0,&
 k11=1.0D0,&
 k12=1.0D0
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=5.56D-2,&
 a2_0=7.9536766850719148946376991261623370731431038799704329121008D-3,&
 a2_1=9.46240962785329377916707930748775242786766534857141424778472D-2,&
 a3_0=3.84666648613518197573656845753899480069324090121317157712305D-2,&
 a3_1=0.0D0,&
 a3_2=0.115399994584055459272097053726169844020797227036395147313692D0,&
 a4_0=0.384391795249995742463414694148930014496173869183551650147263D0,&
 a4_1=0.0D0,&
 a4_2=-1.44137549675338955770139574274596702168699721151180375962618D0,&
 a4_3=1.44158370150339381523798104859703700719082334232825210947892D0,&
 a5_0=4.61679927225245821615407392990543100173216372016370799947874D-2,&
 a5_1=0.0D0,&
 a5_2=0.0D0,&
 a5_3=0.230766671478580119613626743161835287145952544968066036170903D0,&
 a5_4=0.18456533579889529822483251753911040283672581783029688383431D0,&
 a6_0=5.98340656981684945502396474355291084616704020934744687941026D-2,&
 a6_1=0.0D0,&
 a6_2=0.0D0,&
 a6_3=0.111070988365806962422446377950496848082610655711386443707871D0,&
 a6_4=-3.42143109151919315850070051859858123538614336343020983284916D-2,&
 a6_5=1.71092568512164746123209797999598558095803758294411858265184D-2,&
 a7_0=-0.537950077527873022614747154629469286425624922015506921940171D0,&
 a7_1=0.0D0,&
 a7_2=0.0D0,&
 a7_3=-6.93764821309832024684858921678391929611997237568286028940071D0,&
 a7_4=-4.66245382097333412870032588547370265523745361906265846338828D0,&
 a7_5=3.99515211159952739816366225688709123778305091676102567472916D0,&
 a7_6=9.0D0,&
 a8_0=-1.63242744079865907890386589369373383078760689280470097367972D0,&
 a8_1=0.0D0,&
 a8_2=0.0D0,&
 a8_3=-10.827155649128676300106788103992512405570868288293796026025D0,&
 a8_4=-12.4127702165295564262541799808777575076545169312169710357507D0,&
 a8_5=9.72736897958029723333620569705898647991679730938513196401604D0,&
 a8_6=16.1993509145171614713050347911446705109706857016123854149929D0,&
 a8_7=-0.103844308090668356112325763483761031502406924627386132984838D0,&
 a9_0=0.437969506182387858999383959550753086479857907336724342215547D0,&
 a9_1=0.0D0,&
 a9_2=0.0D0,&
 a9_3=3.93953188040780609055484179804938065477461949555863432014746D0,&
 a9_4=2.86077034685671490951716275075256969722509129259119025387242D0,&
 a9_5=-1.77431070886751594245402944936881917177302571479706064516338D0,&
 a9_6=-4.89539051776420117810751382001944291072178515795222398299526D0,&
 a9_7=0.213024858819919749451668262421798635017462468974251987239988D0,&
 a9_8=-5.93953656351114879615135013862399910022202917115162753167701D-2,&
 a10_0=-1.4741971530084051096783684495220401762953156744718612221416D0,&
 a10_1=0.0D0,&
 a10_2=0.0D0,&
 a10_3=-10.9940045687738318564364554515656679159650335742933712527465D0,&
 a10_4=-11.3471035955505584070334810289625429867738425052504299382301D0,&
 a10_5=8.95698732805847648932199138433047040420587034365927857662264D0,&
 a10_6=15.8937788726768491199870127047315726092390614749240868516326D0,&
 a10_7=-9.87525742052314108953828866850888018409927677978676107605402D-2,&
 a10_8=4.88850458495220344526104177023129683434865807354905894219063D-3,&
 a10_9=-4.09681378225102871057731409693442940409595484338446331874127D-3,&
 a11_0=-2.63005933276361313355873621707188422161143965273535160089581D0,&
 a11_1=0.0D0,&
 a11_2=0.0D0,&
 a11_3=-9.17421805116351687719619947160513254615202693435764376311181D0,&
 a11_4=-19.1813926276355863096077509140025946445486735091023102723013D0,&
 a11_5=14.6425586936966487212372354845143649157868001572824958943201D0,&
 a11_6=17.5293194641808432064264792329037561372134510518264860784855D0,&
 a11_7=-0.371917560177255645679931480406461832599860350508258554936629D0,&
 a11_8=-0.700996153831514513111202618029958451946965525303009842001575D0,&
 a11_9=5.10160166123419814594373796567007747993680583372732066087135D-2,&
 a11_10=0.835689551081652570030668604041209869059346704560318853832715D0,&
 a12_0=0.215760322561474994411124388359746470315016189993958868463228D0,&
 a12_1=0.0D0,&
 a12_2=0.0D0,&
 a12_3=8.34514732667823366560737998962040201361659061600363822667194D0,&
 a12_4=2.18566238546511021703127210271410825362935054094744855941151D0,&
 a12_5=-1.68723648027586133071717517625660373873881561040340817393158D0,&
 a12_6=-8.71189790098447646994987320054742934406497913055284564417639D0,&
 a12_7=2.44414593429891254272241354773491646097520831036964631654327D-2,&
 a12_8=8.46378799450539085004935994965690483699123804486085782773911D-2,&
 a12_9=0.543485007267475889689554161135858132263172930458903122118464D0,&
 a12_10=0.0D0,&
 a12_11=0.0D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=4.39177036443990264767390168245612938204165409329946310979955D-2,&
 b1=0.0D0,&
 b2=0.0D0,&
 b3=0.0D0,&
 b4=0.0D0,&
 b5=0.351024625301198059003711947067197522471081955372891881652834D0,&
 b6=0.246142826354923978927496217611357561294520506174049687891934D0,&
 b7=0.900324493052912782292485743663086167980166191074574843765667D0,&
 b8=4.54941872725474687205992745995487359662269870465745627778562D0,&
 b9=4.80250151923706131251435406595354490562463095633152715752232D-3,&
 b10=-4.74105435212004600458540958859609881032020688796717875147351D0,&
 b11=-0.354576525007371775487465150590930876774301641201120097878066D0,&
 b12=0.0D0,&
 b_0=4.33210538122143048906481448990561225500794142960488184912532D-2,&
 b_1=0.0D0,&
 b_2=0.0D0,&
 b_3=0.0D0,&
 b_4=0.0D0,&
 b_5=0.338299786239324428483566710567298525845838208179552347378196D0,&
 b_6=0.248479842819169689968723882111796748697722834015831878169552D0,&
 b_7=0.223788296726381915717640915256391015765072451399129220461423D0,&
 b_8=-4.02016286250246552956552838417823843277019404399301309199296D-2,&
 b_9=0.123292316936198306185965165238573766500855475391670951932863D0,&
 b_10=0.0D0,&
 b_11=0.0D0,&
 b_12=6.30203320917360100491104657686662049681335571576969144866432D-2
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !

CONTAINS

SUBROUTINE rk87EnrightVernerEachStep(NEQ,y0,yn,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6,y7,y8,y9, &
                                      y10,y11,y12
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8, &
                                      dy9,dy10,dy11,dy12
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

  y7=y0+h*(a7_0*dy0+a7_1*dy1+a7_2*dy2+a7_3*dy3+a7_4*dy4+a7_5*dy5 + &
        &  a7_6*dy6)
  ! use y7 to get dy7
  CALL geo_eqns(NEQ,y7,C_t,C_phi,dy7)

  y8=y0+h*(a8_0*dy0+a8_1*dy1+a8_2*dy2+a8_3*dy3+a8_4*dy4+a8_5*dy5 + &
        &  a8_6*dy6+a8_7*dy7)
  ! use y8 to get dy8
  CALL geo_eqns(NEQ,y8,C_t,C_phi,dy8)

  y9=y0+h*(a9_0*dy0+a9_1*dy1+a9_2*dy2+a9_3*dy3+a9_4*dy4+a9_5*dy5 + &
        &  a9_6*dy6+a9_7*dy7+a9_8*dy8)
  ! use y9 to get dy9
  CALL geo_eqns(NEQ,y9,C_t,C_phi,dy9)

  y10=y0+h*(a10_0*dy0+a10_1*dy1+a10_2*dy2+a10_3*dy3+a10_4*dy4+a10_5*dy5 + &
         &  a10_6*dy6+a10_7*dy7+a10_8*dy8+a10_9*dy9)
  ! use y10 to get dy10
  CALL geo_eqns(NEQ,y10,C_t,C_phi,dy10)

  y11=y0+h*(a11_0*dy0+a11_1*dy1+a11_2*dy2+a11_3*dy3+a11_4*dy4+a11_5*dy5 + &
         &  a11_6*dy6+a11_7*dy7+a11_8*dy8+a11_9*dy9+a11_10*dy10)
  ! use y11 to get dy11
  CALL geo_eqns(NEQ,y11,C_t,C_phi,dy11)

  y12=y0+h*(a12_0*dy0+a12_1*dy1+a12_2*dy2+a12_3*dy3+a12_4*dy4+a12_5*dy5 + &
         &  a12_6*dy6+a12_7*dy7+a12_8*dy8+a12_9*dy9+a12_10*dy10+a12_11*dy11)
  ! use y12 to get dy12
  CALL geo_eqns(NEQ,y12,C_t,C_phi,dy12)

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6+b7*dy7 + &
        &  b8*dy8+b9*dy9+b10*dy10+b11*dy11+b12*dy12)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6+b_7*dy7 + &
         &  b_8*dy8+b_9*dy9+b_10*dy10+b_11*dy11+b_12*dy12)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6+(b7-b_7)*dy7+(b8-b_8)*dy8+(b9-b_9)*dy9 + &
       &  (b10-b_10)*dy10+(b11-b_11)*dy11+(b12-b_12)*dy12)
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
    hNew = MAX(0.8D0*err**(-1.D0/8.D0), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-1.D0/8.D0),MAX(0.8D0*err**(-1.D0/8.D0), ReduceAtMost)
  ELSE
    PleaseRerun = .False. ! the error is fine, keep this step and move on
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-1.D0/8.D0))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-1.D0/8.D0),MIN(IncreaseAtMost,0.8D0*err**(-1.D0/8.D0))
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
END SUBROUTINE rk87EnrightVernerEachStep


! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! This subroutine choose the time step based on accuracy required and suggest
!   the next step
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
SUBROUTINE rk87EnrightVerner(NEQ,X0,XN,Y0,C_t,C_phi,EPS,YN,IERR,h)
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
    CALL rk87EnrightVernerEachStep(NEQ,Y0,YN,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
END SUBROUTINE rk87EnrightVerner


END MODULE
