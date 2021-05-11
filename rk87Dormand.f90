MODULE rk87DormandMod

USE GlobalCommonMod
USE DyDtMod

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk87DormandEachStep
! Prince Dormand
! High order embedded Runge-Kutta formulae
! by P.J.Prince and J.R.Dormand, Journal of Computational and Applied Mathematics, vol. 7, 1981, pages 67-75.
! downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK8/RKcoeff8d_1.pdf
! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=5.55555555555555555555555555555555555555555555555555555555556D-2,&
 k2=8.33333333333333333333333333333333333333333333333333333333333D-2,&
 k3=0.125D0,&
 k4=0.3125D0,&
 k5=0.375D0,&
 k6=0.1475D0,&
 k7=0.465D0,&
 k8=0.564865451382259575398358501426168258738567008726413321510786D0,&
 k9=0.65D0,&
 k10=0.924656277640504446745013574318369542649203446702739817769306D0,&
 k11=1.0D0,&
 k12=1.0D0
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=5.55555555555555555555555555555555555555555555555555555555556D-2,&
 a2_0=2.08333333333333333333333333333333333333333333333333333333333D-2,&
 a2_1=6.25D-2,&
 a3_0=3.125D-2,&
 a3_1=0.0D0,&
 a3_2=9.375D-2,&
 a4_0=0.3125D0,&
 a4_1=0.0D0,&
 a4_2=-1.171875D0,&
 a4_3=1.171875D0,&
 a5_0=3.75D-2,&
 a5_1=0.0D0,&
 a5_2=0.0D0,&
 a5_3=0.1875D0,&
 a5_4=0.15D0,&
 a6_0=4.79101371111111111111111111111111111111111111111111111111111D-2,&
 a6_1=0.0D0,&
 a6_2=0.0D0,&
 a6_3=0.112248712777777777777777777777777777777777777777777777777778D0,&
 a6_4=-2.55056737777777777777777777777777777777777777777777777777778D-2,&
 a6_5=1.28468238888888888888888888888888888888888888888888888888889D-2,&
 a7_0=1.69179897872922811814311071360382360655149287954344106876518D-2,&
 a7_1=0.0D0,&
 a7_2=0.0D0,&
 a7_3=0.38784827848604316952654574415937335337072755585260271373015D0,&
 a7_4=3.59773698515003278967008896347723680081587394587496829956692D-2,&
 a7_5=0.196970214215666060156715256072149888128169802168423907433282D0,&
 a7_6=-0.172713852340501838761392997002333845572571026275210714846753D0,&
 a8_0=6.90957533591923006485645489845476785610413724468557040961859D-2,&
 a8_1=0.0D0,&
 a8_2=0.0D0,&
 a8_3=-0.634247976728854151882807874971738546542633603112074744302928D0,&
 a8_4=-0.161197575224604080366876923981817123442223786480849343405336D0,&
 a8_5=0.13865030945882525541986695013301580192766548894958069142448D0,&
 a8_6=0.940928614035756269724239684130258343471181148314502869775847D0,&
 a8_7=0.211636326481943981855372117131902104763536388608398143922536D0,&
 a9_0=0.183556996839045385489806023536880308449764251627732903356782D0,&
 a9_1=0.0D0,&
 a9_2=0.0D0,&
 a9_3=-2.46876808431559245274431575997410745777764975354773249558751D0,&
 a9_4=-0.291286887816300456388002572803951980054338029408172767872218D0,&
 a9_5=-2.64730202331173756884397994659461432596305528267686685125467D-2,&
 a9_6=2.84783876419280044916451825421677377023158218501195315894552D0,&
 a9_7=0.281387331469849792539403641826711782098070545536014017343881D0,&
 a9_8=0.123744899863314657627030212663639720312201353606973852326093D0,&
 a10_0=-1.21542481739588805916051052502966299488002498804411226149293D0,&
 a10_1=0.0D0,&
 a10_2=0.0D0,&
 a10_3=16.6726086659457724322804132885641077485890746007875898190766D0,&
 a10_4=0.915741828416817960595718650450742633159373633429811455667505D0,&
 a10_5=-6.0566058043574709475545055430916340040819670833244136919523D0,&
 a10_6=-16.0035735941561781118417064100788230306807930406390712252868D0,&
 a10_7=14.8493030862976625575453918980266320827229989330274563696348D0,&
 a10_8=-13.3715757352898493182930413961815957908919528946974021431482D0,&
 a10_9=5.13418264817963793317325361165860289871249428616288149527069D0,&
 a11_0=0.258860916438264283815730932231757766729630776630106316325781D0,&
 a11_1=0.0D0,&
 a11_2=0.0D0,&
 a11_3=-4.77448578548920511231011750970604274682939185374656433543205D0,&
 a11_4=-0.435093013777032509440700411810317781932355166161797411630089D0,&
 a11_5=-3.04948333207224150956051286631203161398285491122073524509586D0,&
 a11_6=5.57792003993609911742367663446494185862358894453149867930583D0,&
 a11_7=6.15583158986104009733868912668895448119775493746293746661526D0,&
 a11_8=-5.06210458673693837007740643391039164499022071214167388066848D0,&
 a11_9=2.19392617318067906127491429046580601978826270738903375925111D0,&
 a11_10=0.134627998659334941535726237887323661395585277257194651328493D0,&
 a12_0=0.822427599626507477963168204772666590957230361776585063013017D0,&
 a12_1=0.0D0,&
 a12_2=0.0D0,&
 a12_3=-11.6586732572776642839765530354584147754736908263886424759673D0,&
 a12_4=-0.757622116690936195881116154088244965366375759195411863475444D0,&
 a12_5=0.713973588159581527978269282765054675314224887856630159471613D0,&
 a12_6=12.0757749868900567395661704486006796709570580097219789708483D0,&
 a12_7=-2.12765911392040265639082085896939863542792797327511975033641D0,&
 a12_8=1.99016620704895541832807169834431415217617301735979498279352D0,&
 a12_9=-0.234286471544040292660294691856801531451242581700439470025503D0,&
 a12_10=0.175898577707942265073105105890144818314550863844624383678201D0,&
 a12_11=0.0D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=4.17474911415302462220859284685071151341936028074477848584613D-2,&
 b1=0.0D0,&
 b2=0.0D0,&
 b3=0.0D0,&
 b4=0.0D0,&
 b5=-5.5452328611239308961521894654716718893593281565192263036612D-2,&
 b6=0.239312807201180097046747354248756969660305359311065907177229D0,&
 b7=0.703510669403443023058046410889702151366379972942232998891228D0,&
 b8=-0.759759613814460929884487677085058407655423019215708812172737D0,&
 b9=0.660563030922286341461378594837820639940419712534144275764814D0,&
 b10=0.158187482510123335529614838600685443972856732403464267821114D0,&
 b11=-0.238109538752862804471863555305697193525139079217454159303497D0,&
 b12=0.25D0,&
 b_0=2.95532136763534969819648831120322465773327792500919807673162D-2,&
 b_1=0.0D0,&
 b_2=0.0D0,&
 b_3=0.0D0,&
 b_4=0.0D0,&
 b_5=-0.828606276487797039766805612688719184735403754401812296196554D0,&
 b_6=0.311240900051118327929913751626857051289736321782685558269889D0,&
 b_7=2.46734519059988698196468570406876145856225957547583462841817D0,&
 b_8=-2.54694165184190873912738007541570896178090516311468100432249D0,&
 b_9=1.44354858367677524030187495069010426851106801018924081183454D0,&
 b_10=7.94155958811272872713019541622286771314677863741958767846859D-2,&
 b_11=4.44444444444444444444444444444444444444444444444444444444444D-2,&
 b_12=0.0D0


CONTAINS

SUBROUTINE rk87DormandEachStep(y0,yn,h,hnew,rerun,test)
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
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: y1,y2,y3,y4,y5,y6,y7,y8,y9, &
                                      y10,y11,y12
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8, &
                                      dy9,dy10,dy11,dy12
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

  y8=y0+h*(a8_0*dy0+a8_1*dy1+a8_2*dy2+a8_3*dy3+a8_4*dy4+a8_5*dy5 + &
        &  a8_6*dy6+a8_7*dy7)
  ! use y8 to get dy8
  CALL  dev(y8,dy8,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y9=y0+h*(a9_0*dy0+a9_1*dy1+a9_2*dy2+a9_3*dy3+a9_4*dy4+a9_5*dy5 + &
        &  a9_6*dy6+a9_7*dy7+a9_8*dy8)
  ! use y9 to get dy9
  CALL  dev(y9,dy9,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y10=y0+h*(a10_0*dy0+a10_1*dy1+a10_2*dy2+a10_3*dy3+a10_4*dy4+a10_5*dy5 + &
         &  a10_6*dy6+a10_7*dy7+a10_8*dy8+a10_9*dy9)
  ! use y10 to get dy10
  CALL  dev(y10,dy10,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y11=y0+h*(a11_0*dy0+a11_1*dy1+a11_2*dy2+a11_3*dy3+a11_4*dy4+a11_5*dy5 + &
         &  a11_6*dy6+a11_7*dy7+a11_8*dy8+a11_9*dy9+a11_10*dy10)
  ! use y11 to get dy11
  CALL  dev(y11,dy11,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  y12=y0+h*(a12_0*dy0+a12_1*dy1+a12_2*dy2+a12_3*dy3+a12_4*dy4+a12_5*dy5 + &
         &  a12_6*dy6+a12_7*dy7+a12_8*dy8+a12_9*dy9+a12_10*dy10+a12_11*dy11)
  ! use y12 to get dy12
  CALL  dev(y12,dy12,test)
  IF (test) THEN ! bad dy 
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) RETURN ! stop the program 
    hnew = MAX(MIN(h/2.D0,MaxStepSize),MinStepSize) ! reduce step 
    rerun = .True.
    test = .False.
    RETURN
  END IF

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6+b7*dy7 + &
        &  b8*dy8+b9*dy9+b10*dy10+b11*dy11+b12*dy12)
  ynp=y0+h*(b_0*dy0+b_1*dy1+b_2*dy2+b_3*dy3+b_4*dy4+b_5*dy5+b_6*dy6+b_7*dy7 + &
         &  b_8*dy8+b_9*dy9+b_10*dy10+b_11*dy11+b_12*dy12)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1+(b2-b_2)*dy2+(b3-b_3)*dy3+(b4-b_4)*dy4 + &
       &  (b5-b_5)*dy5+(b6-b_6)*dy6+(b7-b_7)*dy7+(b8-b_8)*dy8+(b9-b_9)*dy9 + &
       &  (b10-b_10)*dy10+(b11-b_11)*dy11+(b12-b_12)*dy12)
  ! Find the max value of y among this step
  DO i = 1, SIZE(y0)
    ymax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))
  END DO
  tolh = rtol*ymax + atol(1:SIZE(y0)) ! atol might be a longer array

  ! using the error to estimate the next step
  err = MAXVAL(ABS(yerr/tolh))
  IF (err.GT.1.D0) THEN
    rerun = .True.
    hnew = MAX(0.8D0*err**(-1.D0/8.D0), 0.1D0)*h ! no less than factor of 0.1
    ! PRINT *, 'Decrease time step by', 0.8D0*err**(-1.D0/8.D0),MAX(0.8D0*err**(-1.D0/8.D0), 0.1D0)
  ELSE
    rerun = .False.
    hnew = MIN(5.D0, 0.8D0*err**(-1.D0/8.D0))*h ! no more than factor of 5
    ! PRINT *, 'Increase time step by', 0.8D0*err**(-1.D0/8.D0),MIN(5.D0,0.8D0*err**(-1.D0/8.D0))
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
    ! if any value is Nan or Inf, reRun this step
    IF (.NOT.IEEE_IS_NORMAL(yn(i))) THEN
      rerun = .True.
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
