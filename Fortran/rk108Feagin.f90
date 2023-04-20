MODULE rk108FeaginMod

USE parameters
USE routines

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: rk108Feagin

! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! THE COEFFICIENTS OF RK10(8)  TO 60 DIGITS
! using the notation of Fehlberg, Bettis, Horn, et alia


! k[i] are the nodes 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.1D0,&
 k2=0.5393578408029817875324851978813024368572734497010090155055D0,&
 k3=0.80903676120447268129872779682195365528591017455151352325825D0,&
 k4=0.30903676120447268129872779682195365528591017455151352325825D0,&
 k5=0.981074190219795268254879548310562080489056746118724882027805D0,&
 k6=0.833333333333333333333333333333333333333333333333333333333333D0,&
 k7=0.354017365856802376329264185948796742115824053807373968324184D0,&
 k8=0.882527661964732346425501486979669075182867844268052119663791D0,&
 k9=0.642615758240322548157075497020439535959501736363212695909875D0,&
 k10=0.357384241759677451842924502979560464040498263636787304090125D0,&
 k11=0.117472338035267653574498513020330924817132155731947880336209D0,&
 k12=0.833333333333333333333333333333333333333333333333333333333333D0,&
 k13=0.30903676120447268129872779682195365528591017455151352325825D0,&
 k14=0.5393578408029817875324851978813024368572734497010090155055D0,&
 k15=0.1D0,&
 k16=1.0D0
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=3.33333333333333333333333333333333333333333333333333333333333D-2,&
 b1=2.5D-2,&
 b2=3.33333333333333333333333333333333333333333333333333333333333D-2,&
 b3=0.0D0,&
 b4=5.0D-2,&
 b5=0.0D0,&
 b6=4.0D-2,&
 b7=0.0D0,&
 b8=0.189237478148923490158306404106012326238162346948625830327194D0,&
 b9=0.277429188517743176508360262560654340428504319718040836339472D0,&
 b10=0.277429188517743176508360262560654340428504319718040836339472D0,&
 b11=0.189237478148923490158306404106012326238162346948625830327194D0,&
 b12=-4.0D-2,&
 b13=-5.0D-2,&
 b14=-3.33333333333333333333333333333333333333333333333333333333333D-2,&
 b15=-2.5D-2,&
 b16=3.33333333333333333333333333333333333333333333333333333333333D-2
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.1D0,&
 a2_0=-0.91517656137529144052001501927534215431895138766436972056466D0,&
 a2_1=1.45453440217827322805250021715664459117622483736537873607016D0,&
 a3_0=0.202259190301118170324681949205488413821477543637878380814562D0,&
 a3_1=0.0D0,&
 a3_2=0.606777570903354510974045847616465241464432630913635142443687D0,&
 a4_0=0.184024714708643575149100693471120664216774047979591417844635D0,&
 a4_1=0.0D0,&
 a4_2=0.197966831227192369068141770510388793370637287463360401555746D0,&
 a4_3=-7.29547847313632629185146671595558023015011608914382961421311D-2,&
 a5_0=8.79007340206681337319777094132125475918886824944548534041378D-2,&
 a5_1=0.0D0,&
 a5_2=0.0D0,&
 a5_3=0.410459702520260645318174895920453426088035325902848695210406D0,&
 a5_4=0.482713753678866489204726942976896106809132737721421333413261D0,&
 a6_0=8.5970050490246030218848022594580840141113261563660022259388D-2,&
 a6_1=0.0D0,&
 a6_2=0.0D0,&
 a6_3=0.330885963040722183948884057658753173648240154838402033448632D0,&
 a6_4=0.48966295730945019284450701113589820117801547843379009721079D0,&
 a6_5=-7.31856375070850736789057580558988816340355615025188195854775D-2,&
 a7_0=0.120930449125333720660378854927668953958938996999703678812621D0,&
 a7_1=0.0D0,&
 a7_2=0.0D0,&
 a7_3=0.0D0,&
 a7_4=0.260124675758295622809007617838335174368108756484693361887839D0,&
 a7_5=3.25402621549091330158899334391231259332716675992700000776101D-2,&
 a7_6=-5.95780211817361001560122202563305121444953672762930724538856D-2,&
 a8_0=0.110854379580391483508936171010218441909425780168656559807038D0,&
 a8_1=0.0D0,&
 a8_2=0.0D0,&
 a8_3=0.0D0,&
 a8_4=0.0D0,&
 a8_5=-6.05761488255005587620924953655516875526344415354339234619466D-2,&
 a8_6=0.32176370560177839010089879904987890408140436860307712925111D0,&
 a8_7=0.51048572560806303157775901228512341674467213703175235406759D0,&
 a9_0=0.112054414752879004829715002761802363003717611158172229329393D0,&
 a9_1=0.0D0,&
 a9_2=0.0D0,&
 a9_3=0.0D0,&
 a9_4=0.0D0,&
 a9_5=-0.144942775902865915672349828340980777181668499748506838876185D0,&
 a9_6=-0.333269719096256706589705211415746871709467423992115497968724D0,&
 a9_7=0.499269229556880061353316843969978567860276816592673201240332D0,&
 a9_8=0.50950460892968610423609869004538625398664323235298960218506D0,&
 a10_0=0.113976783964185986138004186736901163890724752541486831640341D0,&
 a10_1=0.0D0,&
 a10_2=0.0D0,&
 a10_3=0.0D0,&
 a10_4=0.0D0,&
 a10_5=-7.68813364203356938586214289120895270821349023390922987406384D-2,&
 a10_6=0.239527360324390649107711455271882373019741311201004119339563D0,&
 a10_7=0.397774662368094639047830462488952104564716416343454639902613D0,&
 a10_8=1.07558956873607455550609147441477450257136782823280838547024D-2,&
 a10_9=-0.327769124164018874147061087350233395378262992392394071906457D0,&
 a11_0=7.98314528280196046351426864486400322758737630423413945356284D-2,&
 a11_1=0.0D0,&
 a11_2=0.0D0,&
 a11_3=0.0D0,&
 a11_4=0.0D0,&
 a11_5=-5.20329686800603076514949887612959068721311443881683526937298D-2,&
 a11_6=-5.76954146168548881732784355283433509066159287152968723021864D-2,&
 a11_7=0.19478191571210416497630626214738287115614292135440936473809D0,&
 a11_8=0.145384923188325069727524825977071194859203467568236523866582D0,&
 a11_9=-7.82942710351670777553986729725692447252077047239160551335016D-2,&
 a11_10=-0.114503299361098912184303164290554670970133218405658122674674D0,&
 a12_0=0.985115610164857280120041500306517278413646677314195559520529D0,&
 a12_1=0.0D0,&
 a12_2=0.0D0,&
 a12_3=0.330885963040722183948884057658753173648240154838402033448632D0,&
 a12_4=0.48966295730945019284450701113589820117801547843379009721079D0,&
 a12_5=-1.37896486574843567582112720930751902353904327148559471526397D0,&
 a12_6=-0.861164195027635666673916999665534573351026060987427093314412D0,&
 a12_7=5.78428813637537220022999785486578436006872789689499172601856D0,&
 a12_8=3.28807761985103566890460615937314805477268252903342356581925D0,&
 a12_9=-2.38633905093136384013422325215527866148401465975954104585807D0,&
 a12_10=-3.25479342483643918654589367587788726747711504674780680269911D0,&
 a12_11=-2.16343541686422982353954211300054820889678036420109999154887D0,&
 a13_0=0.89508029577163289104961313233658513814815627924156134599171D0,&
 a13_1=0.0D0,&
 a13_2=0.197966831227192369068141770510388793370637287463360401555746D0,&
 a13_3=-7.29547847313632629185146671595558023015011608914382961421311D-2,&
 a13_4=0.0D0,&
 a13_5=-0.851236239662007619739049371445966793289359722875702227166105D0,&
 a13_6=0.398320112318533301719718614174373643336480918103773904231856D0,&
 a13_7=3.63937263181035606029412920047090044132027387893977804176229D0,&
 a13_8=1.54822877039830322365301663075174564919981736348973496313065D0,&
 a13_9=-2.12221714704053716026062427460427261025318461146260124401561D0,&
 a13_10=-1.58350398545326172713384349625753212757269188934434237975291D0,&
 a13_11=-1.71561608285936264922031819751349098912615880827551992973034D0,&
 a13_12=-2.44036405750127452135415444412216875465593598370910566069132D-2,&
 a14_0=-0.91517656137529144052001501927534215431895138766436972056466D0,&
 a14_1=1.45453440217827322805250021715664459117622483736537873607016D0,&
 a14_2=0.0D0,&
 a14_3=0.0D0,&
 a14_4=-0.777333643644968233538931228575302137803351053629547286334469D0,&
 a14_5=0.0D0,&
 a14_6=-9.10895662155176069593203555807484200111889091770101799647985D-2,&
 a14_7=0.0D0,&
 a14_8=0.0D0,&
 a14_9=0.0D0,&
 a14_10=0.0D0,&
 a14_11=0.0D0,&
 a14_12=9.10895662155176069593203555807484200111889091770101799647985D-2,&
 a14_13=0.777333643644968233538931228575302137803351053629547286334469D0,&
 a15_0=0.1D0,&
 a15_1=0.0D0,&
 a15_2=-0.157178665799771163367058998273128921867183754126709419409654D0,&
 a15_3=0.0D0,&
 a15_4=0.0D0,&
 a15_5=0.0D0,&
 a15_6=0.0D0,&
 a15_7=0.0D0,&
 a15_8=0.0D0,&
 a15_9=0.0D0,&
 a15_10=0.0D0,&
 a15_11=0.0D0,&
 a15_12=0.0D0,&
 a15_13=0.0D0,&
 a15_14=0.157178665799771163367058998273128921867183754126709419409654D0,&
 a16_0=0.181781300700095283888472062582262379650443831463199521664945D0,&
 a16_1=0.675D0,&
 a16_2=0.34275815984718983994222055341385087174233873470395891993726D0,&
 a16_3=0.0D0,&
 a16_4=0.259111214548322744512977076191767379267783684543182428778156D0,&
 a16_5=-0.358278966717952089048961276721979397739750634673268802484271D0,&
 a16_6=-1.0459489594088330609505006875640990513158812317237848928608D0,&
 a16_7=0.930327845415626983292300564432428777137601651182965794680397D0,&
 a16_8=1.77950959431708102446142106794824453926275743243327790536D0,&
 a16_9=0.1D0,&
 a16_10=-0.282547569539044081612477785222287276408489375976211189952877D0,&
 a16_11=-0.159327350119972549169261984373485859278031542127551931461821D0,&
 a16_12=-0.145515894647001510860991961081084111308650130578626404945571D0,&
 a16_13=-0.259111214548322744512977076191767379267783684543182428778156D0,&
 a16_14=-0.34275815984718983994222055341385087174233873470395891993726D0,&
 a16_15=-0.675D0
! The estimate of the local truncation error is  (1/360)  h ( f(t1,x1)-f(t15,x15) )
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !

CONTAINS

SUBROUTINE rk108FeaginEachStep(NEQ,y0,yn,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
                                      y10,y11,y12,y13,y14,y15,y16
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: dy0,dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8, &
                                      dy9,dy10,dy11,dy12,dy13,dy14,dy15,dy16
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

  y13=y0+h*(a13_0*dy0+a13_1*dy1+a13_2*dy2+a13_3*dy3+a13_4*dy4+a13_5*dy5 + &
         &  a13_6*dy6+a13_7*dy7+a13_8*dy8+a13_9*dy9+a13_10*dy10+a13_11*dy11 + &
         &  a13_12*dy12)
  ! use y13 to get dy13
  CALL geo_eqns(NEQ,y13,C_t,C_phi,dy13)

  y14=y0+h*(a14_0*dy0+a14_1*dy1+a14_2*dy2+a14_3*dy3+a14_4*dy4+a14_5*dy5 + &
         &  a14_6*dy6+a14_7*dy7+a14_8*dy8+a14_9*dy9+a14_10*dy10+a14_11*dy11 + &
         &  a14_12*dy12+a14_13*dy13)
  ! use y14 to get dy14
  CALL geo_eqns(NEQ,y14,C_t,C_phi,dy14)

  y15=y0+h*(a15_0*dy0+a15_1*dy1+a15_2*dy2+a15_3*dy3+a15_4*dy4+a15_5*dy5 + &
         &  a15_6*dy6+a15_7*dy7+a15_8*dy8+a15_9*dy9+a15_10*dy10+a15_11*dy11 + &
         &  a15_12*dy12+a15_13*dy13+a15_14*dy14)
  ! use y15 to get dy15
  CALL geo_eqns(NEQ,y15,C_t,C_phi,dy15)

  y16=y0+h*(a16_0*dy0+a16_1*dy1+a16_2*dy2+a16_3*dy3+a16_4*dy4+a16_5*dy5 + &
         &  a16_6*dy6+a16_7*dy7+a16_8*dy8+a16_9*dy9+a16_10*dy10+a16_11*dy11 + &
         &  a16_12*dy12+a16_13*dy13+a16_14*dy14+a16_15*dy15)
  ! use y16 to get dy16
  CALL geo_eqns(NEQ,y16,C_t,C_phi,dy16)

  yn=y0+h*(b0*dy0+b1*dy1+b2*dy2+b3*dy3+b4*dy4+b5*dy5+b6*dy6+b7*dy7 + &
        &  b8*dy8+b9*dy9+b10*dy10+b11*dy11+b12*dy12+b13*dy13+b14*dy14+b15*dy15 + &
        &  b16*dy16)
  yerr = (1.D0/3.6D2)*h*ABS(dy1-dy15)
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
    hNew = MAX(0.8D0*err**(-1.D0/9.D0), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-1.D0/9.D0),MAX(0.8D0*err**(-1.D0/9.D0), ReduceAtMost)
  ELSE
    PleaseRerun = .False. ! the error is fine, keep this step and move on
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-1.D0/9.D0))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-1.D0/9.D0),MIN(IncreaseAtMost,0.8D0*err**(-1.D0/9.D0))
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
END SUBROUTINE rk108FeaginEachStep


! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! This subroutine choose the time step based on accuracy required and suggest
!   the next step
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
SUBROUTINE rk108Feagin(NEQ,X0,XN,Y0,C_t,C_phi,EPS,YN,IERR,h)
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
    CALL rk108FeaginEachStep(NEQ,Y0,YN,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
END SUBROUTINE rk108Feagin


END MODULE
