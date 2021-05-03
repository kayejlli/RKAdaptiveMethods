MODULE ODEInterfaceMod

USE DyDtMod
USE GlobalCommonMod 

USE rk54SharpMod
USE rk54DormandMod
USE rk65DormandMod
USE rk87DormandMod
USE rk87EnrightVernerMod
USE rk108FeaginMod
USE rk109LegendreMod
USE rk1210FeaginMod
USE rk1211PeterMod
USE rk1412FeaginMod

IMPLICIT NONE 

CONTAINS

SUBROUTINE OUTPUT(t, yArray, UnitNo)
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: yArray
  REAL(KIND=8), INTENT(IN) :: t
  INTEGER, INTENT(IN) :: UnitNo
  ! general output format 
  WRITE(UNIT=UnitNo, FMT='(50E35.25)') t, yArray
  RETURN 
END SUBROUTINE 

SUBROUTINE Search(y0,rArray,thetaArray,phiArray,rHitboxInput,HitboxtolInput,aInput,SolverName,filenames,mode,atolInput,rtolInput,max_h,min_h,hinit,Print6Input,IntegerOut,RealOut) 

! SUBROUTINE ODE(y0,t0,tfinal,SolverName,filename,atolInput,rtolInput,max_h,min_h,hinit,Print6Input,IntegerOut,RealOut)
  CHARACTER(LEN=*), INTENT(IN) :: mode, SolverName
  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: filenames
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: rArray,thetaArray,phiArray 
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0, atolInput
  REAL(KIND=8), INTENT(IN) :: t0, tfinal, rtolInput, max_h, min_h, hinit 
  REAL(KIND=8), INTENT(IN) :: rHitboxInput,HitboxtolInput,aInput 
  INTEGER, DIMENSION(10), INTENT(OUT) :: IntegerOut
  REAL(KIND=8), DIMENSION(10), INTENT(OUT) :: RealOut 
  REAL(KIND=8) :: CPUTime0, CPUTime1 
  LOGICAL, INTENT(IN) :: Print6Input
  LOGICAL :: test

  ! assign values 
  atol(1:SIZE(y0)) = atolInput(:)  
  rtol = rtolInput
  Print6 = Print6Input ! whether print out to the screen 
  MaxStepSize = max_h
  MinStepSize = min_h

  ! size of the sphere and its tolerance 
  rHitbox = rHitboxInput
  Hitboxtol = HitboxtolInput

  ! spin of the bh
  a = aInput 

  ! initialise global values 
  Rejected = 0
  Accepted = 0
  Evaluated = 0
  TotalSteps = 0
  ReachMax = 0
  ReachMin = 0
  Minh = 1D14
  Maxh = 1D-14
  IntegerOut = 0

  ! CALL RunANDWrite(t0,tfinal,y0,ODEMode,filename,rtolInput,test,SolverName,StepPerOrbit) 
  CALL CPU_Time(CPUTime0) 
  CALL RunANDWrite(t0,tfinal,y0,SolverName,filename,test,hinit) 
  CALL RunANDWrite(y0,rArray,thetaArray,phiArray,SolverName,filenames,mode,test,hinit)
  CALL CPU_Time(CPUTime1) 

  ! assign values 
  RealOut(1) = CPUTime1 - CPUTime0 
  RealOut(2) = Minh
  RealOut(3) = Maxh 

  ! assign values 
  IntegerOut(1) = Rejected
  IntegerOut(2) = Accepted
  IntegerOut(3) = Evaluated
  IntegerOut(4) = TotalSteps
  IntegerOut(5) = ReachMax
  IntegerOut(6) = ReachMin
  IntegerOut(7) = 1 ! Sucess 
  
  IF (test) THEN
    IntegerOut(7) = 0 ! failed
    PRINT *, 'The program ended unexpectedly' 
  END IF
  RETURN
END SUBROUTINE

SUBROUTINE RunANDWrite(y0,rArray,thetaArray,phiArray,SolverName,filenames,mode,test,hinit)
  LOGICAL, INTENT(OUT) :: test
  CHARACTER(LEN=*), INTENT(IN) :: mode, SolverName
  CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: filenames
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0
  REAL(KIND=8), INTENT(IN) :: t0, tfinal, hinit ! initial time step
  REAL(KIND=8), DIMENSION(SIZE(y0)) :: yArray, y1Array, ydotNew, ydotOld
  ! INTEGER :: IterationCounts, i
  LOGICAL :: GoodToWrite, rerun
  REAL(KIND=8) :: time, PrintTimePassed, PrintAtTime, TimePassed, TimeStep, h, hnew, &  
  & hold
  PROCEDURE(ODEMethods), POINTER :: iterations => NULL() ! choose different solvers 
  CHARACTER(LEN=*), PARAMETER :: screen_fmt = "(1X,'t=',F5.2,' x1=',E18.10,' y1=',E18.10,' NSTEP=',I5)"  
  INTEGER :: Iteration, i, rSize

  rSize = SIZE(rArray)  
  
   ! w-> new file, a-> append only
   IF (mode.EQ.'w') THEN
     DO i = 1, rSize 
       OPEN(UNIT=30+i, FILE=filenames(i), STATUS="REPLACE", ACTION="WRITE")
       CLOSE(UNIT=30+i)
     END DO
   END IF
   DO i = 1, rSize 
     OPEN(UNIT=30+i, FILE=filenames(i), STATUS="OLD", ACCESS="SEQUENTIAL",&
        ACTION="WRITE", POSITION="APPEND")
   END DO


   IF (SolverName=='rk54Sharp') THEN
     iterations => rk54SharpEachStep
   ELSE IF (SolverName=='rk54Dormand') THEN
     iterations => rk54DormandEachStep
   ELSE IF (SolverName=='rk65Dormand') THEN
     iterations => rk65DormandEachStep
   ELSE IF (SolverName=='rk87Dormand') THEN
     iterations => rk87DormandEachStep
   ELSE IF (SolverName=='rk87EnrightVerner') THEN
     iterations => rk87EnrightVernerEachStep
   ELSE IF (SolverName=='rk108Feagin') THEN
     iterations => rk108FeaginEachStep
   ELSE IF (SolverName=='rk109Legendre') THEN
     iterations => rk109LegendreEachStep
   ELSE IF (SolverName=='rk1210Feagin') THEN
     iterations => rk1210FeaginEachStep
   ELSE IF (SolverName=='rk1211Peter') THEN
     iterations => rk1211PeterEachStep
   ELSE IF (SolverName=='rk1412Feagin') THEN
     iterations => rk1412FeaginEachStep
   ELSE
     PRINT '(3A)', 'Error: You should choose between [rk54Sharp,rk54Dormand,&
 &rk65Dormand,rk87Dormand,rk87EnrightVerner,rk108Feagin,rk109Legendre,&
 &rk1210Feagin,rk1211Peter,rk1412Feagin] but received [', SolverName, ']' 
   END IF

   ! Initialise 
   yArray(:) = y0(:)
   time = t0 
   h = hinit 

   ! IterationCounts = 0
   test = .False.
   rerun = .False. 
   PrintTimePassed = 0.D0
   TimePassed = 0.D0
   PrintAtTime = tfinal/2D1 ! print about 20 steps
   ! Print data to file at dt = TimeStep 
   TimeStep = tfinal/5.123124D5

   ! MaxStepSize = (2.D0*Pi*ai**1.5)/Sqrt(Mi) / 1D-2
   ! MinStepSize = (2.D0*Pi*ai**1.5)/Sqrt(Mi) / 1D4
   ! MaxStepSize = TimeStep/5.D0
   ! MinStepSize = h*1D-3
   ! PRINT *, h, MaxStepSize, MinStepSize
  
   ! CALL dev(yArray, dot
   ! CALL Getdy(yArray, ydotOld, ODEMode, test, Energy, Jtot)
   CALL dev(yArray, ydotOld, test) 
   IF (LEN(filename)>0) THEN
     CALL OUTPUT(time, yArray, 30)
   END IF

   DO
     ! PRINT *, '---------------'
     Iteration = 0
     DO ! get better time step   
       hold = h
       CALL iterations(yArray, y1Array, h, hnew, rerun, test) 
       IF (test) EXIT       ! Bad value, exit the program 
       ! IF (.NOT.rerun) EXIT ! do not need to rerun 
       IF (rerun) THEN
         Rejected = Rejected + 1
       ELSE
         IF (Iteration==0) Accepted = Accepted + 1
         EXIT
       END IF
       ! PRINT *, 'updating time step',h, hnew, yArray(2), y1Array(2) 
       h = hnew
       Iteration = Iteration + 1
       ! exit when step size is bad 
       IF ((Abs(hnew-hold)/hold).LT.1D-13) test = .True.
       IF ((Abs(hold-MinStepSize)/MinStepSize).LE.1D-13) test = .True.
       IF (test) EXIT       ! Bad value, exit the program 
     END DO
     IF (test) EXIT ! Bad value, exit the program 
     
     ! calculate the max and min time step reached, and how many times 
     Minh = MIN(Minh, hold)
     Maxh = MAX(Maxh, hold) 
     IF (ABS(hold-MinStepSize)/MinStepSize.LE.1D-12) THEN
       ReachMin = ReachMin + 1
     ELSE IF (ABS(hold-MaxStepSize)/MaxStepSize.LE.1D-12) THEN
       ReachMax = ReachMax + 1
     END IF


     ! get dy/dt 
     CALL dev(y1Array, ydotNew, test) 
     IF (test) EXIT ! Bad value, exit the program 

     CALL CheckingHitting(yArray,ydotOld,y1Array,ydotNew,filenames,Hit)  
     IF (Hit.EQ.0) THEN

     time = time + hold
     TimePassed = TimePassed + hold 
     PrintTimePassed = PrintTimePassed + hold
     TotalSteps = TotalSteps + 1
     ydotOld = ydotNew ! save dy/dt in case needed 

     ! update values
     h = hnew
     yArray = y1Array
     GoodToWrite = .False.




     IF (TimePassed.GE.TimeStep) THEN
       GoodToWrite = .True.
       TimePassed = TimePassed - TimeStep
     END IF

     ! print to screen and data 
     IF (GoodToWrite) THEN
       IF (LEN(filename)>0) THEN
         CALL OUTPUT(time, y1Array, 30)
       END IF
       IF ((PrintTimePassed.GE.PrintAtTime).AND.Print6) THEN
         WRITE(6,screen_fmt) time, y1Array(1), y1Array(8), TotalSteps 
         PrintTimePassed = PrintTimePassed - PrintAtTime
       END IF
     END IF

     IF (time+h.GT.tfinal) THEN
       IF (time.GE.tfinal) EXIT 
       h = tfinal - time
       IF (h<1D-15) EXIT
     END IF
  END DO
   
  ! print the last data anyway 
  IF ((.NOT.GoodToWrite).AND.(.NOT.test).AND.(LEN(filename)>0)) THEN
    CALL OUTPUT(time, y1Array, 30)
  END IF

  IF (Print6) THEN
    PRINT *, '--- ENDING ---'
    WRITE(6,screen_fmt) time, y1Array(1), y1Array(8), TotalSteps 
    91 FORMAT(' fcn=',I6,' step=',I6,' accpt=',I6, ' rejct=',I6,' htkmax=',I5,' htkmin=',I5) 
    WRITE(6,91) Evaluated, TotalSteps, Accepted, Rejected, ReachMax, ReachMin 
    PRINT '(A,E11.3,A,E11.3)', ' Maxh=', Maxh, ' UserInput=',  MaxStepSize 
    PRINT '(A,E11.3,A,E11.3)', ' Minh=', Minh, ' UserInput=',  MinStepSize 
    PRINT '(3A)', ' Data saved in [',filename,']' 
  END IF
  
  IF (LEN(filename)>0) THEN
    CLOSE(UNIT=30)
  END IF
END SUBROUTINE RunANDWrite

END MODULE 
