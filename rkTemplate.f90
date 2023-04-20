!# -------------------------------------------------------------
!# -------------------------------------------------------------
!# Comments that starts with "!#" sign will be removed by the python script
!# Comments that starts with "!" will be kept in the final fortran code  
!# Warning: Single quotation marks is NOT allowed in this template,
!#          Please use double quotation marks only
!# -------------------------------------------------------------
!# -------------------------------------------------------------
MODULE {modName}Mod

!# define the modules that your fortran script needs  
USE parameters
USE routines

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: {modName}

! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! --------[start coefficients definitions] (do not touch this line) ------ !
!# anything between [start] and [end] will be removed by python script 
!# there is no need to modify this part, the python script will change these to
!#    the k0,k1,a1_0 to the correct one 
!# here are some examples ...
REAL(KIND=8), PARAMETER, PRIVATE :: &
 k0=0.0D0,&
 k1=0.2D0,&
 etc...
! a[i,j] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 a1_0=0.2D0,&
 a2_0=7.5D-2,&
 a2_1=0.225D0,&
 a3_0=0.977777777777777777777777777777777777777777777777777777777778D0,&
 etc...
! b[i] are coefficients for nodes k[i] 
REAL(KIND=8), PARAMETER, PRIVATE :: &
 b0=9.11458333333333333333333333333333333333333333333333333333333D-2,&
 b1=0.0D0,&
 etc.....
! -------- [end coefficients definitions] (do not touch this line) ------- !
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !

CONTAINS

!# {modName} --> this will be replaced by the python script to, e.g. rk108Feagin
SUBROUTINE {modName}EachStep(NEQ,y0,yn,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
 ! ---- [start construct intermediate steps] (do not touch this line) ----- !
 !# anything between [start] and [end] will be removed by python script 
 !# there is no need to modify this part, the python script will change these to
 !# 
 !#
  CALL  dev(t,y0,dy0,PleaseTerminate)
  y1=y0+h*(a1_0*dy0)
  ! use y1 to get dy1
  CALL  dev(t,y1,dy1,PleaseTerminate)
  y2=y0+h*(a2_0*dy0+a2_1*dy1)
  ! use y2 to get dy2
  CALL  dev(t,y2,dy2,PleaseTerminate)
  ..... and then get y3, y4, etc 

  yn=y0+h*(b0*dy0+b1*dy1+...)
  yerr=h*((b0-b_0)*dy0+(b1-b_1)*dy1...)

 ! ----- [end construct intermediate steps] (do not touch this line) ------ !
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
    hNew = MAX(0.8D0*err**(-{myExp}), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-{myExp}),MAX(0.8D0*err**(-{myExp}), ReduceAtMost)
  ELSE
    PleaseRerun = .False. ! the error is fine, keep this step and move on 
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-{myExp}))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-{myExp}),MIN(IncreaseAtMost,0.8D0*err**(-{myExp}))
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
END SUBROUTINE {modName}EachStep


! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
! This subroutine choose the time step based on accuracy required and suggest
!   the next step 
! ------------------------------------------------------------------------ !
! ------------------------------------------------------------------------ !
SUBROUTINE {modName}(NEQ,X0,XN,Y0,C_t,C_phi,EPS,YN,IERR,h)
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
    CALL {modName}EachStep(NEQ,Y0,YN,h,hNew,EPS,C_t,C_phi,PleaseRerun,PleaseTerminate)
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
END SUBROUTINE {modName}


END MODULE
