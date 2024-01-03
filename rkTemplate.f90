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
USE GlobalCommonMod
USE DyDtMod  

IMPLICIT NONE
PRIVATE

!  Define access to SUBROUTINEs.

PUBLIC :: {modName}EachStep

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
SUBROUTINE {modName}EachStep(t,y0,yn,h,hNew,PleaseRerun,PleaseTerminate)
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
    hNew = MAX(0.8D0*err**(-{myExp}), ReduceAtMost)*h ! no less than factor of ReduceAtMost
    ! PRINT *, "Decrease time step by", 0.8D0*err**(-{myExp}),MAX(0.8D0*err**(-{myExp}), ReduceAtMost)
  ELSE
    ! PleaseRerun = .False. ! the error is fine, keep this step and move on
    ! IncreaseAtMost is suggested to be 5.D0
    hNew = MIN(IncreaseAtMost, 0.8D0*err**(-{myExp}))*h ! no more than factor of IncreaseAtMost
    ! PRINT *, "Increase time step by", 0.8D0*err**(-{myExp}),MIN(IncreaseAtMost,0.8D0*err**(-{myExp}))
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
END SUBROUTINE {modName}EachStep


END MODULE {modName}Mod
