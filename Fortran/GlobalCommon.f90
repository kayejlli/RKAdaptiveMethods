MODULE GlobalCommonMod 
! This module save global common variables, and is accessiable for all subroutines

USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NORMAL

IMPLICIT NONE 

SAVE 

ABSTRACT INTERFACE
  ! Do not touch this ABSTRACT INTERFACE unless you know what are you doing. 
  SUBROUTINE ODEMethods(t,y0,y1,h,hnew,rerun,test)
  REAL(KIND=8), INTENT(IN) :: t
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0
  REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: y1
  REAL(KIND=8), INTENT(IN) :: h
  LOGICAL, INTENT(OUT) :: test, rerun ! if test, quit program, if rerun, change h
  REAL(KIND=8), INTENT(OUT) :: hnew
  END SUBROUTINE
END INTERFACE


  INTEGER,PARAMETER :: DP = selected_real_kind(15, 307) ! Double Precision 
  REAL(KIND=DP), PARAMETER :: ErrorReference = 1E-15_DP
  REAL(KIND=DP), PARAMETER :: G = 1.0_DP, c = 1.0_DP 
  REAL(KIND=DP), PARAMETER :: Pi = 3.1415926535897932384626433832795028841_DP 
  REAL(KIND=DP), PARAMETER :: ReduceAtMost = 0.1D0 ! reduce the step size by at most 0.1
  REAL(KIND=DP), PARAMETER :: IncreaseAtMost = 5.D0 ! increase the step size by at most 5 
  REAL(KIND=DP) :: rtol, MinStepSize, MaxStepSize  
  ! REAL(KIND=DP) :: NumberPerOrbit, hInitial, Period0, MinStepSize, MaxStepSize
  REAL(KIND=DP), DIMENSION(50) :: atol
  INTEGER :: Rejected, Accepted, Evaluated, TotalSteps, ReachMax, ReachMin 
  REAL(KIND=DP) :: Minh, Maxh 
  LOGICAL :: Print6

REAL(KIND=DP) :: a ! Spin of BH
REAL(KIND=DP) :: rT ! torus size
REAL(KIND=DP) :: rEH ! Event horizon 
REAL(KIND=DP) :: phiMax ! max value allowed for phi

END MODULE GlobalCommonMod
