MODULE GlobalCommonMod 
! This module save global common variables, and is accessiable for all subroutines

IMPLICIT NONE 

SAVE 

ABSTRACT INTERFACE
  SUBROUTINE ODEMethods(y0,y1,h,hnew,rerun,test)
  REAL(KIND=16), DIMENSION(:), INTENT(IN) :: y0
  REAL(KIND=16), DIMENSION(SIZE(y0)), INTENT(OUT) :: y1
  REAL(KIND=16), INTENT(IN) :: h
  LOGICAL, INTENT(OUT) :: test, rerun ! if test, quit program, if rerun, change h
  REAL(KIND=16), INTENT(OUT) :: hnew
  END SUBROUTINE
END INTERFACE


  INTEGER,PARAMETER :: QP = selected_real_kind(33, 4931) ! Quad Precision
  REAL(KIND=QP), PARAMETER :: ErrorReference = 1Q-15
  REAL(KIND=QP), PARAMETER :: G = 1.0_QP, c = 1.0_QP 
  REAL(KIND=QP), PARAMETER :: Pi = 3.1415926535897932384626433832795028841_QP 
  REAL(KIND=QP) :: rtol, MinStepSize, MaxStepSize  
  ! REAL(KIND=QP) :: NumberPerOrbit, hInitial, Period0, MinStepSize, MaxStepSize
  REAL(KIND=QP), DIMENSION(50) :: atol
  INTEGER :: Rejected, Accepted, Evaluated, TotalSteps, ReachMax, ReachMin 
  REAL(KIND=QP) :: Minh, Maxh 
  LOGICAL :: Print6

END MODULE GlobalCommonMod
