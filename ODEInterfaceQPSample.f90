MODULE ODEInterfaceQPMod

USE ODEInterfaceMod

IMPLICIT NONE 

CONTAINS


SUBROUTINE ODEQP(SolverName,filename,Print6Input,IntegerOut,RealOut,rtolDP,atolDP)
  CHARACTER(LEN=*), INTENT(IN) :: filename, SolverName
  LOGICAL, INTENT(IN) :: Print6Input
  INTEGER, DIMENSION(10), INTENT(OUT) :: IntegerOut
  REAL(KIND=8), DIMENSION(10), INTENT(OUT) :: RealOut 
  REAL(KIND=8), DIMENSION(28), INTENT(IN) :: atolDP
  REAL(KIND=8), INTENT(IN) :: rtolDP
  
  ! local numbers 
  REAL(KIND=16) :: t0, tfinal, rtolInput, max_h, min_h, hinit 
  REAL(KIND=16), DIMENSION(28) :: y0, atolInput
  INTEGER :: i

----------------- include here ---------------


  ! convert rtol and atol to QP 
  DO i = 1, 28
    atolInput(i) = atolDP(i)*1.Q0
  END DO
  rtolInput = rtolDP*1.Q0 
  PRINT *, 'atolInput', atolInput(1) 
  PRINT *, 'rtolInput', rtolInput 

  CALL ODE(y0,t0,tfinal,SolverName,filename,atolInput,rtolInput,max_h,min_h,hinit,Print6Input,IntegerOut,RealOut)

END SUBROUTINE

END MODULE 
