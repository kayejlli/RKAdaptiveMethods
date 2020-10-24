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

  y0(1) = 3.0Q+0
  y0(2) = 3.0Q+0
  y0(3) = -1.0Q+0
  y0(4) = -3.0Q+0
  y0(5) = 2.0Q+0
  y0(6) = -2.0Q+0
  y0(7) = 2.0Q+0
  y0(8) = 3.0Q+0
  y0(9) = -3.0Q+0
  y0(10) = 2.0Q+0
  y0(11) = 0.0Q+0
  y0(12) = 0.0Q+0
  y0(13) = -4.0Q+0
  y0(14) = 4.0Q+0
  y0(15) = 0.0Q+0
  y0(16) = 0.0Q+0
  y0(17) = 0.0Q+0
  y0(18) = 0.0Q+0
  y0(19) = 0.0Q+0
  y0(20) = 1.75Q+0
  y0(21) = -1.5Q+0
  y0(22) = 0.0Q+0
  y0(23) = 0.0Q+0
  y0(24) = 0.0Q+0
  y0(25) = -1.25Q+0
  y0(26) = 1.0Q+0
  y0(27) = 0.0Q+0
  y0(28) = 0.0Q+0
  t0 = 0.0Q+0
  tfinal = 3.0Q+0
  max_h = 1.0Q+0
  min_h = 1.0Q-10
  hinit = 1.0Q-3


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
