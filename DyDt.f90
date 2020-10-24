MODULE DyDtMod
! Pleiades problem

USE GlobalCommonMod 

IMPLICIT NONE 

CONTAINS

SUBROUTINE checkNanInf(test1,test) 
  REAL(KIND=QP), INTENT(IN) :: test1 
  LOGICAL, INTENT(INOUT) :: test
  test = .False.
  IF ((test1.NE.test1).OR.(test1-1.Q0.EQ.test1)) THEN
    test = .TRUE.
  END IF
  RETURN
END SUBROUTINE checkNanInf

SUBROUTINE dev(y,dy,test)
  REAL(KIND=16), DIMENSION(:), INTENT(IN) :: y
  REAL(KIND=16), DIMENSION(SIZE(y)), INTENT(OUT) :: dy
  LOGICAL, INTENT(OUT) :: test
  REAL(KIND=16) :: mj, dotx, doty, rij
  INTEGER :: i, j
  test = .False.

  dy = 0.Q0

  DO i = 1, 7 !  Particle No
    dotx = 0.Q0 
    doty = 0.Q0 
    DO j = 1, 7 ! the other particle no 
      mj = DFLOAT(j)    ! mass
      IF (j==i) CYCLE
      rij = Sqrt((y(i)-y(j))**2 + (y(i+7)-y(j+7))**2)  
      IF (rij.LE.1Q-15) THEN
        test = .True.
        RETURN 
      END IF
      dotx = dotx + mj*(y(j)-y(i))/rij**3
      doty = doty + mj*(y(j+7)-y(i+7))/rij**3
    END DO
    dy(i+14) = dotx
    dy(i+21) = doty 
  END DO

  ! assign dot(x) and dot(y) 
  dy(1:14) = y(15:28)
 
  Evaluated = Evaluated + 1
  RETURN 
END SUBROUTINE
 

END MODULE
