MODULE DyDtMod
! Pleiades problem

USE GlobalCommonMod 

IMPLICIT NONE 

CONTAINS


SUBROUTINE dev(t,y,dy,test)
  ! Please modify this file and use your own y' = f(t,y)
  ! t    : current time
  ! y(n) : current y (say, position of a particle in a binary system)
  ! dy(n) : dy/dt 
  ! test: logical, set it to .True. if you would like to force the solver to stop right now
  !       you may use it to detect Nan values or bad values

  REAL(KIND=8), INTENT(IN) :: t
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y
  REAL(KIND=8), DIMENSION(SIZE(y)), INTENT(OUT) :: dy
  LOGICAL, INTENT(OUT) :: test
  REAL(KIND=8) :: mj, dotx, doty, rij
  INTEGER :: i, j
  test = .False.

  dy = 0.D0

  DO i = 1, 7 !  Particle No
    dotx = 0.D0 
    doty = 0.D0 
    DO j = 1, 7 ! the other particle no 
      mj = DFLOAT(j)    ! mass
      IF (j==i) CYCLE
      rij = Sqrt((y(i)-y(j))**2 + (y(i+7)-y(j+7))**2)  
      IF (rij.LE.1D-15) THEN
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
