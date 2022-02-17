PROGRAM TestProgram 
! This code is a sample program, which demonstrate the usage of other modules      

USE ODEInterfaceMod ! the ODE interface

IMPLICIT NONE 

! y0         : initial condition (n-dimensional vector)
!              ['x1','x2','x3','x4','x5','x6','x7','y1','y2','y3','y4','y5','y6','y7',
!               'dx1','dx2','dx3','dx4','dx5','dx6','dx7','dy1','dy2','dy3','dy4','dy5','dy6','dy7']
! t0         : initial time (float)
! tfinal     : final time (float)
! SolverName : ode solver, should be in ['rk54Sharp','rk54Dormand','rk65Dormand',
!       'rk87Dormand','rk87EnrightVerner','rk108Feagin','rk109Legendre','rk1210Feagin','rk1211Peter','rk1412Feagin']
! filename   : format:'Data/*.dat', filename='' if you do not want to save the data file 
! atolInput  : an array (n-dimensional) 
! rtolInput  : float
! max_h      : maximum step size (float)
! min_h      : minimum step size (float)
! hinit      : initial time step (float)
! Print6Input : logical, True if you want to see output on screen

!--- OUTPUTs
! IntegerOut : (10-dimensional vector) with each number representing different meanings
!      IntegerOut(1): Rejected (how many time steps have been rejected due to time step being too large)
!      IntegerOut(2): Accepted (how many times the solver accept a time step without being rejected before)  
!      IntegerOut(3): Evaluated (how many times the y'=f(t,y) is Evaluated)  
!      IntegerOut(4): TotalSteps (the total number of time step (np.size(t) \\approx min(TotalSteps,5.123124E5))  
!      IntegerOut(5): ReachMax (how many times h reaches max_step, if this number is large, you may want to change max_step)  
!      IntegerOut(6): ReachMin (how many times h reaches min_step, if this number is large, you may want to change min_step)  
!      IntegerOut(7): Success (0 or 1 # 1=succeeded, 0=the program did not reach the end)
!      IntegerOut(8:10) : nothing, you can save your integers here!
! RealOut    : (10-dimensional vector) 
!      RealOut(1): cpuTime 
!      RealOut(2): Minh (the minimum time step used inside the program) 
!      RealOut(3): Maxh (the maximum time step used inside the program)

CHARACTER(LEN=*), PARAMETER :: filename = 'Data/test.dat' ! save data to this file
CHARACTER(LEN=*), PARAMETER :: SolverName = 'rk65Dormand' ! one of the RK method (check ODEInterface.f90 for details)
LOGICAL :: Print6Input = .True. ! print output to the screen
INTEGER, DIMENSION(10) :: IntegerOut
REAL(KIND=8), DIMENSION(10) :: RealOut
! local numbers 
REAL(KIND=8) :: t0, tfinal, rtolInput, max_h, min_h, hinit 
REAL(KIND=8), DIMENSION(28) :: y0, atolInput
INTEGER :: i

  


  ! initial values
  y0(1) = 3.0D+0
  y0(2) = 3.0D+0
  y0(3) = -1.0D+0
  y0(4) = -3.0D+0
  y0(5) = 2.0D+0
  y0(6) = -2.0D+0
  y0(7) = 2.0D+0
  y0(8) = 3.0D+0
  y0(9) = -3.0D+0
  y0(10) = 2.0D+0
  y0(11) = 0.0D+0
  y0(12) = 0.0D+0
  y0(13) = -4.0D+0
  y0(14) = 4.0D+0
  y0(15) = 0.0D+0
  y0(16) = 0.0D+0
  y0(17) = 0.0D+0
  y0(18) = 0.0D+0
  y0(19) = 0.0D+0
  y0(20) = 1.75D+0
  y0(21) = -1.5D+0
  y0(22) = 0.0D+0
  y0(23) = 0.0D+0
  y0(24) = 0.0D+0
  y0(25) = -1.25D+0
  y0(26) = 1.0D+0
  y0(27) = 0.0D+0
  y0(28) = 0.0D+0

  t0 = 0.0D+0
  tfinal = 3.0D+0
  max_h = 1.0D+0
  min_h = 1.0D-10
  hinit = 1.0D-3

  ! convert rtol and atol to DP 
  DO i = 1, 28
    atolInput(i) = 1D-8 
  END DO
  rtolInput = 1D-7

  CALL ODE(y0,t0,tfinal,SolverName,filename,atolInput,rtolInput,max_h,min_h,hinit,Print6Input,IntegerOut,RealOut)

END PROGRAM
