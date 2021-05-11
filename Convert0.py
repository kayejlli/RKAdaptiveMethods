from Python.RKGeneral import Printdy, PrintOutdys
from Python.Common import GetFortranFloat
import numpy as np
import sys


filename = sys.argv[1] 
stages = int(sys.argv[2])
order = int(sys.argv[3]) 

# python3 Convert0.py rk1412.py 35 12 Feagin 
#filename = 'rk1412.py' 
#stages = 35
#order = 12

try:
  modName = sys.argv[4] 
except IndexError:
  modName = filename.replace('.py','') 

f90 = open(modName+'.f90', 'w')

b_ = True
if filename in ['rk108.py','rk1210.py','rk1412.py','rk1412Long.py']:
  b_ = False

Head = '\
MODULE %sMod\n\
\n\
USE GlobalCommonMod\n\
USE DyDtMod\n\
\n\
IMPLICIT NONE\n\
PRIVATE\n\
\n\
!  Define access to SUBROUTINEs.\n\
\n\
PUBLIC :: %sEachStep\n' % (modName, modName) 
f90.write(Head) 

error = False
with open(filename,'r') as f:
  lines = f.readlines()
  Main = False 
  aStart = False 
  bStart = False 
  kStart = False 
  for line in lines:
    # print('\n' in line[:-1])
    if len(line[:-1]) > 1: # non-empty line  
      if line[0] == '#':
        if 'The estimate of the local truncation error' in line:
          error = line[1:-1].replace('The estimate of the local truncation error is', '') 
        # print(line[:-1]) 
        i = 1
        while line[i] == ' ': # trim spaces into one space
          i += 1				
        out = '! ' + line[i:-1]  
        # print('! ' + line[i:]) 
        # f90.write('! ' + line[i:]) 
        if line[i] == 'k':
          if not Main:
            #print('\n'*2) 
            f90.write('\n'*2) 
          Main = True 
          #print('REAL(KIND=8), PARAMETER, PRIVATE :: &') 
          #f90.write('REAL(KIND=8), PARAMETER, PRIVATE :: &\n') 
        else:
          f90.write('! ' + line[i:])
          # print('! ' + line[i:-1]) 
      elif line[0] in ['a', 'b', 'k'] and '=' in line:
        if not aStart and (line[0] == 'a'):
          f90.write('! a[i,j] \n') 
          f90.write('REAL(KIND=8), PARAMETER, PRIVATE :: &\n') 
          aStart = True
        if not bStart and (line[0] == 'b'):
          f90.write('! b[i] are coefficients for nodes k[i] \n') 
          f90.write('REAL(KIND=8), PARAMETER, PRIVATE :: &\n') 
          bStart = True
        if not kStart and (line[0] == 'k'):
          f90.write('! k[i] are the nodes \n') 
          f90.write('REAL(KIND=8), PARAMETER, PRIVATE :: &\n') 
          kStart = True
        left, right = line[:-1].split('=') 
        out = '%s=%s' % (left, GetFortranFloat(right.strip(), Precision=60, fortran=True))  
        if left[0]=='a' and int(left[1:].split('_')[1]) == stages-2:
          End = '\n'
        elif left[0]=='k' and int(left[1:]) == stages-1:
          End = '\n'
        elif left[0:2]=='b_' and int(left[2:]) == stages-1:
          End = '\n'
        elif not b_ and left[0]=='b' and int(left[1:]) == stages-1:
          End = '\n'
        else:
          End = ',&\n' 
        f90.write(' %s%s' % (out, End)) 
       
SubStart='\n\n\
CONTAINS\n\
\n\
SUBROUTINE %sEachStep(y0,yn,h,hnew,rerun,test)\n\
 REAL(KIND=8), DIMENSION(:), INTENT(IN) :: y0     ! y(t)\n\
 REAL(KIND=8), DIMENSION(SIZE(y0)), INTENT(OUT) :: yn ! y(t+h)\n\
 REAL(KIND=8), INTENT(IN) :: h ! initial step size\n\
 REAL(KIND=8), INTENT(OUT) :: hnew       ! new step size\n\
 LOGICAL, INTENT(OUT) :: test, rerun ! test=stop the program, rerun=re-run this step, reject\n\
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: tolh ! the tolerance, determined by the problem\n\
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: yerr ! the error between embedded method\n\
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ynp  ! the embedded\n\
 REAL(KIND=8), DIMENSION(SIZE(y0)) :: ymax ! the max value among y0 and yn\n\
 REAL(KIND=8) :: err\n\
' % (modName) 
f90.write(SubStart) 

breakNo = 9
Printdy(stages,breakNo,f90,mode='y')  
Printdy(stages,breakNo,f90,mode='dy')  
# declare more var 
print(' INTEGER :: i', file=f90)
print('  test = .False.', file=f90)
print('  rerun = .False. ', file=f90)

if error:
  PrintOutdys(f90, Msize=stages, devString = 'dev(y0,dy0,test)', breakNo = 6, test=True, yerr=False)
else:
  PrintOutdys(f90, Msize=stages, devString = 'dev(y0,dy0,test)', breakNo = 6, test=True)

def IntToFloat(intForm):
  raw = '%e' % (intForm)
  first, second = raw.split('e') 
  while first[-1] == '0': # trim the 0
    first = first[:-1] 
  return first + 'D' + '%d' % (int(second)) 

def ConvertFraction(form):
  numerator, denominator = form.split('/') 
  return '(%s/%s)' % (IntToFloat(float(numerator)), IntToFloat(float(denominator))) 

# print error 
if error:
  error = error.replace(' ','') 
  coef, form = error.split('h')
  print('coef=[%s], form=[%s]' % (coef, form)) 
  coef = ConvertFraction(coef.replace('(','').replace(')',''))
  No1 = int(form.split('-')[0].split('x')[-1].replace(')',''))
  No2 = int(form.split('-')[1].split('x')[-1].replace(')','')) 
  print('  yerr = %s*h*ABS(dy%d-dy%d)' % (coef, No1, No2), file=f90) 


Exp1 = '1.D0/%d.D0' % (order+1) 
# Exp1 = '%20.15e' % (1./(order+1))  
# Exp1 = Exp1.replace('e', 'D') 
ErrorAndTimeStep = "\
  ! Find the max value of y among this step\n\
  DO i = 1, SIZE(y0)\n\
    ymax(i) = MAX(MAX(ABS(y0(i)), ABS(yn(i))), h*ABS(dy0(i)))\n\
  END DO\n\
  tolh = rtol*ymax + atol(1:SIZE(y0)) ! atol might be a longer array\n\
\n\
  ! using the error to estimate the next step\n\
  err = MAXVAL(ABS(yerr/tolh))\n\
  IF (err.GT.1.D0) THEN\n\
    rerun = .True.\n\
    hnew = MAX(0.8D0*err**(-%s), 0.1D0)*h ! no less than factor of 0.1\n\
    ! PRINT *, 'Decrease time step by', 0.8D0*err**(-%s),MAX(0.8D0*err**(-%s), 0.1D0)\n\
  ELSE\n\
    rerun = .False.\n\
    hnew = MIN(5.D0, 0.8D0*err**(-%s))*h ! no more than factor of 5\n\
    ! PRINT *, 'Increase time step by', 0.8D0*err**(-%s),MIN(5.D0,0.8D0*err**(-%s))\n\
  END IF\n\
\n\
  ! adjust the step\n\
  hnew = MAX(MIN(hnew,MaxStepSize),MinStepSize)\n\
  IF (rerun) THEN\n\
    IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min\n\
      test = .True. ! stop the program\n\
    ELSE\n\
      hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half\n\
    END IF\n\
    RETURN\n\
  END IF\n\
\n\
  ! check if any value have went crazy (Nan or Inf)\n\
  DO i = 1, SIZE(y0)\n\
    CALL checkNanInf(yn(i), rerun)\n\
    IF (rerun) THEN\n\
      IF (ABS(h-MinStepSize)/MinStepSize.LE.1D-13) THEN ! h is already the min\n\
        test = .True. ! stop the program\n\
      ELSE\n\
        hnew = MAX(MinStepSize,MIN(hnew, h*0.5D0)) ! if have not reduced, reduce to half\n\
      END IF\n\
      RETURN\n\
    END IF\n\
  END DO\n\
\n\
  RETURN\n\
END SUBROUTINE\n\
\n\
END MODULE\n\
" % (Exp1, Exp1, Exp1, Exp1, Exp1, Exp1) 

f90.write(ErrorAndTimeStep) 






f90.close()
 
