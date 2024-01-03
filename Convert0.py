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


b_ = True
if filename.split('/')[1] in ['rk108.py','rk1210.py','rk1412.py','rk1412Long.py']:
  b_ = False

real = True # for complex y & dy

# the lines that will be written in between [start coefficients definitions]
#     and [end coefficients definitions]
coefficients_definitions = []


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
        # coefficients_definitions.append('! ' + line[i:])
        if line[i] == 'k':
          if not Main:
            #print('\n'*2) 
            coefficients_definitions.append('\n'*2)
          Main = True 
          #print('REAL(KIND=8), PARAMETER, PRIVATE :: &') 
          #coefficients_definitions.append('REAL(KIND=8), PARAMETER, PRIVATE :: &\n')
        else:
          coefficients_definitions.append('! ' + line[i:])
          # print('! ' + line[i:-1]) 
      elif line[0] in ['a', 'b', 'k'] and '=' in line:
        if not aStart and (line[0] == 'a'):
          coefficients_definitions.append('! a[i,j] \n')
          coefficients_definitions.append('REAL(KIND=8), PARAMETER, PRIVATE :: &\n')
          aStart = True
        if not bStart and (line[0] == 'b'):
          coefficients_definitions.append('! b[i] are coefficients for nodes k[i] \n')
          coefficients_definitions.append('REAL(KIND=8), PARAMETER, PRIVATE :: &\n')
          bStart = True
        if not kStart and (line[0] == 'k'):
          coefficients_definitions.append('! k[i] are the nodes \n')
          coefficients_definitions.append('REAL(KIND=8), PARAMETER, PRIVATE :: &\n')
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
        coefficients_definitions.append(' %s%s' % (out, End))
       

breakNo = 9

# the lines that will be written in between [start construct intermediate steps]
#     and [end construct intermediate steps]
construct_intermediate_steps = []

# determine if dy-s are real or not
if real:
  construct_intermediate_steps += Printdy(stages,breakNo,mode='y',Head=' REAL(KIND=8), DIMENSION(SIZE(y0)) ::')
  construct_intermediate_steps += Printdy(stages,breakNo,mode='dy',Head=' REAL(KIND=8), DIMENSION(SIZE(y0)) ::')
else:
  construct_intermediate_steps += Printdy(stages,breakNo,mode='y',Head=' COMPLEX(KIND=8), DIMENSION(SIZE(y0)) ::')
  construct_intermediate_steps += Printdy(stages,breakNo,mode='dy',Head=' COMPLEX(KIND=8), DIMENSION(SIZE(y0)) ::')

# declare more var 
construct_intermediate_steps.append(' INTEGER :: i\n')
construct_intermediate_steps.append('  ! To initialise the logical variables\n')
construct_intermediate_steps.append('  PleaseTerminate = .False.\n')
construct_intermediate_steps.append('  PleaseRerun = .False.\n')
construct_intermediate_steps.append('  ! ------------------------------------------------------------------------ !\n')
construct_intermediate_steps.append('  ! ------------------------------------------------------------------------ !\n')


if error:
  more = PrintOutdys(Msize=stages, devString = 'dev(t,y0,dy0,PleaseTerminate)', breakNo = 6, test=True, yerr=False)
else:
  more = PrintOutdys(Msize=stages, devString = 'dev(t,y0,dy0,PleaseTerminate)', breakNo = 6, test=True)
construct_intermediate_steps += more


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
  construct_intermediate_steps.append('  yerr = %s*h*ABS(dy%d-dy%d)\n' % (coef, No1, No2))

myExp = '1.D0/%d.D0' % (order+1)


f90 = open('Fortran/'+modName+'.f90', 'w')

skip = False 

# load the template 
with open('rkTemplate.f90', 'r') as f:
  linesTemplate = f.readlines()

for line in linesTemplate:
  ########################################
  # inset the lines at the right place
  if '[start coefficients definitions]' in line:
    skip = True
    for newline in coefficients_definitions:
      f90.write(newline)
    continue
  if '[end coefficients definitions]' in line:
    skip = False  
    continue
  ########################################
  # inset the lines at the right place
  if '[start construct intermediate steps]' in line:
    skip = True
    for newline in construct_intermediate_steps:
      f90.write(newline)
    continue
  if '[end construct intermediate steps]' in line:
    skip = False  
    continue
  # remove comments 
  if line.lstrip()[:2] == '!#': continue 
  if not skip:
    exec("newline = f'%s'" % (line.rstrip()))
    f90.write(newline + '\n')

f90.close()
