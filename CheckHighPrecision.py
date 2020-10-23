#  this script check the coef up to 60 digits 
from RKGeneral import Printdy, PrintOutdys
import numpy as np
import sys
from Common import GetFortranFloat, rmrf

for filename in ['rk1412.py', 'rk1210.py', 'rk108.py', 'RK8713M.py', 'rk109_21.py', 'rk1211_31.py', 'rk65Dormand8.py', 'rk87EnrightVerner_13.py','rk1412Long.py','rk54Dormand7.py', 'rk54Sharp7.py']:
  print('--'*10) 
  print('------ %s ------' % (filename)) 
  new = filename.replace('.py', '_HighPrecision.py')
  f1 = open(new, 'w')
  
  print('from mpmath import mp, mpf, nstr', file=f1)
  print('mp.dps = 100', file=f1)
  
  with open(filename, 'r') as f0:
    lines = f0.readlines()
    for line in lines:
      if len(line[:-1]) <= 1 or line[0] == '#': 
        out = line[:-1] 
      elif 'Print' in line[:-1]:
        out = line[:-1].replace('())', '(), mpmode=True)') 
      elif line[0] in ['a', 'b', 'k']:
        left, right = line[:-1].split('=')
        if '/' in right:
          out = "%s=mpf('%s')" % (left.strip(), GetFortranFloat(right, Precision=60, fortran=False)) 
        else:
          out = "%s=mpf('%s')" % (left.strip(), right.strip())
      elif not 'print' in line: 
        out = line[:-1] 
      else:
        line = line[:-1].replace(' ','')
        left = line.split('%(')[0][7:-1].replace('%e', '%s') 
        right = line.split('%(')[1][:-2].split(',')
        List = [left] + ['nstr(%s, 10)' % element for element in right]
        out = 'print(\'%s\' %% (%s, %s, %s))' % (tuple(List))
      print(out, file=f1)
  
  f1.close()
  __import__(filename.replace('.py', '_HighPrecision'))
  rmrf(new) 
  
