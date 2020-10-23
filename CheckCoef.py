# from CheckCoef import *
# Print(35,globals()) 
import numpy as np
from mpmath import *
mp.dps = 100 

def mynsum(Array):
  Sum = mpf('0.')
  for No in Array:
    Sum = Sum + No
  return Sum

def Print(Rank, InputDict, mpmode=False):
  bSum = mpf('0.') 
  kArray = []
  bArray = []
  for i in range(Rank):
    aSum = mpf('0.') 
    bSum += InputDict['k%d' % i]*InputDict['b%d' % i] 
    kArray.append(InputDict['k%d' % i])
    bArray.append(InputDict['b%d' % i]) 
    for j in range(i):
      aSum += InputDict['a%d_%d' % (i, j)] 
    err = aSum-InputDict['k%d' % i] 
    # print('k%d=%e, %e, dif=%e' % (i,InputDict['k%d' % i], aSum, abs(aSum-InputDict['k%d' % i]))) 
    # print('dif=%e, ratio=%e' % (abs(err), abs(err/aSum))) 
    eps = 1E-15
    if mpmode: 
      eps = 1E-58
    if abs(err) > eps:
      print('k%d, dif=%s' % (i, nstr(abs(err), 10))) 
  print('b*k=%e, err=%s' % (bSum,nstr(abs(bSum-mpf('0.5')), 10))) 
  print('bSum=%e, err=%s' % (mynsum(bArray), nstr(abs(mynsum(bArray)-mpf('1')), 10))) 
  return 

def Print2(Rank, InputDict, mpmode=False): 
  bSum = mpf('0.') 
  b_Sum = mpf('0.') 
  bArray = []
  b_Array = []
  for i in range(Rank):
    aSum = mpf('0.') 
    bSum += InputDict['k%d' % i]*InputDict['b%d' % i] 
    b_Sum += InputDict['k%d' % i]*InputDict['b_%d' % i] 
    bArray.append(InputDict['b%d' % i]) 
    b_Array.append(InputDict['b_%d' % i]) 
    for j in range(i):
      aSum += InputDict['a%d_%d' % (i, j)] 
    err = aSum-InputDict['k%d' % i] 
    # print('k%d=%e, %e, dif=%e' % (i,InputDict['k%d' % i], aSum, abs(aSum-InputDict['k%d' % i]))) 
    # print('dif=%e, ratio=%e' % (abs(err), abs(err/aSum))) 
    eps = 1E-15
    if mpmode: 
      eps = 1E-58
    if abs(err) > eps:
      print('k%d, dif=%s' % (i, nstr(abs(err), 10))) 
  print('b*k=%e, err=%s' % (bSum,nstr(abs(bSum-mpf('0.5')), 10))) 
  print('b_*k=%e, err=%s' % (b_Sum,nstr(abs(b_Sum-mpf('0.5')), 10))) 
  print('bSum=%e, err=%s' % (mynsum(bArray), nstr(abs(mynsum(bArray)-mpf('1')), 10))) 
  print('b_Sum=%e, err=%s' % (mynsum(b_Array), nstr(abs(mynsum(b_Array)-mpf('1')), 10))) 
  return 
