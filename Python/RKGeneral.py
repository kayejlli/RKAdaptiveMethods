# output *f90 file with given order 
import numpy as np
import sys
# from RKGeneral import Printdy, PrintOutdys 


def Printdy(jsize, breakNo, Head=' REAL(KIND=8), DIMENSION(SIZE(y0)) :: ', mode='y'):
  # addLines.append(Head) # REAL(KIND=DP), DIMENSION(SIZE(y0)) :: 
  addLines = []
  if mode == 'y':
    left = jsize - 1 
    maxSize = jsize - 1
    Add = 1
  elif mode == 'dy':
    left = jsize
    maxSize = jsize
    Add = 0
  start = 0
  if left < breakNo:
    addLines.append(Head + ','.join('%s%d' % (mode,j+Add) for j in range(left)) + '\n') 
    return addLines
  while left > 0: # more lines to write 
    if start > 0:
      Head = ' '* len(Head)
    End = ', &'
    if left - breakNo <= 0:
      End = ''
    line = ','.join('%s%d' % (mode,j+Add) for j in range(start, min(start+breakNo,maxSize)))      
    addLines.append(Head + line + End + '\n')
    left = left - breakNo
    start += breakNo
  return addLines 

def BreakLines(List, joinSign, FirstLineHead, LinkSign, EndString, breakNo, spaces):
  addLines = []
  jsize = np.size(List)  
  addSingleLine = '%s%s%s' % (FirstLineHead, joinSign.join(List[i] for i in range(jsize)), EndString) 
  # print('addSingleLine', addSingleLine) 
  left = jsize
  start = 0 
  printed = False
  if left < breakNo:
    addLines.append(spaces+addSingleLine + '\n') 
    printed = True
  addBreak = ''
  while left > 0: 
    line = '%s%s' % (FirstLineHead, joinSign.join(List[i] for i in range(start, min(start+breakNo,jsize))))   
    if left - breakNo <= 0:
      if not printed:
        addLines.append(spaces+line + EndString + '\n') 
      addBreak = addBreak + line + EndString
    else:
      addLines.append(spaces+line + LinkSign + '\n') 
      addBreak = addBreak + line + LinkSign 
    FirstLineHead = ' '*(len(FirstLineHead)-3) + '&  ' 
    left = left - breakNo 
    start += breakNo
  # check that these two matches !!! 
  addBreakClean = addBreak.replace(' ','').replace('&','')
  if not addBreakClean == addSingleLine:
    addLines.append('addSingleLine = [%s]' % addSingleLine + '\n') 
    addLines.append('addBreakClean = [%s]' % addBreakClean + '\n') 
    raise ValueError('add break does not match add') 
  return addLines


def PrintOutdys(Msize=13, devString = 'dev(y0,dy0,test)', breakNo = 6, spaceNo=2, test=False, yerr=True): 
  # if (test), return with smaller h
  addLines = []
  spaces = ' '*spaceNo 
  for i in range(Msize): 
    # print out y2=y0+h*(a2_0*dy0+a2_1*dy1) when i > 0 
    if i > 0:
      addLines += BreakLines(['a%d_%d*dy%d' % (i,j,j) for j in range(i)],'+','y%d=y0+h*('%i,' + &',')',breakNo,spaces)
    addLines.append('%s! use y%d to get dy%d' % (spaces, i, i) + '\n') 
    addLines.append(spaces + 'CALL '+ devString.replace('0', '%d' % i) + '\n')
    addLines.append('\n') 

  # addLines.append('delta=h*(%s)' % ('+'.join('(b%d-b_%d)*dy%d' % (j, j, j) for j in range(Msize))) + '\n') 
  addLines += BreakLines(['b%d*dy%d' % (j,j) for j in range(Msize)],'+','yn=y0+h*(',' + &',')',breakNo+2,spaces)
  if yerr:
    addLines += BreakLines(['b_%d*dy%d' % (j,j) for j in range(Msize)],'+','ynp=y0+h*(',' + &',')',breakNo+2,spaces)
    addLines += BreakLines(['(b%d-b_%d)*dy%d' % (j, j, j) for j in range(Msize)],'+','yerr=h*(',' + &',')',int((breakNo+2)/1.5),spaces) 
  #    y5(:)=y0(:)+a51*h*dy1(:)+a52*h*dy2(:)+a53*h*dy3(:)+a54*h*dy4(:)+a55*h*dy5(:)
  #    CALL Getdy(y5,dy6,ODEMode,rerun,Energy,Jtot)
  return addLines
    


if __name__ == '__main__':
  # order = 8
  Msize = 13 
  devString = 'dev(y0,dy0,h,hnew)'
  breakNo = 6 # ...dy6 + & 
  addLines = []
  #f = open('RKGeneralTemp.dat', 'w')
  addLines += PrintOutdys(Msize=13, devString = 'dev(y0,dy0,test)', breakNo = 6)
  #PrintOutdys(f, Msize=13, devString = 'dev(y0,dy0,h,hnew)', breakNo = 6)
  #f.close()
  #with open('RKGeneralTemp.dat', 'r') as f:
  #  lines = f.readlines()
  #  for line in lines:
  #    print(line[:-1]) 
  addLines += Printdy(13, breakNo-1, Head=' REAL(KIND=8), DIMENSION(SIZE(y0)) :: ')
  addLines += Printdy(13, breakNo-1, Head=' REAL(KIND=8), DIMENSION(SIZE(y0)) :: ', mode='dy')

  print(addLines)


