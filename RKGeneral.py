# output *f90 file with given order 
import numpy as np
import sys
# from RKGeneral import Printdy, PrintOutdys 


def Printdy(jsize, breakNo, OpenFile, Head=' REAL(KIND=16), DIMENSION(SIZE(y0)) :: ', mode='y'):
  # print(Head, file=OpenFile) # REAL(KIND=DP), DIMENSION(SIZE(y0)) :: 
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
    print(Head + ','.join('%s%d' % (mode,j+Add) for j in range(left)), file=OpenFile) 
    return 
  while left > 0: # more lines to write 
    if start > 0:
      Head = ' '* len(Head)
    End = ', &'
    if left - breakNo <= 0:
      End = ''
    line = ','.join('%s%d' % (mode,j+Add) for j in range(start, min(start+breakNo,maxSize)))      
    print(Head + line + End, file=OpenFile)
    left = left - breakNo
    start += breakNo
  return 

def BreakLines(List, joinSign, FirstLineHead, LinkSign, EndString, breakNo, OpenFile, spaces):
  jsize = np.size(List)  
  addSingleLine = '%s%s%s' % (FirstLineHead, joinSign.join(List[i] for i in range(jsize)), EndString) 
  # print('addSingleLine', addSingleLine) 
  left = jsize
  start = 0 
  printed = False
  if left < breakNo:
    print(spaces+addSingleLine, file=OpenFile) 
    printed = True
  addBreak = ''
  while left > 0: 
    line = '%s%s' % (FirstLineHead, joinSign.join(List[i] for i in range(start, min(start+breakNo,jsize))))   
    if left - breakNo <= 0:
      if not printed:
        print(spaces+line + EndString, file=OpenFile) 
      addBreak = addBreak + line + EndString
    else:
      print(spaces+line + LinkSign, file=OpenFile) 
      addBreak = addBreak + line + LinkSign 
    FirstLineHead = ' '*(len(FirstLineHead)-3) + '&  ' 
    left = left - breakNo 
    start += breakNo
  # check that these two matches !!! 
  addBreakClean = addBreak.replace(' ','').replace('&','')
  if not addBreakClean == addSingleLine:
    print('addSingleLine = [%s]' % addSingleLine, file=OpenFile) 
    print('addBreakClean = [%s]' % addBreakClean, file=OpenFile) 
    raise ValueError('add break does not match add') 
  return 


def PrintOutdys(OpenFile, Msize=13, devString = 'dev(y0,dy0,test)', breakNo = 6, spaceNo=2, test=False, yerr=True): 
  # if (test), return with smaller h
  spaces = ' '*spaceNo 
  for i in range(Msize): 
    # print out y2=y0+h*(a2_0*dy0+a2_1*dy1) when i > 0 
    if i > 0:
      BreakLines(['a%d_%d*dy%d' % (i,j,j) for j in range(i)],'+','y%d=y0+h*('%i,' + &',')',breakNo,OpenFile,spaces)
    print('%s! use y%d to get dy%d' % (spaces, i, i), file=OpenFile) 
    print(spaces + 'CALL ', devString.replace('0', '%d' % i), file=OpenFile)
    if test:
      testString = '\
%sIF (test) THEN ! bad dy \n\
%s  IF (ABS(h-MinStepSize)/MinStepSize.LE.1Q-13) RETURN ! stop the program \n\
%s  hnew = MAX(MIN(h/2.Q0,MaxStepSize),MinStepSize) ! reduce step \n\
%s  rerun = .True.\n\
%s  test = .False.\n\
%s  RETURN\n\
%sEND IF' % (spaces, spaces, spaces, spaces, spaces, spaces, spaces) 
      print(testString, file=OpenFile) 

    print('', file=OpenFile) 
  # print('delta=h*(%s)' % ('+'.join('(b%d-b_%d)*dy%d' % (j, j, j) for j in range(Msize))), file=OpenFile) 
  BreakLines(['b%d*dy%d' % (j,j) for j in range(Msize)],'+','yn=y0+h*(',' + &',')',breakNo+2,OpenFile,spaces)
  if yerr:
    BreakLines(['b_%d*dy%d' % (j,j) for j in range(Msize)],'+','ynp=y0+h*(',' + &',')',breakNo+2,OpenFile,spaces)
    BreakLines(['(b%d-b_%d)*dy%d' % (j, j, j) for j in range(Msize)],'+','yerr=h*(',' + &',')',int((breakNo+2)/1.5),OpenFile,spaces) 
  #    y5(:)=y0(:)+a51*h*dy1(:)+a52*h*dy2(:)+a53*h*dy3(:)+a54*h*dy4(:)+a55*h*dy5(:)
  #    CALL Getdy(y5,dy6,ODEMode,rerun,Energy,Jtot)
  return 
    


if __name__ == '__main__':
  # order = 8
  Msize = 13 
  devString = 'dev(y0,dy0,h,hnew)'
  breakNo = 6 # ...dy6 + & 
  #f = open('RKGeneralTemp.dat', 'w')
  PrintOutdys(sys.stdout, Msize=13, devString = 'dev(y0,dy0,test)', breakNo = 6)
  #PrintOutdys(f, Msize=13, devString = 'dev(y0,dy0,h,hnew)', breakNo = 6)
  #f.close()
  #with open('RKGeneralTemp.dat', 'r') as f:
  #  lines = f.readlines()
  #  for line in lines:
  #    print(line[:-1]) 
  Printdy(13, breakNo-1, sys.stdout, Head=' REAL(KIND=16), DIMENSION(SIZE(y0)) :: ')
  Printdy(13, breakNo-1, sys.stdout, Head=' REAL(KIND=16), DIMENSION(SIZE(y0)) :: ', mode='dy')




