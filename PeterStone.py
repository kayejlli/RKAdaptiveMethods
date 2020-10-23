# http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK10/RKcoeff10r_1.pdf
# from PeterStone import * 
# PrintFixed(a)
from Common import GetFortranFloat

def PrintFixed(a,Mode='fortran',FortranBreak=False,HeadSpace=2):
  # print(a.split('|')[0])
  if FortranBreak:
    Head = ' '* HeadSpace
    End = ' ,&' 
  else:
    Head, End = '', ''
  for name in a.split('|'):
    if '=' in name:
      first = name.split('=')[0]
      second = name.split('=')[1] 
      
      #if '/' in second: # in fraction 
      second = GetFortranFloat(second, fortran=False, Precision=85)
      if not ',' in first:
        no = int(first.split('[')[1][:-1])
        ab = first.split('[')[0] 
        if ab == 'c':
          ab = 'k'
        ab = ab.replace('*', '_') 
        if Mode == 'bracket':
          print('%s[%d]=%s' % (ab, no-1, second))
        elif Mode == 'fortran':
          print('%s%s%d=%s%s' % (Head,ab, no-1, second,End))
      else:
        no0 = int(first.split('[')[1][:-1].split(',')[0]) 
        no1 = int(first.split('[')[1][:-1].split(',')[1]) 
        ab = first.split('[')[0] 
        if Mode == 'bracket':
          print('%s[%d,%d]=%s' % (ab,no0-1,no1-1,second)) 
        elif Mode == 'fortran':
          print('%s%s%d_%d=%s%s' % (Head,ab,no0-1,no1-1,second,End)) 
       
 
