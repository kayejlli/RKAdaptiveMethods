
CopyRightInfo = '\
# modification of the Prince-Dormand 6 stage\n\
# downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5c_2.pdf'

a = '\
c[1]=0|\
'


print(CopyRightInfo) 
print('\n'*2)
from PeterStone import * 
PrintFixed(a)
print('\n'*2)
