
CopyRightInfo = '\
#  A family of embedded Runge-Kutta formulae, by J. R. Dormand and P. J. Prince, \n\
# Journal of Computational and Applied Mathematics, Vol 6, No. 1, 1980, pages 19 to 26.\n\
# downloaded from http://www.peterstone.name/Maplepgs/Maple/nmthds/RKcoeff/Runge_Kutta_schemes/RK5/RKcoeff5d_1.pdf'

a = '\
c[1]=0|\
c[2]=1/5|\
c[3]=3/10|\
c[4]=4/5|\
c[5]=8/9|\
c[6]=1|\
c[7]=1|\
a[2,1]=1/5|\
a[3,1]=3/40|\
a[3,2]=9/40|\
a[4,1]=44/45|\
a[4,2]=-56/15|\
a[4,3]=32/9|\
a[5,1]=19372/6561|\
a[5,2]=-25360/2187|\
a[5,3]=64448/6561|\
a[5,4]=-212/729|\
a[6,1]=9017/3168|\
a[6,2]=-355/33|\
a[6,3]=46732/5247|\
a[6,4]=49/176|\
a[6,5]=-5103/18656|\
a[7,1]=35/384|\
a[7,2]=0|\
a[7,3]=500/1113|\
a[7,4]=125/192|\
a[7,5]=-2187/6784|\
a[7,6]=11/84|\
b[1]=35/384|\
b[2]=0|\
b[3]=500/1113|\
b[4]=125/192|\
b[5]=-2187/6784|\
b[6]=11/84|\
b[7]=0|\
b*[1]=5179/57600|\
b*[2]=0|\
b*[3]=7571/16695|\
b*[4]=393/640|\
b*[5]=-92097/339200|\
b*[6]=187/2100|\
b*[7]=1/40|\
'

print(CopyRightInfo) 
print('\n'*2)
from PeterStone import * 
PrintFixed(a)
print('\n'*2)

