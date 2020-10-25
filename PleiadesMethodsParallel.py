from odefortran import * 
from Common import *
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d
 
FileName = sys.argv[0].replace('.py','').replace('Parallel','')
datafilename = sys.argv[1] 

filename = datafilename.replace('npz', 'dat') 

method = filename.split('_')[1] 
rtol = float(filename.split('_')[2].replace('rtol',''))
atol = float(filename.split('_')[3].replace('atol','').replace('.dat',''))
# Data/PleiadesMethods_rk108Feagin_rtol1E-32_atol1E-32.npz
print('method = [%s]' % (method)) 
print('rtol = [%e]' % (rtol)) 
print('atol = [%e]' % (atol)) 

# inital and final values 
y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
y0 = np.array([mpf(str(yEach)) for yEach in y0])

methods = ['rk54Sharp','rk54Dormand','rk1412Feagin','rk108Feagin','rk1210Feagin','rk87EnrightVerner','rk65Dormand','rk1412Long','rk1211Peter','rk109Legendre','rk87Dormand']

t, y, outDict = solve_ivp_(mpf('0.'),mpf('3.'),y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=mpf('1E-3'),max_step=mpf('1.'),min_step=mpf('1E-10'),Print6=True,load=True)
if outDict['Success']:
  print('Success') 
