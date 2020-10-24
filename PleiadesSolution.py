from odefortran import * 
from Common import *
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d
 
FileName = sys.argv[0].replace('.py','')

def mynsum(Array):
  Sum = mpf('0.')
  for No in Array:
    Sum = Sum + No
  return Sum

# inital and final values 
y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
y0 = np.array([mpf(str(yEach)) for yEach in y0])


tolArray = [1E-50, 1E-40, 1E-38, 1E-35, 1E-34, 1E-32, 1E-30, 1E-25, 1E-20, 1E-16, 1E-15, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]
tolArray = np.array(tolArray) 

method = 'rk1412Feagin'
for i in range(1, 8):
  globals()['x%d' % i] = np.full(np.size(tolArray), mpf(0))
  globals()['y%d' % i] = np.full(np.size(tolArray), mpf(0))
  globals()['dx%d' % i] = np.full(np.size(tolArray), mpf(0))
  globals()['dy%d' % i] = np.full(np.size(tolArray), mpf(0))

for j, atol in enumerate(tolArray):
  rtol = atol 
  filename = 'Data/%s_%s_rtol%s_atol%s.dat' % (FileName, method, PrintExponentialBetter(rtol), PrintExponentialBetter(atol)) 
  t, y, outDict = solve_ivp_(mpf('0.'),mpf('3.'),y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=mpf('1E-3'),max_step=mpf('1.'),min_step=mpf('1E-10'),Print6=False,load=True)
  if outDict['Success']:
    for i, name in enumerate(IndexList): 
      globals()[name][j] = y[-1][i] 
  else:
    t, y, outDict = solve_ivp_(mpf('0.'),mpf('3.'),y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=mpf('1E-3'),max_step=mpf('1.'),min_step=mpf('1E-10'),Print6=True,load=True)

yfinal = np.full(28, mpf(1))   
fig0, ax0 = plt.subplots()
for i in range(1,8):
  fig, ax = plt.subplots() 
  for j, atol in enumerate(tolArray):
    if j < 2 or atol > 1E-10:
      continue 
    arg = tolArray < atol
    xArray = globals()['x%d' % i].copy()
    yArray = globals()['y%d' % i].copy()
    xArrayMean = mynsum(xArray[arg])/np.size(tolArray[arg]) 
    yArrayMean = mynsum(yArray[arg])/np.size(tolArray[arg]) 
    line, = ax.loglog(tolArray[arg], abs(xArray[arg]-xArrayMean), marker='o', label=r'$atol <  %s$' % (LabelExponentialBetter(atol))) 
    color = line.get_color()
    ax.loglog(tolArray[arg], abs(yArray[arg]-yArrayMean), color=color, marker='v') #, label=r'$y_1$')  
    if j == 5:
      arg[0] = False 
      dyArray = globals()['dy%d' % i].copy()
      dxArray = globals()['dx%d' % i].copy()
      xArrayMean = mynsum(xArray[arg])/np.size(tolArray[arg]) 
      dxArrayMean = mynsum(dxArray[arg])/np.size(tolArray[arg]) 
      yArrayMean = mynsum(yArray[arg])/np.size(tolArray[arg]) 
      dyArrayMean = mynsum(dyArray[arg])/np.size(tolArray[arg]) 
      yfinal[i-1] = xArrayMean
      yfinal[i-1+7] = yArrayMean
      yfinal[i-1+14] = dxArrayMean
      yfinal[i-1+21] = dyArrayMean
      line, = ax0.semilogx(tolArray[arg], xArray[arg]-xArrayMean, label=r'$x_%d$' % (i), marker='o')  
      color = line.get_color()
      ax0.semilogx(tolArray[arg], yArray[arg]-yArrayMean, marker='v') 
  ax.legend(loc='best',fontsize=10.)
  ax.set_xlabel(r'atol')
  ax.set_ylabel(r'$\delta x_%d, \delta y_%d$' % (i,i)) 
  fig.tight_layout()
  fig.savefig('Plots/%s_error_particle%d.png' % (FileName, i), dpi=300)
  fig.clf()

ax0.legend(loc='best',fontsize=10.)
ax0.set_xlabel(r'atol')
ax0.set_ylabel(r'$\delta x, \delta y$') 
fig0.tight_layout()
fig0.savefig('Plots/%s_error.png' % (FileName), dpi=300)
fig0.clf()
np.savez('Data/%s_QP_rk1412Feagin_yfinal.npz' % (FileName), yfinal=yfinal) 

