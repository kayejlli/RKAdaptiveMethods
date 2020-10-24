import numpy as np
import ode
import sys, os
from mpmath import mpf, nstr, nprint, mp
import mpmath
from Common import rmrf
mp.dps = 100  
# from odefortran import * 
# solve_ivp_(0., 3., y0, filename='Data/test.dat', atol=1E-6, rtol=1E-7, first_step=1E-2, max_step=1., min_step=1E-4, Print6=True) 
 
IndexList = []
for Name in ['x', 'y', 'dx', 'dy']:
  for i in range(7):
    IndexList.append('%s%d' % (Name, i+1))


LocalPrint = lambda t0: nstr(t0,60,show_zero_exponent=True,min_fixed=-2, strip_zeros=True).replace('e','Q') 

def WriteQP(f1,y0,t0,tfinal,max_step,min_step,first_step):
  for i in range(np.size(y0)):
    print('  y0(%d) = %s' % (i+1, LocalPrint(y0[i])), file=f1) 
  Dict = {'t0': t0, 'tfinal': tfinal, 'max_h': max_step, 'min_h': min_step,\
          'hinit': first_step} 
  for name in Dict:  
    print('  %s = %s' % (name, LocalPrint(Dict[name])), file=f1) 
  return 


def solve_ivp_(t0,tfinal,y0,method='rk1412Feagin',filename='',atol=mpf('1E-6'),rtol=mpf('1E-3'),first_step=mpf('1E-3'),max_step=np.inf, min_step=mpf('0.'),Print6=True,load=True):
  """'solve_ivp_' wrapped the fortran code 

  Usage: 
  y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
  solve_ivp_(0., 3., y0,method='rk1412Feagin',filename='Data/test.dat', atol=1E-6, rtol=1E-7, first_step=1E-2, max_step=1., min_step=1E-4, Print6=True) 

  t0         : initial time
  tfinal     : final time 
  y0         : initial condition
  method     : ode solver, should be in ['rk1412Feagin'] 
  filename   : data file, filename='' if you do not want to save it
  atol       : input a float or array of np.size(y0), default value 1E-6
  rtol       : float, default value 1E-3
  first_step : float, default value 1E-3
  max_step   : float, default value np.inf 
  min_step   : float, default value 0.
  Print6     : logical, True if you want to see output on screen, default True  

  """
  # fix atol
  if np.size(atol) == 1:
    atol = np.full(np.size(y0), atol) 
  if load and os.path.isfile(filename.replace('dat', 'npz')): # load data, do not calculate 
    data = np.load(filename.replace('dat', 'npz'),allow_pickle=True)
    return data['t'], data['y'], data['msg'].item() 

  '''
  # generate fortron file ODEInterfaceQP.f90 
  SampleFile = 'ODEInterfaceQPSample.f90' # include here
  f1 = open(SampleFile.replace('Sample', ''), 'w') 
  with open(SampleFile) as f0:
    lines = f0.readlines()
    for line in lines:
      if not '-----' in line:
        f1.write(line) # copy the line to f1 
      else:
        WriteQP(f1,y0,t0,tfinal,max_step,min_step,first_step) 
  f1.close()
  # recompile the program 
  bashCommand = 'make odeqp.so'
  import subprocess
  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
  output, error = process.communicate()
  #print('output', output) 
  #print('error', error) 
  '''

  # import 
  import odeqp 
  IntArray, RealArray = odeqp.odeinterfaceqpmod.odeqp(method,filename,Print6,rtol,atol)  

  RealArrayIndexes = ['cpuTime', 'Minh', 'Maxh'] 
  IntArrayIndexes = ['Rejected','Accepted','Evaluated','TotalSteps','ReachMax','ReachMin', 'Success'] 
  Dict = {}
  for i, name in enumerate(RealArrayIndexes):
    Dict[name] = RealArray[i]
  for i, name in enumerate(IntArrayIndexes):
    Dict[name] = IntArray[i]
  if Print6:
    print(' CPU_Time=%E' % (RealArray[0])) 
  if filename == '':
    return 
  else:
    data = np.loadtxt(filename, dtype=mpmath.ctx_mp_python.mpf)
    t = data[:, 0]
    y = data[:, 1:] 
    print(np.shape(data), np.shape(data[0]), np.shape(data[:, 0]))   
    print(np.shape(t), np.shape(y)) 
    t = np.array([mpf(tt) for tt in t]) 
    ycopy = np.full(np.shape(y), mpf('1.')) 
    for i in range(np.shape(y)[0]):
      for j in range(np.shape(y)[1]):
        ycopy[i, j] = mpf(y[i,j])  
    np.savez(filename.replace('dat', 'npz'), t=t, y=ycopy, msg=Dict)  
    rmrf(filename) 
    if load:
      return t, ycopy, Dict
    else:
      return filename.replace('dat', 'npz'), Dict  

if __name__ == '__main__':
  FileName = sys.argv[0].replace('.py','')
  # The Pleiades problem  
  # from https://archimede.dm.uniba.it/~testset/report/plei.pdf 
  # ref https://archimede.dm.uniba.it/~testset/problems/plei.php
  y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
  y0 = np.array([mpf(str(yEach)) for yEach in y0]) 
  yfinal = [0.3706139143970502,0.3237284092057233E1,-0.3222559032418324E1,0.6597091455775310,\
  0.3425581707156584,0.1562172101400631E1,-0.7003092922212495,-0.3943437585517392E1,\
  -0.3271380973972550E1,0.5225081843456543E1,-0.2590612434977470E1,0.1198213693392275E1,\
  -0.2429682344935824,0.1091449240428980E1,0.3417003806314313E1,0.1354584501625501E1,\
  -0.2590065597810775E1,0.2025053734714242E1,-0.1155815100160448E1,-0.8072988170223021,\
  0.5952396354208710,-0.3741244961234010E1,0.3773459685750630,0.9386858869551073,\
  0.3667922227200571,-0.3474046353808490,0.2344915448180937E1,-0.1947020434263292E1]
  yfinal = [mpf('%20.15e' % yEach) for yEach in yfinal]
  rmrf('Data/Reference.npz') 
  t, y, outputDict = solve_ivp_(mpf('0.'), mpf('3.'), y0, filename='Data/Reference.dat', atol=mpf('1E-15'), rtol=mpf('1E-15'), first_step=mpf('1E-3'), max_step=mpf('1.'), min_step=mpf('1E-10'), Print6=True, load=True) 
  # data = np.load(filename) 
  # t = data['t']
  # y = data['y']
  import pylab as plt
  import matplotlib.gridspec as gridspec

  fig = plt.figure(figsize=(10, 8)) 
  gs = gridspec.GridSpec(2, 2)
  ax = []
  for i in range(2):
    ax.append(fig.add_subplot(gs[0, i])) 
  ax.append(fig.add_subplot(gs[1, :])) 
  ax[0].plot(t, y[:, IndexList.index('x1')], 'k-')
  ax[1].plot(t, y[:, IndexList.index('y1')], 'k-')
  ax[0].set_xlim(0., 3.)
  ax[1].set_xlim(0., 3.)
  ax[0].set_ylim(-1, 3)
  ax[1].set_ylim(-4, 4)
  ax[2].plot(y[:, IndexList.index('x1')], y[:, IndexList.index('y1')], 'rx', label='p1')
  ax[2].plot(y[:, IndexList.index('x3')], y[:, IndexList.index('y3')], 'bo', label='p3')
  ax[2].set_xlim(-1, 3)
  ax[2].set_ylim(0, 4)
  ax[2].legend(loc='best')
  plt.show()
  #fig.savefig('Plots/%s_xy.png' % (FileName), dpi=300)
  #print('Plots/%s_xy.png' % (FileName)) 
  #fig.clf()





