import numpy as np
import ode
import sys, os
from Python.Common import rmrf
from mpmath import mpf, nstr, nprint, mp
import mpmath
mp.dps = 50 # 50 digits  
 
IndexList = []
for Name in ['x', 'y', 'dx', 'dy']:
  for i in range(7):
    IndexList.append('%s%d' % (Name, i+1))


def solve_ivp_(t0,tfinal,y0,method='rk1412Feagin',filename='',atol=1E-6,rtol=1E-3,first_step=1E-3,max_step=np.inf, min_step=0.,Print6=True,rm=False):
  """'solve_ivp_' wrapped the fortran code ODEInterface.f90 (=>ode.so) 
  

  Usage: 
  t0, tfinal = 0., 3.
  y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
  t, y, msgDict, filename = solve_ivp_(t0, tfinal, y0,method='rk1412Feagin',filename='Data/test.dat', atol=1E-6, rtol=1E-7, first_step=1E-2, max_step=1., min_step=1E-4, Print6=True) 

  t0         : initial time
  tfinal     : final time 
  y0         : initial condition
  method     : ode solver, should be in ['rk54Sharp','rk54Dormand','rk65Dormand','rk87Dormand','rk87EnrightVerner','rk108Feagin','rk109Legendre','rk1210Feagin','rk1211Peter','rk1412Feagin']
  filename   : format:'Data/*.dat', filename='' if you do not want to save it
  atol       : input a float or array of np.size(y0), default value 1E-6
  rtol       : float, default value 1E-3
  first_step : float, default value 1E-3
  max_step   : float, default value np.inf 
  min_step   : float, default value 0.
  Print6     : logical, True if you want to see output on screen, default True  
  rm         : rm the Data/*.npz even if it already exist, default False  

  return
  t0         : array of time
  y          : array of data, IndexList=['x1','x2',....]
  msg        : a dictionary that contain all the output 
  filename   : Data/*.npz

  more 
  msg = {
  'cpuTime': float # how much cpuTime is used in fortran code
  'Minh'   : float # the min(h) used in the fortran code, h = time_step
  'Maxh'   : float # the max(h) used in the fortran code, h = time_step
  'Rejected'  : int # how many times the solver reject a time step 
  'Accepted'  : int # how many times the solver accept a time step without being rejected before 
  'Evaluated' : int # =fcn, how many times the y'=f(t,y) is Evaluated
  'TotalSteps': int # the total number of time step (np.size(t) \\approx min(TotalSteps,5.123124E5)   
  'ReachMax'  : int # how many times h reaches max_step, if this number is large, you may want to change max_step 
  'ReachMin'  : int # how many times h reaches min_step, if this number is large, you may want to change min_step 
  'Success'   : 0 or 1 # 1=succeeded, 0=the program did not reach the end 
  }

  """
  if rm: # force removal 
    rmrf(filename.replace('.dat', '.npz'))
    rmrf(filename) 
  if os.path.isfile(filename.replace('.dat', '.npz')): # load data, do not calculate
    data = np.load(filename.replace('.dat', '.npz'),allow_pickle=True)
    try:
      return data['t'], data['y'], data['msg'].item(), filename.replace('.dat', '.npz')
    except KeyError:
      print('The file is probably old, removing it anyway')  
      rmrf(filename.replace('.dat', '.npz')) 
      rmrf(filename) 
  if np.size(atol) == 1:
    atol = np.full(np.size(y0), atol) 
  IntArray, RealArray = ode.odeinterfacemod.ode(y0,t0,tfinal,method,filename,atol,rtol,max_step,min_step,first_step,Print6)
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
    data = np.loadtxt(filename, dtype=np.float64)
    t = data[:, 0]
    y = data[:, 1:] 
    # print(np.shape(data), np.shape(data[0]), np.shape(data[:, 0]))   
    # print(np.shape(t), np.shape(y)) 
    np.savez(filename.replace('.dat', '.npz'), t=t, y=y, msg=Dict)  
    rmrf(filename) 
    return t, y, Dict, filename.replace('.dat', '.npz') 

if __name__ == '__main__':
  FileName = sys.argv[0].replace('.py','')
  # The Pleiades problem  
  # from https://archimede.dm.uniba.it/~testset/report/plei.pdf 
  # ref https://archimede.dm.uniba.it/~testset/problems/plei.php
  y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
  yfinal = np.load('DataSaved/PleiadesSolution_QP_rk1412Feagin_yfinal.npz',allow_pickle=True)['yfinal']
  t, y, outputDict, filename = solve_ivp_(0., 3., y0, filename='Data/test.dat', atol=1E-8, rtol=1E-8, first_step=1E-2, max_step=1., min_step=1E-4, Print6=True, rm=True) 
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
  ax[0].set_xlabel(r'$t\,[{\rm s}]$') 
  ax[1].set_xlabel(r'$t\,[{\rm s}]$') 
  ax[0].set_ylabel(r'$x_1$') 
  ax[1].set_ylabel(r'$y_1$') 

  ax[2].plot(y[:, IndexList.index('x1')], y[:, IndexList.index('y1')], 'rx', label='p1')
  ax[2].plot(y[:, IndexList.index('x3')], y[:, IndexList.index('y3')], 'bo', label='p3')
  ax[2].set_xlabel(r'$x$')
  ax[2].set_ylabel(r'$y$')
  ax[2].set_xlim(-1, 3)
  ax[2].set_ylim(0, 4)
  ax[2].legend(loc='best')

  # plt.show()
  fig.tight_layout()
  fig.savefig('Plots/%s_xy.png' % (FileName), dpi=300)
  print('Plots/%s_xy.png' % (FileName)) 
  fig.clf()
