import numpy as np
import sys, os
from mpmath import mpf, nstr, nprint, mp
import mpmath
mp.dps = 50 # 50 digits  
 
IndexList = []
for Name in ['x', 'y', 'dx', 'dy']:
  for i in range(7):
    IndexList.append('%s%d' % (Name, i+1))

# The Pleiades problem  
# from https://archimede.dm.uniba.it/~testset/report/plei.pdf 
# ref https://archimede.dm.uniba.it/~testset/problems/plei.php

if __name__ == '__main__':
  FigFileName = sys.argv[0].replace('.py','')
  filename = 'Data/test.dat'
  data = np.loadtxt(filename, dtype=np.float64)
  t = data[:, 0]
  y = data[:, 1:] 

  # the accurate final results
  yfinal = np.load('DataSaved/PleiadesSolution_QP_rk1412Feagin_yfinal.npz',allow_pickle=True)['yfinal']

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
  fig.savefig('Plots/%s_xy.png' % (FigFileName), dpi=300)
  print('Plots/%s_xy.png' % (FigFileName)) 
  fig.clf()

  print('----- The error in the y[t=tfinal] -----')
  for i, var in enumerate(IndexList):
    print('%4s=%23.15e -- err=%.4e' % (var, y[-1, i], abs(y[-1,i]-yfinal[i]))) 
 
