# graphic setting
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
# matplotlib.style.use('default')

import pylab as plt
from matplotlib import rc
import numpy as np
#from matplotlib import ticker
#from mpl_toolkits.mplot3d import Axes3D

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)

params = {'legend.fontsize': 15,'ytick.labelsize': 18, 
          'axes.labelsize' : 18, 'xtick.labelsize': 18,
          'lines.linewidth': 1, 'font.size': 18}
plt.rcParams.update(params)
plt.rcParams['agg.path.chunksize'] = 10000

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Calibri'


NameColorList = ['purple', 'green', 'blue', 'pink', 'brown', 'red', 'sky blue', 'teal', 'orange', 'magenta', 'grey', 'tan', 'salmon',\
                 'royal blue', 'spring green', 'coral', 'blue violet', 'dark cyan']

Color3 = ['b', 'r', 'g']
Color4 = ['b', 'r', 'g', 'm']
Color5 = ['b', 'r', 'g', 'm', 'orange']
Color6 = ['b', 'r', 'g', 'm', 'orange', 'blueviolet'] 
ls4 = ['-', '--', '-.', ':'] 


def SaveAndClose(fig,name,ax,xlabel='',ylabel='',zlabel='',legend=False,loc='best'):
  if zlabel == '':
    if np.size(xlabel) == 1 and not xlabel == '':
      if np.size(ax) == 1:
        ax.set_xlabel(xlabel)
      else:
        for i, ax1 in enumerate(ax):
          ax1.set_xlabel(xlabel)
    elif np.size(xlabel) > 1:
      for i, ax1 in enumerate(ax):
        if not xlabel[i] == '':
          ax1.set_xlabel(xlabel[i])
    if np.size(ylabel) == 1 and not ylabel == '':
      if np.size(ax) == 1:
        ax.set_ylabel(ylabel)
      else:
        for i, ax1 in enumerate(ax):
          ax1.set_ylabel(ylabel)
    elif np.size(ylabel) > 1:
      for i, ax1 in enumerate(ax):
        if not ylabel[i] == '':
          ax1.set_ylabel(ylabel[i])
       
    if np.size(legend) > 1:
      for i, legendValue in enumerate(legend):
        if legendValue:
          ax[i].legend(loc=loc)
    else:
      if legend:
        try:
          for ax1 in ax:
            ax1.legend(loc=loc)
        except TypeError:
          ax.legend(loc=loc)
  else:
    try:
      ax.set_xlabel(xlabel)
      ax.set_ylabel(ylabel)
    except AttributeError:
      ax[0].set_xlabel(xlabel)
      ax[0].set_ylabel(ylabel)

  fig.tight_layout()
  fig.savefig('%s.png' % (name), dpi=300)
  # fig.savefig('%s.pdf' % (name), dpi=300)
  fig.clf()
  return




