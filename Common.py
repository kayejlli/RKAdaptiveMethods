# graphic setting
from mpmath import mp, mpf, nstr
mp.dps = 100 
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import sys, os


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
Marker5 = ['o', 'v', '+', 'v', 'x'] 


def rmrf(filename):
  try:
    os.remove(filename)
  except OSError:
    pass
  return 



def PrintExponentialBetter(number):
  try:
    # clean up number greater than 10
    string0 = '%.2E' % (number)
    exp = int(string0.split('E')[-1].replace('+',''))
    num = float(string0.split('E')[0])
    dif = abs(num-round(num))
    while dif > 1E-3:
      num = num*10.
      exp = exp - 1
      dif = abs(num-round(num))
    if exp == 0 or exp == 1:
      return '%d' % (num*10**exp)
    else:
      return '%dE%d' % (num, exp)
  except TypeError:
    return str(number)


def LabelExponentialBetter(number):
  label0 = PrintExponentialBetter(number)
  if not 'E' in label0:
    return label0
  elif label0.split('E')[0] == "1": # 1E-3 
    return r'10^{%s}' % (label0.split('E')[1]) 
  elif len(label0.split('E')[0]) > 1: # 123E1 
    part0 = label0.split('E')[0]
    part1 = int(label0.split('E')[1]) + (len(part0) - 1) 
    if part1 == 0:
      return r'%s.%s' % (part0[0], part0[1:]) 
    else:
      return r'%s.%s \times 10^{%d}' % (part0[0], part0[1:], part1)
  else: # 5E-4 
    return r'%s \times 10^{%s}' % (label0.split('E')[0], label0.split('E')[1]) 


def addReference(ax, ref, mode='x', ls='k--', alpha=0.7):
  if mode == 'x':
    xmin, xmax = ax.get_xlim()
    ax.plot([xmin, xmax], [ref, ref], ls, alpha=alpha)
    ax.set_xlim(xmin, xmax)
  elif mode == 'y':
    ymin, ymax = ax.get_ylim()
    ax.plot([ref, ref], [ymin, ymax], ls, alpha=alpha) 
    ax.set_ylim(ymin, ymax)
  return 


def TrimRound(String, Length=60):
  # if the String == '9'*len(String), then error will appears 
  cut = 1
  # print('String = [%s]' % (String)) 
  # print('len=%d, Pre=%d, cut=%d' % (len(String), Length, cut)) 
  while String[Length-cut] == '9':
    cut += 1
  #print('String=[%s]' % ( String))  
  left = String[:Length-cut]
  #print('left=[%s]' % ( left))  
  key = int(String[Length-cut: Length])    
  right = String[Length]
  #print('right=[%s]' % ( right))  
  if int(right) >= 5:
    return '%s%d' % (left, key+1)
  else:
    return '%s%d' % (left, key)


def mpfToString(valuempf,fortran=True, Precision=60):
  valuestring = nstr(valuempf, Precision, show_zero_exponent=True, min_fixed=-2, strip_zeros=True) 
  if not 'e' in valuestring:
    raise ValueError(valuestring) 
  if fortran:
    return valuestring.replace('e+', 'D').replace('e','D')
  else:
    return valuestring.replace('e+', 'E').replace('e','E') 


def GetFortranFloat(String, Precision=60, fortran=True):
  String = String.replace('e', 'E').strip() 
  if len(String) > 500:
    raise ValueError('The string is too long! len=%d' % (len(String))) 
  # check if it is a fraction 
  if '/' in String:
    numerator, denominator = String.split('/')
    valuempf = mpf(numerator)/mpf(denominator)  
    return mpfToString(valuempf, Precision=Precision, fortran=fortran) 
  else:
    if len(String) > mp.dps:
      raise ValueError('Your precision choice is %d, but your string [%d] long = [%s]' % (mp.dps, len(String), String)) 
    valuempf = mpf(String) # convert to multiple precision 
    return mpfToString(valuempf, Precision=Precision, fortran=fortran) 

'''
  elif 'E' in String:
    if not String == '0E0':
      return String.replace('E', 'D')
    else:
      return '0.D0' 
  elif '.' in String:
    return String + 'D0' 
  else: # integer 
    return String + '.D0' 
'''

if __name__ == '__main__':
  bList = '\
a8_4=-991136638972626168678903371416456100093900405535164924683058122802429707354033382826947398158683765324439618282500000000/79848142008846002789298925227605775190331269194726743910364273272231784282184770467794155269096224513726772081370189773|\
a8_5=99411279821210413387149352497211785829547149358696952646781033905129048593757052549024957474512389085654050445280000000/10219750071154355529199360565479179548497227516418867455143546824991883616697162903707552409044338358355837306818391433|\
a8_6=194327599672380134095898291719912961363678073793023525007081328425098431574448809779310732532821200046895000000000/11996011488227722649656673931136490891256463292620473601841875115170259043531987881275965965525693289412424400177|\
a8_7=-4738143867684122189593816244199450540483384372163549951990525387550768038015218275414120082248510000000000/45627381556119209828916528385326434273376137158228892503158792567081801761628344925618992749059885819540261|\
a9_0=100509763879264306824096153463041174636629364248095333923106653001873/229490324007644628042361756217436155832461488260089524115475000000000|\
a9_3=35076261889213578261995286390053983221937920015616/8903662403052468890234321680409039895089390971875|\
a9_4=29877053472248545227782869189767925950557009/10443709158645362958089905740134206606110000|\
a9_5=-72602025182798889442893966553844286012776770019588838776297031451/40918439380405007926071673834718276106178863851548244800289187500|\
a9_6=-15322912063864370512130145988492098605486502994107816190/3130069400645604876183669549605539076333579290524919889|\
a9_7=66085154677219418645471125072555541174985695924/310222736648235062495097951962638800384203417327|\
a9_8=-33475654618965607625490266678231230366345830527265525310030016230875755239420324728600957368877132012320553021/563607181486505082775419237614873016043169125130359330269370345097328575740250457475506596993780051575000000000|\
b8=63998419659074502960979467027044380533513499562179716145788153981258354882170183557294261050806789914954901252698438375102943237148428019253408419/14067383878224159980676935999500324267727723947516309536464061683130063437616452951326819824531063499969576146304645311203502485409773697867038000|\
b9=594562755257530592552224345703125000000/123802720910327682301431417435953442122031|\
b10=-1886691133979705639959153870454656/397947585885835383951563312487675|\
b11=-50061468875139778913910254637881/141186641936036184819986782313781|\
b_0=3635543992294021475202312550589/83920950031667636967840240807600|\
b_5=91472553308336221233020750122000000000/270389036674192633283981559366098375979|\
b_6=31245710879106859500854236329453125000000/125747467177230198813995913261775760852127|\
b_7=7580382785455138868782239796250000000000/33873008089978031564549954803189465564389|\
'
  for b in bList.split('|'):
    left, right = b.split('=') 
    print(GetFortranFloat(right))     







