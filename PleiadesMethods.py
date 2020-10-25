from odefortran import * 
from Common import *
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d
mp.dps = 40
 
FileName = sys.argv[0].replace('.py','')

# inital and final values 
y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
y0 = np.array([mpf(str(yEach)) for yEach in y0])
# t, y, outDict = solve_ivp_(0.,3.,y0,load=True,filename='Data/Reference.dat') 
# yfinal = y[-1]
yfinal = np.load('Data/PleiadesSolution_QP_rk1412Feagin_yfinal.npz',allow_pickle=True)['yfinal']
# mv Data/PleiadesSolution_yfinal.npz Data/PleiadesSolution_QP_rk1412Feagin_yfinal.npz

print('yfinal', type(yfinal), type(yfinal[0]))
# print(yfinal) 


fig, ax = plt.subplots() 
fig1, ax1 = plt.subplots() 
figfcn, axfcn = plt.subplots() 
figxy, axxy = plt.subplots(2)

tolArray = [1E-35, 1E-32, 1E-30, 1E-25, 1E-20, 1E-16, 1E-15, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]
tolArray = [1E-32, 1E-25, 1E-20, 1E-16, 1E-15, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]
tolArray = [1E-35, 1E-33, 1E-32, 1E-31, 1E-30, 1E-29, 1E-28, 1E-26, 1E-25, 1E-20, 1E-16, 1E-15, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]


methods = ['rk1412Feagin','rk108Feagin','rk1210Feagin','rk87EnrightVerner','rk65Dormand','rk1412Long','rk1211Peter','rk109Legendre','rk87Dormand']
methods = ['rk54Sharp','rk54Dormand','rk1412Feagin','rk108Feagin','rk1210Feagin','rk87EnrightVerner','rk65Dormand','rk1412Long','rk1211Peter','rk109Legendre','rk87Dormand']

methods = ['rk54Sharp','rk54Dormand','rk65Dormand','rk87Dormand','rk87EnrightVerner','rk108Feagin','rk109Legendre','rk1210Feagin','rk1211Peter','rk1412Long','rk1412Feagin']

# methods = ['rk1412Feagin','rk108Feagin','rk1210Feagin'] 

for i, method in enumerate(methods):
  scd = np.full(np.size(tolArray), -1E10) 
  mescd = np.full(np.size(tolArray), -1E10) 
  cpuTime = np.full(np.size(tolArray), -1E10) 
  fcn = np.full(np.size(tolArray), -1E10)  
  for j, atol in enumerate(tolArray):
    rtol = atol 
    filename = 'Data/%s_%s_rtol%s_atol%s.dat' % (FileName, method, PrintExponentialBetter(rtol), PrintExponentialBetter(atol)) 
    t, y, outDict = solve_ivp_(mpf('0.'),mpf('3.'),y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=mpf('1E-2'),max_step=mpf('1.'),min_step=mpf('1E-10'),Print6=True,load=True)
    if outDict['Success']:
      scd[j] = (np.max(abs(yfinal[:14]-y[-1][:14])/abs(yfinal[:14]))) 
      mescd[j] = np.max(abs(yfinal-y[-1])/(atol + rtol*abs(yfinal))) 
      cpuTime[j] = outDict['cpuTime'] 
      fcn[j] = outDict['Evaluated']
    #else:
    #  t, y, outDict = solve_ivp_(mpf('0.'),mpf('3.'),y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=mpf('1E-2'),max_step=mpf('1.'),min_step=mpf('1E-10'),Print6=True,load=True)
    print('rtol=%s, atol=%s, scd=%5.2f, mescd=%5.2f' % (PrintExponentialBetter(rtol), PrintExponentialBetter(atol), -np.log(scd[j])/np.log(10.), -np.log(mescd[j])/np.log(10.)))  
    if j == 0:
      #axxy[0].plot(y[:, IndexList.index('x1')], y[:, IndexList.index('y1')], color=color)
      #axxy[1].plot(y[:, IndexList.index('x3')], y[:, IndexList.index('y3')], color=color)
      #axxy[0].plot(t, y[:, IndexList.index('x1')], color=color)
      line, = axxy[0].plot(t, y[:, IndexList.index('x1')])
      color = line.get_color()
      axxy[1].plot(t, y[:, IndexList.index('y1')], color=color)
  marker = Marker13[np.mod(i,13)]
  ax.loglog(tolArray, scd, label=method, marker=marker, color=color) 
  ax1.semilogy(-np.log(scd)/np.log(10.), cpuTime, color=color, label=method, marker=marker) 
  # axfcn.loglog(fcn, scd, label=method, marker=marker, color=color) 
  log = lambda x: np.log(x)/np.log(10.)
  axfcn.plot(log(fcn),  log(scd), label=method, marker=marker, color=color)

ax.axis('equal')
addReference(ax, 1E-14, mode='x', ls='k:', alpha=0.7)
addReference(ax, 1E-14, mode='y', ls='k:', alpha=0.7)

ax.legend(loc='best',fontsize=10.)
ax.set_xlabel(r'atol')
ax.set_ylabel(r'$\left\Vert \delta y/y \right\Vert_{\infty}$')
fig.tight_layout()

fig.savefig('Plots/%s_error.png' % (FileName), dpi=300)
fig.clf()

ax1.legend(loc='best',fontsize=10.)
ax1.set_xlabel(r'scd')
ax1.set_ylabel(r'cpu (sec)')

fig1.tight_layout()
fig1.savefig('Plots/%s_scd_cpu.png' % (FileName), dpi=300)
fig1.clf()

for i in range(2):    
  axxy[i].legend(loc='best')
axxy[0].set_xlim(0., 3.)
axxy[1].set_xlim(0., 3.)
axxy[0].set_ylim(-1, 3)
axxy[1].set_ylim(-4, 4)

figxy.tight_layout()
figxy.savefig('Plots/%s_x1y1.png' % (FileName), dpi=300)
figxy.clf()


axfcn.set_xlim(3.03, 9.21) 
axfcn.set_ylim(-30.5,0.22) 
axfcn.xaxis.set_major_locator(MultipleLocator(1))
axfcn.xaxis.set_major_formatter(FormatStrFormatter('%d'))
axfcn.xaxis.set_minor_locator(AutoMinorLocator(n=5)) 
axfcn.yaxis.set_major_locator(MultipleLocator(5)) 
axfcn.yaxis.set_minor_locator(AutoMinorLocator(n=5)) 
axfcn.grid(True, which='major', color='k', alpha=0.4) 
axfcn.grid(True, which='minor', color='k', alpha=0.1)


axfcn.legend(loc='best',fontsize=10.)
axfcn.set_ylabel(r'$\log \left\Vert \delta y/y \right\Vert_{\infty}$')
axfcn.set_xlabel(r'$\log ( {\rm fcn})$')

figfcn.tight_layout()
figfcn.savefig('Plots/%s_scd_fcn.png' % (FileName), dpi=300)
figfcn.clf()

  
  


  
 

