from odefortran import * 
from Common import *
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d
 
FileName = sys.argv[0].replace('.py','')

# inital and final values 
y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
yfinal = [0.3706139143970502,0.3237284092057233E1,-0.3222559032418324E1,0.6597091455775310,\
0.3425581707156584,0.1562172101400631E1,-0.7003092922212495,-0.3943437585517392E1,\
-0.3271380973972550E1,0.5225081843456543E1,-0.2590612434977470E1,0.1198213693392275E1,\
-0.2429682344935824,0.1091449240428980E1,0.3417003806314313E1,0.1354584501625501E1,\
-0.2590065597810775E1,0.2025053734714242E1,-0.1155815100160448E1,-0.8072988170223021,\
0.5952396354208710,-0.3741244961234010E1,0.3773459685750630,0.9386858869551073,\
0.3667922227200571,-0.3474046353808490,0.2344915448180937E1,-0.1947020434263292E1]
yfinal = np.array(yfinal)

yfinal = np.load('Data/PleiadesSolution_QP_rk1412Feagin_yfinal.npz',allow_pickle=True)['yfinal']

fig, ax = plt.subplots() 
fig1, ax1 = plt.subplots() 
figfcn, axfcn = plt.subplots() 
figxy, axxy = plt.subplots(2)

tolArray = [1E-16, 3E-16, 1E-15, 3E-15, 1E-14, 3E-14, 1E-13, 3E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]

methods = ['rk54Sharp','rk54Dormand','rk65Dormand','rk87Dormand','rk87EnrightVerner','rk108Feagin','rk109Legendre','rk1210Feagin','rk1211Peter','rk1412Long','rk1412Feagin']

for i, method in enumerate(methods):
  scd = np.full(np.size(tolArray), -1E10) 
  mescd = np.full(np.size(tolArray), -1E10) 
  cpuTime = np.full(np.size(tolArray), -1E10) 
  fcn = np.full(np.size(tolArray), -1E10) 
  for j, atol in enumerate(tolArray):
    rtol = 1E-20
    filename = 'Data/%s_%s_rtol%s_atol%s.dat' % (FileName, method, PrintExponentialBetter(rtol), PrintExponentialBetter(atol)) 
    t, y, outDict = solve_ivp_(0.,3.,y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=1E-2,max_step=1.,min_step=1E-10,Print6=False) 
    if outDict['Success']:
      scd[j] = (np.max(abs(yfinal[:14]-y[-1][:14])/abs(yfinal[:14]))) 
      mescd[j] = np.max(abs(yfinal-y[-1])/(atol + rtol*abs(yfinal))) 
      cpuTime[j] = outDict['cpuTime'] 
      fcn[j] = outDict['Evaluated'] 
    else:
      t, y, outDict = solve_ivp_(0.,3.,y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=1E-2,max_step=1.,min_step=1E-5,Print6=True) 
    print('rtol=%s, atol=%s, scd=%5.2f, mescd=%5.2f' % (PrintExponentialBetter(rtol), PrintExponentialBetter(atol), -np.log(scd[j])/np.log(10.), -np.log(mescd[j])/np.log(10.)))  
    if j == 0:
      #axxy[0].plot(y[:, IndexList.index('x1')], y[:, IndexList.index('y1')], color=color)
      #axxy[1].plot(y[:, IndexList.index('x3')], y[:, IndexList.index('y3')], color=color)
      #axxy[0].plot(t, y[:, IndexList.index('x1')], color=color)
      line, = axxy[0].plot(t, y[:, IndexList.index('x1')])
      color = line.get_color()
      axxy[1].plot(t, y[:, IndexList.index('y1')], color=color)
  marker = Marker13[np.mod(13, i)]
  ax.loglog(tolArray, scd, label=method, marker=marker, color=color) 
  ax1.semilogy(-np.log(scd)/np.log(10.), cpuTime, color=color, label=method, marker=marker) 
  log = lambda x: np.log(x)/np.log(10.) 
  # axfcn.loglog(fcn,  scd, label=method, marker=marker, color=color) 
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

axfcn.set_xlim(3.03, 6.2)
axfcn.set_ylim(-12.5, -1) 
axfcn.xaxis.set_major_locator(MultipleLocator(0.5))
axfcn.xaxis.set_minor_locator(AutoMinorLocator(n=5)) 
axfcn.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axfcn.set_yticks([-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1], minor=False) 
axfcn.yaxis.set_minor_locator(AutoMinorLocator(n=2)) 
axfcn.grid(True, which='major', color='k', alpha=0.4) 
axfcn.grid(True, which='minor', color='k', alpha=0.1)

axfcn.legend(loc='best',fontsize=10.)
axfcn.set_xlabel(r'fcn')
axfcn.set_ylabel(r'$\log \left\Vert \delta y/y \right\Vert_{\infty}$')

figfcn.tight_layout()
figfcn.savefig('Plots/%s_scd_fcn.png' % (FileName), dpi=300)
figfcn.clf()
 
  


  
 

