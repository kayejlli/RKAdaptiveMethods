# ref: https://archimede.dm.uniba.it/~testset/report/plei.pdf
from odef2py import *
from Python.Common import *
import time 
 
start_time = time.time()
FileName = sys.argv[0].replace('.py','')

# set inital and final values 
y0 = np.array([3,3,-1,-3,2,-2,2,3,-3,2,0,0,-4,4,0,0,0,0,0,1.75,-1.5,0,0,0,-1.25,1,0,0])
yfinal = np.load('DataSaved/PleiadesSolution_QP_rk1412Feagin_yfinal.npz',allow_pickle=True)['yfinal']

fig, ax = plt.subplots() 
fig1, ax1 = plt.subplots() 
figfcn, axfcn = plt.subplots() 
fign, axn = plt.subplots() 

tolArray = [1E-17, 1E-16, 3E-16, 1E-15, 3E-15, 1E-14, 3E-14, 1E-13, 3E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]

methods = ['rk54Sharp','rk54Dormand','rk65Dormand','rk87Dormand','rk87EnrightVerner','rk108Feagin','rk109Legendre','rk1210Feagin','rk1211Peter','rk1412Feagin']

for i, method in enumerate(methods):
  scd = np.full(np.size(tolArray), -1E10) 
  mescd = np.full(np.size(tolArray), -1E10) 
  cpuTime = np.full(np.size(tolArray), -1E10) 
  fcn = np.full(np.size(tolArray), -1E10) 
  n = np.full(np.size(tolArray), -1E10)
  for j, atol in enumerate(tolArray):
    rtol = atol 
    filename = 'Data/%s_%s_rtol%s_atol%s.dat' % (FileName, method, PrintExponentialBetter(rtol), PrintExponentialBetter(atol)) 
    t, y, outDict, filename = solve_ivp_(0.,3.,y0,method=method,filename=filename,atol=atol,rtol=rtol,first_step=1E-3,max_step=1.,min_step=1E-10,Print6=False) 
    if outDict['Success']:
      scd[j] = (np.max(abs(yfinal[:14]-y[-1][:14])/abs(yfinal[:14]))) 
      mescd[j] = np.max(abs(yfinal-y[-1])/(atol + rtol*abs(yfinal))) 
      cpuTime[j] = outDict['cpuTime'] 
      fcn[j] = outDict['Evaluated'] 
      n[j] = outDict['TotalSteps'] 
    else:
      for name in outDict:
        print(name, outDict[name])  
      raise ValueError('The program did not end properly') 
    print('rtol=%s, atol=%s, scd=%5.2f, mescd=%5.2f' % (PrintExponentialBetter(rtol), PrintExponentialBetter(atol), -np.log(scd[j])/np.log(10.), -np.log(mescd[j])/np.log(10.)))  
  color = 'C%d' % (np.mod(i, 10)) 
  marker = Marker13[np.mod(i, 13)]
  ax.loglog(tolArray, scd, label=method, marker=marker, color=color) 
  ax1.semilogy(-np.log(scd)/np.log(10.), cpuTime, color=color, label=method, marker=marker) 
  log = lambda x: np.log(x)/np.log(10.) 
  # axfcn.loglog(fcn,  scd, label=method, marker=marker, color=color) 
  axfcn.plot(log(fcn),  log(scd), label=method, marker=marker, color=color) 
  axn.plot(log(n),  log(scd), label=method, marker=marker, color=color) 

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



axfcn.set_xlim(3.03, 9.21) 
axfcn.set_ylim(-30.5,0.22) 
axfcn.xaxis.set_major_locator(MultipleLocator(1))
axfcn.xaxis.set_major_formatter(FormatStrFormatter('%d'))
axfcn.xaxis.set_minor_locator(AutoMinorLocator(n=5)) 
axfcn.yaxis.set_major_locator(MultipleLocator(5)) 
axfcn.yaxis.set_minor_locator(AutoMinorLocator(n=5)) 

axfcn.set_xlim(3.03, 6.2)
axfcn.set_xlim(3.4, 5.5)
axfcn.set_ylim(-12.5, -1) 
axfcn.xaxis.set_major_locator(MultipleLocator(0.5))
axfcn.xaxis.set_minor_locator(AutoMinorLocator(n=5)) 
axfcn.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
axfcn.set_yticks([-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1], minor=False) 
#axfcn.set_yticks([-12,-11,-10,-9,-8,-7,-6], minor=False)
axfcn.yaxis.set_minor_locator(AutoMinorLocator(n=2)) 
axfcn.grid(True, which='major', color='k', alpha=0.4) 
axfcn.grid(True, which='minor', color='k', alpha=0.1)

axfcn.legend(loc=1,fontsize=10.)
axfcn.set_xlabel(r'$\log({\rm fcn})$')
axfcn.set_ylabel(r'$\log \left\Vert \delta y/y \right\Vert_{\infty}$')

figfcn.tight_layout()
figfcn.savefig('Plots/%s_scd_fcn.png' % (FileName), dpi=300)
print('Saving to [Plots/%s_scd_fcn.png]' % (FileName)) 
figfcn.clf()
 
axn.legend(loc=1,fontsize=10.)
axn.set_xlabel(r'$\log({\rm N})$')
axn.set_ylabel(r'$\log \left\Vert \delta y/y \right\Vert_{\infty}$')

fign.tight_layout()
fign.savefig('Plots/%s_scd_n.png' % (FileName), dpi=300)
print('Saving to [Plots/%s_scd_n.png]' % (FileName)) 
fign.clf()


end_time = time.time()
print('The time it takes is %.2f sec' % (end_time-start_time)) 
