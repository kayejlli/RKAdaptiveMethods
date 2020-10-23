from odefortran import * 
from Common import *
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d

# solve_ivp_(0., 3., y0, filename='Data/test.dat', atol=1E-6, rtol=1E-7, first_step=1E-2, max_step=1., min_step=1E-4, Print6=True) 
 
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

fig, ax = plt.subplots() 
fig1, ax1 = plt.subplots() 

tolArray = [1E-16, 1E-15, 1E-14, 1E-13, 1E-12, 1E-11, 1E-10, 1E-9, 1E-8, 1E-6, 1E-5]
# tolArray = tolArray[::3]
scd = np.zeros((np.size(tolArray),np.size(tolArray))) 
mescd = np.zeros((np.size(tolArray),np.size(tolArray))) 
cpuTime = np.zeros((np.size(tolArray),np.size(tolArray))) 
for i, atol in enumerate(tolArray):
  for j, rtol in enumerate(tolArray):
    filename = 'Data/%s_rtol%s_atol%s.dat' % (FileName, PrintExponentialBetter(rtol), PrintExponentialBetter(atol)) 
    t, y, outDict = solve_ivp_(0.,3.,y0,method='rk108Feagin',filename=filename,atol=atol,rtol=rtol,first_step=1E-2,max_step=1.,min_step=1E-5,Print6=False) 
    #print('tfinal=%e, err=%e' % (t[-1], abs(t[-1]-3.))) 
    #x1Inter = newinterp1d(t, y[:, 0], k=5)
    #print(y[0])
    #print('err=%e' % (np.max(abs(yfinal-y[-1])/abs(yfinal)))) 
    scd[i, j] = (np.max(abs(yfinal[:14]-y[-1][:14])/abs(yfinal[:14]))) 
    mescd[i, j] = np.max(abs(yfinal-y[-1])/(atol + rtol*abs(yfinal))) 
    cpuTime[i, j] = outDict['cpuTime'] 
    # print('mescd[i, j]=%e' % (mescd[i, j])) 
    print('rtol=%s, atol=%s, scd=%5.2f, mescd=%5.2f' % (PrintExponentialBetter(rtol), PrintExponentialBetter(atol), -np.log(scd[i, j])/np.log(10.), -np.log(mescd[i, j])/np.log(10.)))  

for i in range(np.size(tolArray)):
  color = NameColorList[i].replace(' ','')
  ax.loglog(tolArray, scd[i, :], color=color, label=r'$atol=%s$' % (LabelExponentialBetter(tolArray[i]))) 
  ax.loglog(tolArray, scd[:, i], color=color, ls='--')
  ax1.semilogy(-np.log(scd[i,:])/np.log(10.), cpuTime[i,:], color=color)  
  ax1.semilogy(-np.log(scd[:,i])/np.log(10.), cpuTime[:,i], color=color, ls='--')  

ax.loglog([], [], 'k--', label=r'$rtol$') 
ax.axis('equal')
# ax.loglog(tolArray, tolArray, 'k--', alpha=0.7) 
addReference(ax, 1E-14, mode='x', ls='k:', alpha=0.7)
addReference(ax, 1E-14, mode='y', ls='k:', alpha=0.7)

ax.legend(loc='best',fontsize=10.)
ax.set_xlabel(r'tol')
ax.set_ylabel(r'tol')
fig.tight_layout()

fig.savefig('Plots/%s_error.png' % (FileName), dpi=300)
fig.clf()

ax1.set_xlabel(r'scd')
ax1.set_ylabel(r'cpu (sec)')

fig1.tight_layout()
fig1.savefig('Plots/%s_scd_cpu.png' % (FileName), dpi=300)
fig1.clf()

#    ls = ['-', '--', '-.'][j]
    

  
  


  
 

