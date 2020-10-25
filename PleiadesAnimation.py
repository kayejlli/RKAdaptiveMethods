"""
Matplotlib Animation Example

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from scipy.interpolate import InterpolatedUnivariateSpline as newinterp1d
FileName = sys.argv[0].replace('.py','')
Color7 = ['b', 'r', 'g', 'm', 'orange', 'blueviolet','salmon']

IndexList = []
for Name in ['x', 'y', 'dx', 'dy']:
  for i in range(7):
    IndexList.append('%s%d' % (Name, i+1))


# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
ax = fig.add_subplot(111, aspect='equal', autoscale_on=False,
       xlim=(-5.4, 5.4), ylim=(-5.4, 5.4))
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('equal')
ax.text(-5.1, 4.9, r'${\rm The \, Pleiades}$')  
ms = 2
for i in range(1,8):
  globals()['particles%d'%i], = ax.plot([], [], color=Color7[i-1], marker='o', ms=i, markevery=-1) 
timeLabel = ax.text(-5.1, -5.1, r'$t=%5.2f \, {\rm s}$' % (0.))
# particles, = ax.plot([], [], 'bo', ms=ms)

bounds = [-5.24, 5.24, -5.24, 5.24] 
# rect is the box edge
rect = plt.Rectangle(bounds[::2],
       bounds[1] - bounds[0],
       bounds[3] - bounds[2],
       ec='none', lw=1.5, fc='none')
ax.add_patch(rect)

def init():
  """initialize animation"""
  global box, rect
  for i in range(1,8):
    globals()['particles%d'%i].set_data([], []) 
  rect.set_edgecolor('none')
  timeLabel.set_text('') 
  # return particles, rect
  return particles1, particles2, particles3, particles4, particles5, particles6, particles7, rect, timeLabel

# load data 
data = np.load('Data/test.npz',allow_pickle=True)
t = data['t']
y = data['y']
for i, Name in enumerate(IndexList):
  if 'd' in Name:
    break
  globals()['%s_fit' % (Name)] = newinterp1d(t, y[:, i], k=5)
  print('min = %5.2f, max = %5.2f' % (min(y[:, i]), max(y[:, i]))) 
 
  
# animation function.  This is called sequentially
def animate(i):
  time = np.linspace(0., 3., 1000)[i]
  timeLabel.set_text(r'$t=%5.2f \, {\rm s}$' % (time)) 
  if i == 0:
    xArray, yArray = [], []
    for i in range(1,8): 
      color = Color7[i-1] 
      x_pos, y_pos = globals()['x%d_fit' % i](time), globals()['y%d_fit' % i](time)
      xArray.append(x_pos) 
      yArray.append(y_pos) 
      globals()['particles%d'%i].set_data(x_pos, y_pos) 
      globals()['particles%d'%i].set_markersize(i) 
    # update pieces of the animation
    rect.set_edgecolor('k')
    # particles.set_data(xArray, yArray)
    # particles.set_markersize(ms)
    return particles1, particles2, particles3, particles4, particles5, particles6, particles7, rect, timeLabel
    #return particles, rect
  else:
    tArray = np.linspace(0., time, 5000) 
    for i in range(1,8): 
      color = Color7[i-1] 
      x_pos, y_pos = globals()['x%d_fit' % i](tArray), globals()['y%d_fit' % i](tArray)
      globals()['particles%d'%i].set_data(x_pos, y_pos) 
      globals()['particles%d'%i].set_markersize(i) 
      globals()['particles%d'%i].set_markevery([4999]) 
    # update pieces of the animation
    rect.set_edgecolor('k')
    # particles.set_data(xArray, yArray)
    # particles.set_markersize(ms)
    return particles1, particles2, particles3, particles4, particles5, particles6, particles7, rect, timeLabel
       


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
       frames=1000, interval=3E3/1000., blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('%s.mp4' % (FileName), fps=1000, extra_args=['-vcodec', 'libx264'],dpi=300)

# plt.show()
