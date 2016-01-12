import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
from scipy.special import sph_harm
import os

#from matplotlib.colors import LightSource

output_dir = 'plots/sph_harm'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

cmap = cm.seismic

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

def plot_sph_harm(l, m, cross_out=False, output_dir=output_dir):
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.axis('off')
    
    s = sph_harm(m, l, theta, phi).real
    
    fc = [[cmap(b) for b in a] for a in colors.Normalize()(s)]
    #ls = LightSource(azdeg=315, altdeg=45)
    #fc = ls.shade(colors.Normalize()(s), cmap=cmap)
    
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, antialiased=False, 
        linewidth=0, facecolors=fc)
    
    fname = str(l)+'_'+str(m)
    if cross_out:
        ax.text(x=0, y=0, z=0, s='x', horizontalalignment='center',
            verticalalignment='center', fontsize=240)
        fname+='x'
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, fname+'.png'), 
        transparent=True, bbox_inches='tight', dpi=400)
    plt.close()

fig = plt.figure(frameon=False)
ii = 1
for l in range(0, 4):
    for m in range(-l, l+1):
        plot_sph_harm(l, m)
        print(l, m)
        if (m != 0):
            plot_sph_harm(l, m, cross_out=True)
            print(l, m, 'x')

