import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
from scipy.special import sph_harm
import os

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

fig = plt.figure(frameon=False)
ii = 1
for l in range(5, 8):
    for m in range(2, 5):
        ax = fig.add_subplot(3,3,ii, projection='3d')
        ax.axis('off')
        s = sph_harm(m, l, theta, phi).real
        ax.plot_surface(x, y, z,  rstride=1, cstride=1, antialiased=False, 
            facecolors=[[cm.hot(b) for b in a] for a in colors.Normalize()(s)])
        ii = ii + 1

output_dir = 'plots/sph_harm'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'sph_harm.pdf'), 
      transparent=True, bbox_inches='tight')

