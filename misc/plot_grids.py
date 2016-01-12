import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import numpy as np
import os

from sys import path
path.append('../scripts')
from sobol_lib import i4_sobol

np.random.seed(12)

output_dir = 'plots/grids'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

points = [100, 400, 2000]

### Random grid
fig = plt.figure(figsize=(2.35037, 4.17309), dpi=400)
xs = np.random.random(points[-1])
ys = np.random.random(points[-1])
zs = np.random.random(points[-1])
for ii, n_points in enumerate(points[::-1]):
    if ii>0:
        xs = np.random.choice(xs, size=n_points, replace=False)
        ys = np.random.choice(ys, size=n_points, replace=False)
        zs = np.random.choice(zs, size=n_points, replace=False)
    ax = fig.add_subplot(3,1,3-ii)
    plt.scatter(xs, ys, c=zs, 
                cmap=cm.gist_heat,
                linewidth=0.25,
                s=10/(3-ii))
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.xlim([0, 1])
    plt.ylim([0, 1])

fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'grid-random.png'), 
      bbox_inches='tight', dpi=400)

### Quasi-random grid
fig = plt.figure(figsize=(2.35037, 4.17309), dpi=400)
for ii, n_points in enumerate(points):
    ax = fig.add_subplot(3,1,ii+1)
    grid = np.transpose([np.array(i4_sobol(3, i)[0]) 
                         for i in range(20000, 20000+n_points)])
    plt.scatter(grid[0], grid[1], c=grid[2],
                cmap=cm.gist_heat,
                #linewidth=0.5/(1+ii),
                linewidth=0.25,
                s=10/(1+ii))
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.xlim([0, 1])
    plt.ylim([0, 1])

fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'grid-quasirandom.png'), 
      bbox_inches='tight', dpi=400)

### Linear
fig = plt.figure(figsize=(2.35037, 4.17309), dpi=400)
for ii, n_points in enumerate(points):
    ax = fig.add_subplot(3,1,ii+1)
    nmax = round(n_points**(1/3))
    grid = np.linspace(0.01, 0.99, nmax)
    x, y, z = np.meshgrid(*[grid]*3)
    #plt.scatter(np.repeat([grid], nmax, axis=0).flatten(), 
    #            np.repeat([grid], nmax, axis=1).flatten(), 
    #            color="black", marker=".", s=0.2)
    plt.scatter(x.flatten(), 
                y.flatten(), 
                c='black', 
                #cmap=cm.gist_heat, 
                #linewidth=0.5/(1+ii), 
                linewidth=0.25,
                s=10/(1+ii))
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.xlim([0, 1])
    plt.ylim([0, 1])

fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'grid-linear.png'), 
      bbox_inches='tight', dpi=400)

