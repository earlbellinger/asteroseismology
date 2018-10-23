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

#cmap = cm.seismic
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "black", "white"])

r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

def plot_sph_harm(l, m, cross_out=False, output_dir=output_dir):
    plt.rcParams['figure.figsize'] = 4.8, 4.8
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.axis('off')
    
    s = sph_harm(m, l, theta, phi).real
    
    fc = [[cmap(b) for b in a] for a in colors.Normalize()(s)]
    #ls = LightSource(azdeg=315, altdeg=45)
    #fc = ls.shade(colors.Normalize()(s), cmap=cmap)
    
    #ax.plot_surface(x, y, z,  rstride=1, cstride=1, antialiased=True, 
    #    linewidth=0, facecolors=fc)
    #ax.contour(x, y, z)
    ax.plot_surface(x, y, z, rstride=1, cstride=1, antialiased=True, 
        linewidth=1, facecolors=fc)
    
    fname = str(l)+'_'+str(m)
    #if cross_out:
    #    ax.text(x=0, y=0, z=0, s='x', horizontalalignment='center',
    #        verticalalignment='center', fontsize=240)
    #    fname+='x'
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, fname+'.pdf'), 
        transparent=False, bbox_inches='tight')#, dpi=400)
    plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.close()

#fig = plt.figure(frameon=False)
ii = 1
m = 0
for l in [0, 1, 2, 3, 20, 25, 75]:
    #for m in range(-l, l+1):
    plot_sph_harm(l, m)
    #print(l, m)
        #if (m != 0):
        #    plot_sph_harm(l, m, cross_out=True)
        #    print(l, m, 'x')



output_dir = 'plots/sph_harm'
theta = np.linspace(0, np.pi, 20) # polar angle
phi = np.linspace(0, 2*np.pi, 20) # azimuth angle
phi, theta = np.meshgrid(phi, theta)
for degree in range(5):
    s = sph_harm(0, degree, theta, phi).real
    #Ymn = legendre(degree, np.cos(theta(:,1)))
    #Ymn = Ymn(order+1,:)'
    #yy = Ymn
    
    #for kk = 2: size(theta,1)
    #    yy = [yy Ymn];
    #end
    
    #yy = yy.*cos(order*phi);
    
    #order2 = max(max(abs(yy)));
    #rho = radius + amplitude*yy/order2;
    
    r = s*np.sin(theta)
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    z = s*np.cos(theta)
    
    plt.rcParams['figure.figsize'] = 4.8, 4.8
    fig = plt.figure(frameon=False)
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.axis('off')
    
    cmap = #matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "black", "white"])
    fc = [[cmap(b) for b in a] for a in colors.Normalize()(s)]
    ax.plot_surface(x, y, z, rstride=1, cstride=1, antialiased=True, 
        linewidth=1, facecolors=fc)
    
    fname = str(degree)
    
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, fname+'.pdf'), 
        transparent=False, bbox_inches='tight')#, dpi=400)
    plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.close()










import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm 
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm

theta_1d = np.linspace(0,   np.pi,  91) # 2 GRAD Schritte
phi_1d   = np.linspace(0, 2*np.pi, 181) # 2 GRAD Schritte
theta_2d, phi_2d = np.meshgrid(theta_1d, phi_1d)
xyz_2d = np.array([np.sin(theta_2d) * np.sin(phi_2d),
                  np.sin(theta_2d) * np.cos(phi_2d),
                  np.cos(theta_2d)]) 

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#F97100", "white", "#0571b0"])
colormap = cm.ScalarMappable( cmap=cmap )#plt.get_cmap("seismic"))
colormap.set_clim(-.45, .45)
limit = .5

plt.rcParams['figure.figsize'] = 4.8, 4.8
fig = plt.figure(frameon=False)
for l in range(0,4):
    print("Y_%i_%i" % (l,0)) 
    ax = fig.add_subplot(1, 4, l+1, projection='3d')
    ax.axis('off')
    #plt.figure(frameon=False)
    #ax = plt.gca(projection = "3d")
    #plt.title("$Y^{%i}_{%i}$" % (m,l))
    Y_lm = sph_harm(0,l, phi_2d, theta_2d)
    r = np.abs(Y_lm.real)*xyz_2d
    ax.plot_surface(r[0], r[1], r[2], 
                    facecolors=colormap.to_rgba(Y_lm.real), 
                    rstride=2, cstride=2)
    ax.set_xlim(-limit,limit)
    ax.set_ylim(-limit,limit)
    ax.set_zlim(-limit,limit)
    ax.set_aspect("equal")
    #ax.set_axis_off()

fig.tight_layout()
fig.savefig(os.path.join('plots', 'sph.png'), 
    transparent=False, bbox_inches='tight', dpi=400)
plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
plt.axis([-1.1, 1.1, -1.1, 1.1])
plt.close()








import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm 
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm


theta_1d = np.linspace(0,   np.pi,  91) # 2 GRAD Schritte
phi_1d   = np.linspace(0, 2*np.pi, 181) # 2 GRAD Schritte
theta_2d, phi_2d = np.meshgrid(theta_1d, phi_1d)
xyz_2d = np.array([np.sin(theta_2d) * np.sin(phi_2d),
                   np.sin(theta_2d) * np.cos(phi_2d),
                   np.cos(theta_2d)]) 

cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", 
    ["#F97100", "white", "#0571b0"])
colormap = cm.ScalarMappable( cmap=cmap )#plt.get_cmap("seismic"))
colormap.set_clim(-.45, .45)
limit = 1#2#.5

scales = np.concatenate((np.linspace(0, 1, 10), 
    np.linspace(1, -1, 20),
    np.linspace(-1, 0, 10)))
for i in range(len(scales)):
    print(i)
    scale = scales[i]
    plt.rcParams['figure.figsize'] = 16, 9 #4.8, 4.8
    fig = plt.figure(frameon=False)
    for l in range(0,4):
        print("Y_%i_%i" % (l,0)) 
        ax = fig.add_subplot(1, 4, l+1, projection='3d')
        ax.axis('off')
        #plt.figure(frameon=False)
        #ax = plt.gca(projection = "3d")
        #plt.title("$Y^{%i}_{%i}$" % (m,l))
        Y_lm = sph_harm(0, l, phi_2d, theta_2d)
        #rho = 3+scales[i]
        #xyz_2d = np.array([np.sin(theta_2d) * np.sin(phi_2d),
        #                   np.sin(theta_2d) * np.cos(phi_2d),
        #                   np.cos(theta_2d)]) 
        rho = 0.3 + scales[i] * np.abs(Y_lm.real) #/ np.max(np.abs(Y_lm.real))
        r = rho*xyz_2d
        ax.plot_surface(r[0], r[1], r[2], 
                        facecolors=colormap.to_rgba(Y_lm.real),# * scales[i]), 
                        rstride=2, cstride=2)
        ax.set_xlim(-limit,limit)
        ax.set_ylim(-limit,limit)
        ax.set_zlim(-limit,limit)
        ax.set_aspect("equal")
        #ax.set_axis_off()
    
    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'anim_sph_py', 'sph'+str(i)+'.png'), 
        transparent=False, bbox_inches='tight', dpi=400)
    plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.close()







phi = np.linspace(0, np.pi, 10)#0)
theta = np.linspace(0, 2*np.pi, 10)#0)
phi, theta = np.meshgrid(phi, theta)

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

m, l = 0, 2

# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
Yml = sph_harm(m, l, theta, phi).real
fmax, fmin = Yml.max(), Yml.min()
fcolors = (Yml - fmin)/(fmax - fmin)

# Set the aspect ratio to 1 so our sphere looks spherical
for i in range(len(scales)):
    fig = plt.figure(figsize=plt.figaspect(1.))
    ax = fig.add_subplot(111, projection='3d')
    
    rho = 7 + scales[i]*30*Yml #/order2
    
    r = rho * np.sin(theta);
    x = r * np.cos(theta) #np.sin(phi) * np.cos(theta)
    y = r * np.sin(theta) #np.sin(phi) * np.sin(theta)
    z = rho * np.cos(phi)
    
    ax.plot_surface(x, y, z,  rstride=1, cstride=1, 
        facecolors=cm.seismic(fcolors))# * scales[i]))
    # Turn off the axis planes
    ax.set_axis_off()
    fig.tight_layout()
    
    fig.savefig(os.path.join('plots', 'anim_sph_py', 'sph'+str(i)+'.png'), 
        transparent=False, bbox_inches='tight', dpi=50)#400)
    plt.subplots_adjust(top=1,bottom=0,left=0,right=1)
    plt.axis([-1.1, 1.1, -1.1, 1.1])
    plt.close()


