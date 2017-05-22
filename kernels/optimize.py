import numpy as np
import pandas as pd
from scipy import integrate, optimize
from scipy.interpolate import interp1d
import os

observed_freqs = pd.DataFrame(np.loadtxt('Sun-freqs.dat', skiprows=1), 
                     columns=['n', 'l', 'nu_s', 'sigma'])
observed_freqs.head()

model_freqs = pd.DataFrame(np.loadtxt(os.path.join('modelS', 'modelS.freq')),
                   columns=['l', 'n', 'nu_m', 'Q'])
model_freqs.head()

ell_0 = model_freqs[model_freqs['l']==0]
Q_0 = interp1d(ell_0['nu_m'], ell_0['Q'], 
               fill_value='extrapolate')
Q_norm = model_freqs['Q'] / Q_0(model_freqs['nu_m'])
model_freqs['Q_norm'] = Q_norm
model_freqs.head()

nu = pd.merge(observed_freqs, model_freqs)
nu.head()

dnu = (nu['nu_m'] - nu['nu_s']) / nu['nu_m']
nu['dnu'] = dnu
nu.head()

f1 = r'c^2'
f2 = r'\rho'
#f2 = r'%5Crho'

nls = [(int(n), int(l)) 
       for (n,l) in np.array(nu[['n', 'l']])]
       
K = {(int(n),int(l)) :
     pd.DataFrame(np.loadtxt(os.path.join('modelS_ker',
         f1+'-'+f2+'_l='+str(int(l))+'_n='+str(int(n))+'.dat')),
                  columns=['x', f1, f2])
     for (n,l) in nls}

n,l = nls[10]
K[(n,l)].head()

model = pd.DataFrame(np.loadtxt(f1+'_'+f2+'.dat'), 
                      columns=['r', f1, f2])
x = model['r']
f1_m = model[f1]
f2_m = model[f2]
model.head()

def A(n, l, df1, df2):
    kernel = K[n,l]
    fx = kernel[f1] * df1 / f1_m + \
         kernel[f2] * df2 / f2_m
    return integrate.simps(fx, x=x)

def F_surf(n, l, a, nu_ac=5000):
    nu_per_ac = nu['nu_m'][nls.index((n,l))] / nu_ac
    Q_norm = nu['Q_norm'][nls.index((n,l))]
    return ( a[0] * nu_per_ac**-1 + a[1] * nu_per_ac**3 ) / Q_norm

def chi_sq(df1, df2, a):
    return sum(
        ( nu['dnu'][nls.index((n,l))] - A(n,l,df1,df2) - F_surf(n,l,a) )**2 \
        / nu['sigma'][nls.index((n,l))]
               for (n,l) in nls)

def L_2(df1, df2):
    dx = np.diff(x, 2)
    d2f1_dr2 = np.diff(df1, n=2) / dx
    d2f2_dr2 = np.diff(df2, n=2) / dx
    return integrate.simps(d2f1_dr2**2 + d2f2_dr2**2, x[2:])

def objective(xs):
    df1 = xs[:len(f1_m)]
    df2 = xs[len(f1_m):len(f1_m)+len(f2_m)]
    a = xs[-3:-1]
    alpha = xs[-1]
    return chi_sq(df1, df2, a) + alpha * L_2(df1, df2)

df1 = np.zeros(len(x))
df2 = np.zeros(len(x))
a = [0,0]
alpha = 0 

A(nls[0][0], nls[0][1], df1, df2) # should be 0 
chi_sq(f1_m, f2_m, a) # should be close to 0 
L_2(df1, df2) # should be 0 

x0 = list(df1) + list(df2) + list(a) + [alpha]
print(objective(x0)) # should be small

bounds = [(None, None) for i in range(len(x0)-1)] + [(0, None)]
res = optimize.minimize(objective, x0, bounds=bounds, 
    options={'disp':True})
print(res)
