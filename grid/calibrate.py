#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

import numpy as np
import pandas as pd
from scipy import optimize
from os import path
from subprocess import Popen
from math import log10

Z_div_X_solar = 0.02293 # GS98 # 0.0245 # GN93 #
log10_Z_div_X_solar = np.log10(Z_div_X_solar) 
constraint_names = ("log L", "log R", "Fe/H") 
param_names = ("Y", "alpha", "Z") 
param_init = [0.273449170177157, 1.83413390909832, 0.0197444964340224] 
directory = 'calibrate_py'
print(directory)

def objective():
    ## minimize sum(log(model values / solar values)**2) 
    # searches in LOGS_MS subdirectory of the global 'directory' variable 
    
    hstry_file = path.join(directory, 'LOGS_MS', 'history.data')
    if (not path.exists(hstry_file)): 
        return np.inf 
    hstry = pd.read_table(hstry_file, header=0, skiprows=5, delimiter='\s+') #header=1, 
    mdl = hstry.loc[hstry.shape[0]-1] #hstry[nrow(hstry),]
    
    # [Fe/H] = log10 ( Z / X / (Z/X)_Sun )
    mdl_Fe_H = mdl['log_surf_cell_z']-np.log10(mdl['surface_h1'])-log10_Z_div_X_solar 
    mdl_vals = [mdl['log_L'], mdl['log_R'], mdl_Fe_H]
    
    print("*** Model values") 
    print(constraint_names, mdl_vals)
    print('L', 10**mdl['log_L'], 'R', 10**mdl['log_R']) 
    
    result = sum([ii**2 for ii in mdl_vals]) 
    if np.isfinite(result): 
        return log10(result) 
    return 10**10 

### SEARCH
iteration = 0
best_val = np.inf
best_param = param_init

#run = function(params) {
def run(params): 
    global iteration
    global best_val
    global best_param
    
    iteration = iteration + 1
    print("**** iter:", iteration)
    
    Y, alpha, Z = params 
    print(param_names, (Y, alpha, Z))
    
    if (Y < 0.2 or Y > 0.4 or Z < 0 or Z > 0.04 or alpha < 1 or alpha > 3):
        return 10**10
    
    #if (Y < 0.23):
    #    Y = 0.23
    #if (Y > 0.33):
    #    Y = 0.33
    #if (Z < 0.01):
    #    Z = 0.01
    #if (Z > 0.04):
    #    Z = 0.04
    #if (alpha < 1):
    #    alpha = 1
    #if (alpha > 3):
    #    alpha = 3
    
    command = "./dispatch.sh" + \
        ' -Y ' + str(Y) + \
        ' -a ' + str(alpha) + \
        ' -o ' + '0' + \
        ' -Z ' + str(Z) + \
        ' -D ' + '1' + \
        ' -g ' + '1' + \
        ' -e ' + '0' + \
        ' -c ' + "4572000000" + \
        ' -d ' + directory
    print(command)
    #system(command)
    process = Popen(command.split(), shell=False)
    process.wait()
    
    obj_val = objective() 
    print("**** objective value =", obj_val) 
    
    if (obj_val < best_val):
        best_val = obj_val
        print("*****", param_names, params)
        best_param = params
        print("***** New record!")
    
    return obj_val

result = optimize.minimize(fun=run, x0=param_init, method='Nelder-Mead', 
    options={'disp': True,
             'maxiter': 10000}) #,
    #bounds=((0.25, 0.32), (1, 3), (0.012, 0.03))) 

print("Optimization terminated. Saving best result")
Y, alpha, Z = result.x
command = "./dispatch.sh" + \
    ' -Y ' + str(Y) + \
    ' -a ' + str(alpha) + \
    ' -o ' + '0' + \
    ' -Z ' + str(Z) + \
    ' -D ' + '1' + \
    ' -g ' + '1' + \
    ' -e ' + '0' + \
    ' -c ' + "4572000000" + \
    ' -d ' + directory
print(command)
process = Popen(command.split(), shell=False)
process.wait()
print(result)

