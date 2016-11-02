import numpy as np
import os

kernel_info = np.loadtxt('kernels', skiprows=1)
ell = kernel_info[:,0]
nnn = kernel_info[:,1]

fname = 'rho-Gamma1' # 'c2-rho' # 'rho-c2' # 

os.makedirs(fname)

mode = 0
kernel = []
with open(fname+'.dat', 'r') as f:
    while True:
        line = f.readline()
        if line == '':    # done
            break
        elif line == '&\n': # save file
            np.savetxt(os.path.join(fname,
                           str(int(ell[mode]))+'_'+str(int(nnn[mode]))),
                       np.array(kernel))
            kernel = []
            mode += 1
        else:
            nums = line.split(' ')
            while '' in nums:
                nums.remove('')
            x = float(nums[0])
            K = float(nums[1].split('\n')[0])
            kernel += [[x, K]]

