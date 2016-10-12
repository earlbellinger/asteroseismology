import numpy as np
from scipy.integrate import cumtrapz, simps

np.seterr(divide='ignore')

def integrate(y, x):
    return np.hstack((0., cumtrapz(y, x)))

def complement(y, x):
    return np.trapz(y, x) - integrate(y, x)


def load_fgong(filename, N=16):
    '''Given an FGONG file, returns a Python dictionary containing
    NumPy arrays that correspond to the structures in the
    specification of the FGONG format:
    
    https://www.astro.up.pt/corot/ntools/docs/CoRoT_ESTA_Files.pdf
    
    That is, the dictionary has arrays indexed by 'nn', 'iconst',
    'ivar', 'ivers', 'glob' and 'var'.'''
    f = open(filename, 'r')
    
    fgong = {'header':[]}
    for i in range(4):
        fgong['header'].append(f.readline())
    
    tmp = [int(i) for i in f.readline().split()]
    fgong['nn'] = tmp[0]
    fgong['iconst'] = tmp[1]
    fgong['ivars'] = tmp[2]
    fgong['ivers'] = tmp[3]
    
    lines = f.readlines()
    tmp = []
    for line in lines:
        for i in range(len(line)//N):
            s = line[i*N:i*N+N]
            if s[-9:] == '-Infinity':
                s = '-Inf'
            elif s[-9:] == ' Infinity':
                s = 'Inf'
            elif s[-4].lower() != 'e':
                s = s[:-4] + 'e' + s[-4:]
            
            tmp.append(float(s))
    
    fgong['glob'] = np.array(tmp[:fgong['iconst']])
    fgong['var'] = np.array(tmp[fgong['iconst']:]).reshape((-1,fgong['ivars']))
    
    f.close()
    
    return fgong



def c2_rho(eig, fgong):
    P = fgong['var'][::-1,3]                 # pressure
    rho = fgong['var'][::-1,4]               # density
    Gamma1 = fgong['var'][::-1,9]            # first adiabatic index
    cs2 = Gamma1*P/rho                       # square of the sound speed
    return (cs2, rho)


def kernel(var1, var2, ell, nu, eig, fgong):
    """Returns a structural kernel.  I have tried to make this as
    notationally similar to Gough & Thompson (1991) as possible.

    Parameters
    ----------
    var1: str
        The first variable of the kernel pair. 
    var2: str
        The second variable of the kernel pair. 
    ell: int
        The angular degree of the mode.
    nu: float
        The cyclic frequency of the mode.
    eig: np.array, shape(N,7)
        Eigenfrequency data for the mode, as produced by ADIPLS.
    fgong: dict
        Stellar model data in a dictionary, as per load_fgong()
        above.

    Returns
    -------
    kernel: np.array, length N
        The density or sound speed structure kernel.
    """
    omega = 2.*np.pi*nu                      # convert to cyclic frequency
    G = 6.672e-8                             # gravitational constant
    L2 = ell*(ell+1)
    L = np.sqrt(L2)
    M, R = fgong['glob'][:2]                 # mass and radius from FGONG
    sigma = np.sqrt(R**3/G/M)*omega          # dimensionless frequency

    ## unpack fgong file (c.f. page 6 section 2.1.3 of above doc)
    r = fgong['var'][::-1,0]                 # radial co-ordinate
    m = M*np.exp(fgong['var'][::-1,1])       # mass co-ordinate
    P = fgong['var'][::-1,3]                 # pressure
    rho = fgong['var'][::-1,4]               # density
    Gamma1 = fgong['var'][::-1,9]            # first adiabatic index
    cs2 = Gamma1*P/rho                       # square of the sound speed
    
    ## equilibrium model (c.f. ADIPLS - The Aarhus adi. osc. pack., section 2.1)
    A = fgong['var'][::-1,14] # 1/Gamma_1 (dln p / dln r) - (dln rho / dln r)
    A1 = (m/M)/(r/R)**3                      # fractional volume
    Vg = (G*m*rho)/(Gamma1*P*r)
    drho_dr = -(Vg+A)*rho/r                  # density gradient
    
    ## unpack eigenfunction (c.f. page 7 of Notes on adi. osc. prog.)
    x = eig[:,0]                             # dimensionless radius (i.e. r/R)
    y1 = eig[:,1]                            # xi_r / R
    y2 = eig[:,2]                            # l(l+1)/R * xi_h
    y3 = eig[:,3]                            # -x * Phi' / (g * r)
  # y4 = eig[:,4]                            # x^2 * d/dx (y_3 / x)
    
    xi_r = y1*R                              # radial component of eigenfunction
    
    dxi_r_dr = np.hstack((0., np.diff(xi_r)/np.diff(r)))
    
    # chi is the "dilatation"
    if ell == 0:
        xi_h = 0.*xi_r 
      # eta = G*m/r**3/omega**2
        chi = Vg/x*(y1-sigma**2/A1/x*y2)
    elif ell > 0:
        xi_h = y2*R/L2
      # eta = L2*G*m/r**3/omega**2
        chi = Vg/x*(y1-y2*sigma**2/(A1*L2)-y3)
    else:
        raise ValueError('ell must be non-negative')
    
    # most stellar models include the central point, at which the
    # numerical value of chi and drho_dr might be buggy.
    chi[0] = 0.
    drho_dr[0] = 0.
    
    S = np.trapz((xi_r**2+L2*xi_h**2)*rho*r**2, r)
    
    # c.f. Gough & Thompson 1991 equation (60)
    K_c_rho = rho*cs2*chi**2*r**2 / (S*omega**2)
    
    # d(c^2)/c^2 = 2 * d(c)/c
    # so K_c2_rho = K_c_rho / 2
    K_c2_rho = K_c_rho / 2.
    
    # first compute the huge bracketed terms 
    # in last two lines of equation (61) 
    K_rho_c = (ell+1.)/r**ell*(xi_r-ell*xi_h) \
        *integrate((rho*chi+xi_r*drho_dr)*r**(ell+2.), r) \
        - ell*r**(ell+1.)*(xi_r+(ell+1.)*xi_h) \
        *complement((rho*chi+xi_r*drho_dr)*r**(1.-ell), r)
    # then combine it with the rest
    K_rho_c = -0.5*(xi_r**2+L2*xi_h**2)*rho*omega**2*r**2 \
        +0.5*rho*cs2*chi**2*r**2 \
        -G*m*(chi*rho+0.5*xi_r*drho_dr)*xi_r \
        -4*np.pi*G*rho*r**2*complement((chi*rho+0.5*xi_r*drho_dr)*xi_r, r) \
        +G*m*rho*xi_r*dxi_r_dr \
        +0.5*G*(m*drho_dr+4.*np.pi*r**2*rho**2)*xi_r**2 \
        -4*np.pi*G*rho/(2.*ell+1.)*K_rho_c
    K_rho_c = K_rho_c / (S*omega**2)
    
    if var1 == 'Gamma1' and var2 == 'rho':
        # K_Gamma1_rho = K_c2_rho (c.f. Basu, Studying Stars 2011, eq 3.20)
        # K_rho_Gamma1 
        # dP/dr = 
        # K_rho_Gamma1 = 
        return K_c2_rho, K_rho_Gamma1
        
    elif var1 == 'c' and var2 == 'rho':
        return K_c_rho, K_rho_c
        
    elif var1 == 'c^2' and var2 == 'rho':
        return K_c2_rho, K_rho_c
        
    else:
        raise ValueError('variable pair not found')

