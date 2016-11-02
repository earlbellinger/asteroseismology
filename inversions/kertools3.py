import numpy as np
from scipy.integrate import cumtrapz, simps, trapz, odeint, solve_bvp
from scipy.optimize import least_squares #minimize
from scipy.interpolate import interp1d, UnivariateSpline

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

def kernel(ell, nu, eig, fgong):
    """Returns a dict of structural kernels.  I have tried to make this as
    notationally similar to Gough & Thompson (1991) as possible.
    
    Parameters
    ----------
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
    
    G = 6.67232e-8                           # gravitational constant
    M, R = fgong['glob'][:2]                 # mass and radius from FGONG
    
    ## unpack eigenfunction (c.f. page 7 of Notes on adi. osc. prog.)
    x = eig[:,0]                             # dimensionless radius (i.e. r/R)
    y1 = eig[:,1]                            # xi_r / R
    y2 = eig[:,2]                            # l(l+1)/R * xi_h
    y3 = eig[:,3]                            # -x * Phi' / (g * r)
  # y4 = eig[:,4]                            # x^2 * d/dx (y_3 / x)
    
    xi_r = y1*R                              # radial component of eigenfunction
    xi_h = y2*R/L2 if ell > 0 else 0. * xi_r # horiz. component of eigenfunction
    
    ## unpack fgong file (c.f. page 6 section 2.1.3 of above doc)
    r = fgong['var'][::-1,0]                 # radial co-ordinate
    m = M*np.exp(fgong['var'][::-1,1])       # mass co-ordinate
    P_0 = fgong['var'][::-1,3]               # pressure
    rho_0 = fgong['var'][::-1,4]             # density
    Gamma1 = fgong['var'][::-1,9]            # first adiabatic index
    c02 = Gamma1*P_0/rho_0                   # square of the sound speed
    Y = 1 - fgong['var'][::-1,5] - fgong['var'][::-1,16] # helium abundance
    # Gamma_1,Y = ( partial ln Gamma_1 / partial ln Y ) _ {P, rho} etc
    #Gamma_1rho = fgong['var'][::-1,25]
    #Gamma_1p = fgong['var'][::-1,26]
    #Gamma_1Y = fgong['var'][::-1,27]
    
    ## equilibrium model (c.f. ADIPLS - The Aarhus adi. osc. pack., section 2.1)
    #A = fgong['var'][::-1,14]    # 1/Gamma_1 (dln p / dln r) - (dln rho / dln r)
    #A1 = (m/M)/(r/R)**3                      # fractional volume
    #Vg = (G*m*rho_0)/(Gamma1*P_0*r)
    #drho_dr = -(Vg+A)*rho_0/r                # density gradient
    
    omega = 2.*np.pi*nu                      # convert to cyclic frequency
    sigma = np.sqrt(R**3/G/M)*omega          # dimensionless frequency
    
    L2 = ell*(ell+1)
    L = np.sqrt(L2)
    
    # eta = m_0 / r**3
    # m_0 = 4*pi*int_s=0^r p_0(s) s^2 ds
    m_0 = 4 * np.pi * integrate( P_0 * r**2 , r )
    eta = m_0 / r**3
    
    # chi ("dilatation") = (d(xi)/dr) + 2(xi)/r - l(l+1)*eta/r
    dxi_r_dr = np.hstack((0., np.diff(xi_r)/np.diff(r)))
    chi = dxi_r_dr + 2*xi_r/r - L2*eta/r
    chi[0] = 0.
    
    # rho = -(d(rho_0)/dr)*xi_r - p_0*chi)
    drho_0_dr = np.hstack((0., np.diff(rho_0)/np.diff(r)))
    rho = -drho_0_dr * xi_r - rho_0 * chi
    
    # psi = -4*pi*G / (2*l+1) * [ int_s=0^r p(s) s^(l+2)/r^(l+1) ds +
    #                         l * int_s=r^R rho(s) r^l / s^(l-1) ds ]
    psi = -4*np.pi*G / (2*ell+1) * \
           ( integrate( rho * r**(ell+2) / r[-1]**(ell+2), r ) + \
        L * complement( rho * r[0]**(ell) / r**(ell-1), r ) )
    dpsi_dr = np.hstack((0., np.diff(psi)/np.diff(r)))
    
    g_0 = G * m_0 / r**2
    
    # I = int_r=0^R rho_0(r) (xi_r**2 + l(l+1) eta^2) r^2 dr
    I = trapz( rho_0 * ( xi_r**2 + L2 * eta**2 ) * r**2 , r )
    
    #S = np.trapz((xi_r**2+L2*xi_h**2)*rho*x**2, x)
    
    ### Calculate c, rho pair
    K_c2_rho = rho_0 * c02 * chi**2 * r**2 / ( 2 * I * sigma**2 )
    K_rho_c2 = K_c2_rho
    
    return { 
        #('c',      'rho'): (K_c_rho,      K_rho_c2),
        ('c2',     'rho'): (K_c2_rho,     K_rho_c2),
        #('Gamma1', 'rho'): (K_Gamma1_rho, K_rho_Gamma1)
        #('u',      'Y'):   (K_u_Y,        K_Y_u),
        #('psi', 'psi'):    (psi,          psi)
    }

