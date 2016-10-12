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
    Y = 1 - fgong['var'][::-1,5] - fgong['var'][::-1,16] # helium abundance
    
    # Gamma_1,Y = ( partial ln Gamma_1 / partial ln Y ) _ {P, rho} etc
    Gamma_1rho = fgong['var'][::-1,25]
    Gamma_1p = fgong['var'][::-1,26]
    Gamma_1Y = fgong['var'][::-1,27]
    
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
    
    
    ### Calculate c, rho pair
    # c.f. Gough & Thompson 1991 equation (60)
    K_c_rho = rho*cs2*chi**2*r**2 / (S*omega**2)
    
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
    
    
    ### Calculate c^2, rho pair
    # d(c^2)/c^2 = 2 * d(c)/c
    # so K_c2_rho = K_c_rho / 2
    K_c2_rho = K_c_rho / 2.
    K_rho_c2 = K_rho_c
    
    
    ### Calculate Gamma_1, rho pair
    # K_Gamma1_rho = K_c2_rho (c.f. Basu, Studying Stars 2011, eq 3.20)
    # K_rho_Gamma1 (c.f. InversionKit v. 2.2, eq. 105)
    K_rho_Gamma1 = K_rho_c2
    integrand = Gamma1*chi**2*r**2 / (2*S*omega**2)
    first_int = integrate(integrand , r)
    second_int = complement(4*np.pi*G*rho/r**2*integrate(integrand, r), r)
    K_Gamma1_rho = K_rho_c2 - K_c2_rho + G*m*rho/r**2 * \
        first_int + rho*r**2 * second_int
    
    
    ### Calculate u, Y pair
    # K_Y_u = Gamma_1,Y * K_Gamma1_rho (c.f. Basu & Chaplin 2016, eq 10.40)
    # K_u_Y = P * d(P^-1 psi)/dr + Gamma_1,p * K_Gamma1_rho
    # where psi is obtained by solving
    # d[psi(r)]/dr - 4 pi G rho r^2 int_r^R rho/(r^2 P) psi dr = -z(r)
    # z(r) = K_rho_Gamma1 + (Gamma_1,p + Gamma_1,rho) K_Gamma1_rho
    K_Y_u = Gamma_1Y * K_Gamma1_rho

    z = K_rho_Gamma1 + ( Gamma_1p + Gamma_1rho ) * K_Gamma1_rho
    
    rho_spl = UnivariateSpline(r, rho)
    P_spl = UnivariateSpline(r, P)
    z_spl = UnivariateSpline(r, z)
    
    drho = rho_spl.derivative()
    dz = z_spl.derivative()
    
    def bvp(t, psi):
        # psi'' + a(r) psi' + b(r) psi = f(r)
        # where
        # a(r) =  P/rho (r rho' - 2 rho)
        # b(r) = -P/rho [ 1 / ( 4 pi G r^2 rho ) ]
        # f(r) =  P/rho ( z r rho' - rho r z' + 2 z rho )
        # rearrange into y_1' = y_2
        #                y_2' = f(t) - a(t) y_2 - b(t) psi 
        
        # remesh onto variable x from r
        rhot = rho_spl(t) 
        Pt = P_spl(t) 
        zt = z_spl(t) 
        
        # differentiate with respect to x
        drhot = drho(t) 
        dzt = dz(t) 
        
        # calculate variables on mesh x
        #Prhot = Pt / rhot
        #at = Prhot * ( t * rhot - 2 * rhot )
        #bt = - Prhot * ( 1 / ( 4 * np.pi * G * t**2 * rhot ) )
        #ft = Prhot * ( zt * t * drhot - \
        #               rhot * x * dzt + \
        #               2 * zt * rhot )
        
        drhorho = drhot / rhot
        
        at = drhorho - 2 / t
        bt = 4 * np.pi * G * t * rhot**2 / Pt
        ft = zt * (2/t + drhorho) - dzt
        
        return np.vstack(( psi[1], 
                           ft - at * psi[1] - bt * psi[0] ))
    
    result = solve_bvp(bvp, 
        lambda ya, yb: np.array([ ya[0]-1, yb[0]-1 ]), 
        r, 
        np.array([ np.ones(len(r)), np.ones(len(r)) ]), 
        max_nodes=1000000)#, tol=1e-5)
    
    
    
    
    
    #def dpsi_dr(psi, s):
    #    return 4*np.pi*G * np.interp(s, r, rho) * s**2 \
    #        * trapz(rho[r>=s] / (r[r>=s]**2 * P[r>=s]) * psi, r[r>=s]) \
    #        - ( np.interp(s, r, K_rho_Gamma1) \
    #            + ( np.interp(s, r, Gamma_1p) + np.interp(s, r, Gamma_1rho) ) \
    #              * np.interp(s, r, K_Gamma1_rho) \
    #          )
    #psi = odeint(dpsi_dr, 1.0, r).T[0]
    #print(psi)
    
    
    
    #def dpsi_dr(r2, psi):
    #    print('r2 == r', np.all(r2 == r))
    #    psi = psi[0]
    #    print('len(psi) == len(r)', len(psi) == len(r))
    #    return np.vstack((4*np.pi*G*rho*r**2 \
    #        * complement(rho/(r**2*P)*psi, r) \
    #        + ( K_rho_Gamma1 + ( Gamma_1p + Gamma_1rho ) * K_Gamma1_rho ))).T
    
    #def dpsi_dr(r2, psi):
    #    psi = psi[0]
    #    return np.vstack((4*np.pi*G*np.interp(r2, r, rho)*r2**2 \
    #        * complement(np.interp(r2, r, rho)/(r2**2*np.interp(r2,r,P))*psi, r2) \
    #        + ( np.interp(r2, r, K_rho_Gamma1) + ( np.interp(r2, r, Gamma_1p) \
    #        + np.interp(r2, r, Gamma_1rho) ) * np.interp(r2, r, K_Gamma1_rho) ))).T
    
    ##psi = np.zeros(len(r))
    ##psi[0] = r[0] * f(0, 0)
    ##for ii in range(1, len(r)):
    ##    h = r[ii] - r[ii-1]
    ##    psi[ii] = psi[ii-1] + h * f(psi[ii-1], ii-1)
    ##print(psi)
    
    #def resid(psi):
    #    return np.gradient(psi, r) - f(psi)
    
    
    
    #result = solve_bvp(dpsi_dr, 
    #    lambda ya,yb: np.array([ (ya[0]-1)**2 + (yb[0]-1)**2 ]),
    #    r, [np.ones(len(r))], max_nodes=100000, tol=1e-9)
    
    
    # dPsi1_dr = 4 pi G r^2 rho Psi_2 + 
    #            [ K_rho_Gamma1 + ( Gamma_1P + Gamma_1rho ) * K_Gamma1_rho ]
    # dPsi2_dr = -rho/(r^2 P) Psi_1
    
    #def bvp(t, psi):
    #    rhot = rho_spl(t) 
    #    psi1 = 4 * np.pi * G * rhot * t**2 * psi[1] + z_spl(t)
    #    psi2 = - rhot / (t * P_spl(t)) * psi[0]
    #    return np.vstack(( psi1, psi2 ))
    
    #result = solve_bvp(bvp, 
    #    lambda ya, yb: np.array([ ya[0]-1, yb[0]-1 ]), 
    #    r, 
    #    np.array([ np.ones(len(r)), np.zeros(len(r)) ]), 
    #    max_nodes=100000000, tol=1e-4)
    
#    def bvp(t, psi):
#        rhot = rho_spl(t) 
#        dpsi = 4 * np.pi * G * rhot * t**2 \
#            * complement(rhot / ( t**2 * P_spl(t) ) * psi[0], t) - z_spl(t)
#        print('dpsi', dpsi)
#        return np.vstack(( dpsi )).T
    
#    result = solve_bvp(bvp, 
#        lambda ya, yb: np.array([ yb[0]-1 ]), 
#        r, 
#        np.array([ np.ones(len(x)) ]), 
#        max_nodes=100000000, tol=1e-3) 
    
    #result = solve_bvp(dPsi2_dr2, 
    #    lambda ya, yb: np.array([ (ya[0]-1)**2, (yb[0]-1)**2 ]),
    #    r, [ np.ones(len(r)), np.zeros(len(r)) ], 
    #    max_nodes=1e6)#, tol=1e-5)
    
    print(result)
    psi = result.sol(r)[0] #interp1d(result['x'], result['y'][0])(r)
    print('psi', psi)
    
    #result = least_squares(resid, np.ones(len(r)))
    #print(result)
    #psi = result['x']
    #psi = minimize(sqdiff, np.ones(len(r)), method='Nelder-Mead')['x']
    
    K_u_Y = P * UnivariateSpline(P**-1 * psi, r).derivative()(r) \
        + Gamma_1p * K_Gamma1_rho
    
    
    return { 
        ('c',      'rho'): (K_c_rho,      K_rho_c2),
        ('c2',     'rho'): (K_c2_rho,     K_rho_c2),
        ('Gamma1', 'rho'): (K_Gamma1_rho, K_rho_Gamma1),
        ('u',      'Y'):   (K_u_Y,        K_Y_u),
        ('psi', 'psi'):    (psi,          psi)
    }

