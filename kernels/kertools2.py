import numpy as np
from numpy.linalg import inv
from scipy.integrate import cumtrapz, simps, trapz, odeint, solve_bvp
from scipy.optimize import least_squares #minimize
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
#from scipy.special import sph_harm
#import scikits.bvp_solver

np.seterr(divide='ignore')

def integrate(y, x):
    return np.hstack((0., cumtrapz(y, x)))

def complement(y, x):
    return np.trapz(y, x) - integrate(y, x)

def derivative(y, x, k=2):
    return InterpolatedUnivariateSpline(x, y, k=k).derivative()(x)

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
G = 6.67232e-8                           # gravitational constant
L2 = ell*(ell+1)
L = np.sqrt(L2)
M, R = fgong['glob'][:2]                 # mass and radius from FGONG
sigma = np.sqrt(R**3/G/M)*omega          # dimensionless frequency

## unpack fgong file (c.f. page 6 section 2.1.3 of above doc)
r = fgong['var'][::-1,0]                 # radial co-ordinate
log_q = fgong['var'][::-1,1]             # dimensionless mass
m = M*np.exp(log_q)                      # mass co-ordinate
P = fgong['var'][::-1,3]                 # pressure
rho = fgong['var'][::-1,4]               # density
Gamma1 = fgong['var'][::-1,9]            # first adiabatic index
cs2 = Gamma1*P/rho                       # square of the sound speed
Y = 1 - fgong['var'][::-1,5] - fgong['var'][::-1,16] # helium abundance

# Gamma_1,Y = ( partial ln Gamma_1 / partial ln Y ) _ {P, rho} etc
Gamma_1rho = fgong['var'][::-1,25]
Gamma_1P = fgong['var'][::-1,26]
Gamma_1Y = fgong['var'][::-1,27]

## unpack eigenfunction (c.f. page 7 of Notes on adi. osc. prog.)
x = eig[:,0]                             # dimensionless radius (i.e. r/R)
y1 = eig[:,1]                            # xi_r / R
y2 = eig[:,2]                            # l(l+1)/R * xi_h
y3 = eig[:,3]                            # -x * Phi' / (g * r)
# y4 = eig[:,4]                            # x^2 * d/dx (y_3 / x)

## equilibrium model (c.f. ADIPLS - The Aarhus adi. osc. pack., section 2.1)
#A = fgong['var'][::-1,14] # 1/Gamma_1 (dln p / dln r) - (dln rho / dln r)
#A1 = (np.exp(log_q))/(x)**3              # fractional volume
#Vg = (G*m*rho)/(Gamma1*P*r)
#drho_dr = -(Vg+A)*rho/r                  # density gradient

xi_r = y1*R                              # radial component of eigenfunction
if ell != 1: xi_r[0] = 0.
xi_h = y2*R/L2 if ell > 0 else xi_r*0.   # horiz. component of eigenfunction

#dxi_r_dr = np.hstack((0., np.diff(y1)/np.diff(x)))
#drho_dr = np.hstack((0., np.diff(rho)/np.diff(r)))
drho_dr = derivative(rho, r)
dxi_r_dr = derivative(y1, x)

#print(dxi_r_dr)

# most stellar models include the central point, at which the
# numerical value of chi and drho_dr might be buggy.
dxi_r_dr[0] = 0.
drho_dr[0] = 0.
    
    # chi is the "dilatation"
    #if ell == 0:
    #    xi_h = 0.*xi_r 
    #    eta = G*m/r**3/omega**2
    #   chi = np.hstack((0, Vg[1:]/x[1:]*(y1[1:]-sigma**2/A1[1:]/x[1:]*y2[1:])))
    #elif ell > 0:
    #    xi_h = y2*R/L2
    #    eta = L2*G*m/r**3/omega**2
    #    chi = np.hstack((0, Vg[1:]/x[1:]*(y1[1:]-y2[1:]*sigma**2/(A1[1:]*L2)-y3[1:])))
    #else:
    #    raise ValueError('ell must be non-negative')
    
    
    # eta = m_0 / r**3
    #eta = xi_h #m_0 / r**3
    
    # chi ("dilatation") = (d(xi)/dr) + 2(xi)/r - l(l+1)*eta/r
chi = dxi_r_dr + np.hstack((0., 2*y1[1:]/x[1:])) #2*xi_r/r
if ell > 0:
    chi -= np.hstack((0., y2[1:]/x[1:]))
    #chi[0] = 3*dxi_r_dr[0] # Daniel Reese does this 


S = 1 / ( omega**2 * np.trapz((xi_r**2 + L2*xi_h**2) * rho * r**2, r) )
#if ell > 0:
#    S = np.trapz((y1**2 + L2*(y2/L2)**2) * rho * x**2, x)
#else:
#    S = np.trapz(y1**2 * rho * x**2, x)


#S = np.trapz((xi_r**2 + L2*eta**2)*rho*r**2, r)
#S = np.trapz((y1**2 + L2*eta**2)*rho*x**2, r)


### Calculate c, rho pair
# c.f. Gough & Thompson 1991 equation (60)
#K_c_rho = rho * cs2 * chi**2 * r**2 / (S * omega**2)
K_c_rho = rho * cs2 * chi**2 * r**2 * S #/ S #/ (S * sigma**2)
K_c_rho[0] = 0.

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
K_rho_c = K_rho_c * S #(S*omega**2)
K_rho_c[0] = 0.
    
    ### Reese implementation
    # m_0 = 4*pi*int_s=0^r rho_0(s) s^2 ds
    #m_0 = 4 * np.pi * integrate(rho * r**2, r)
    #g = G * m_0 / r**2 #4.*np.pi/3.*G*rho*r
    #rho2 = - drho_dr * xi_r - rho * chi
    
    #psi = np.hstack((0, np.cumsum([
    #        np.trapz(rho2[:i] * r[:i]**(ell+2) / r[i]**(ell+1),  r[:i]) + \
    #        np.trapz(rho2[i:] * r[i]**ell      / r[i:]**(ell-1), r[i:])
    #    for i in range(1, len(r))])))
    
    #dpsi_dr = np.hstack((0, np.cumsum([(ell+1) * \
    #        np.trapz(rho2[:i] * r[:i]**(ell+2) / r[i]**(ell+2), r[:i]) + \
    #        ell * np.trapz(rho2[i:] * r[i]**(ell-1) / r[i:]**(ell-1), r[i:])
    #    for i in range(1, len(r))])))
    
    
    #psi = integrate(rho*r**(ell+2), r)/(r**(ell+1)) + \
    #    r**ell*complement(rho/r**(ell-1), r)
    
    #dpsi_dr = -(ell+1)*integrate(rho*r**(ell+2), r)/(r**(ell+1)) + \
    #    ell*complement(rho/r**(ell-1), r)*r**(ell-1)
            
    #c_psi = -4*np.pi*G / (2*ell+1)
    #psi = psi * c_psi
    #dpsi_dr = dpsi_dr * c_psi
    
    #K_rho_c = rho * r**2 * S/2 * \
    #    ( cs2 * chi**2 - omega**2 * ( xi_r**2 + L2 * xi_h**2 ) - \
    #      2. * g * xi_r * chi - \
    #      4. * np.pi * G * complement(2*rho*xi_r*chi + drho_dr * xi_r**2, r) + \
    #      2. * g * xi_r * dxi_r_dr + \
    #      4. * np.pi * G * rho * xi_r**2 + \
    #      2. * (xi_r * dpsi_dr + L2*xi_h*psi / r) )
    #K_rho_c[0] = 0.
    
    #print(K_rho_c)
    
    ### Calculate c^2, rho pair
    # d(c^2)/c^2 = 2 * d(c)/c
    # so K_c2_rho = K_c_rho / 2
K_c2_rho = K_c_rho / 2.
K_rho_c2 = K_rho_c


### Calculate Gamma_1, rho pair
# K_Gamma1_rho = K_c2_rho (c.f. Basu, Studying Stars 2011, eq 3.20)
# K_rho_Gamma1 (c.f. InversionKit v. 2.2, eq. 105)
K_Gamma1_rho = K_c2_rho
#integrand = Gamma1*chi**2*r**2 / (2*S*omega**2)
first_int = integrate(K_c2_rho/P , r)
second_int = complement(4*np.pi*G*rho/r**2*first_int, r)
K_rho_Gamma1 = K_rho_c2 - K_c2_rho + G*m*rho/r**2 * \
    first_int + rho*r**2 * second_int
K_rho_Gamma1[0] = 0.


### Calculate u, Y pair
K_Y_u = Gamma_1Y * K_Gamma1_rho
    
    # Kosovichev 1999 eq 40
    #V = G * m * rho / (r * P)
    #U = 4 * np.pi * rho * r**3 / m
    
    #V = -derivative(np.log(P), np.log(x))
    #U = derivative(np.log(m), np.log(x))
    #U[0] = 0
    #V[0] = 0
    
    #m_0 = np.cumsum(m) #4 * np.pi * integrate(rho * r**2, r)
    #V = -derivative(np.log(P), t)
    #U = derivative(np.log(m), t)
    #U[0] = 0.
    
g = G*m/r**2
gp = derivative(g, r)
F = 1/r * (2 * xi_r - L2 * xi_h)
Kp = (dxi_r_dr + F)**2
E = omega**2 * np.trapz( rho*r**2 * (xi_r + L2 * xi_h**2), r)
C = -1/r**2 * integrate(Gamma1 * Kp * r**2, r)
Si = -2 * complement( rho * xi_r * F, r)
S1i = complement(rho * C, r)
K_rho_Gamma2 = rho * r**2 / E * ( -omega**2 * (xi_r**2 + L2 * xi_h**2) +\
    2 * xi_r * (gp + 4*np.pi*G*rho*xi_r - F*g) - C*g + 4*np.pi*G*(Si - S1i)) 

q = np.exp(log_q)
t = r/R
V = -derivative(np.log(P), np.log(t), k=1)
U = derivative(log_q, np.log(t), k=1)

#print(np.vstack((V, G * m * rho / (r * P), U, 4 * np.pi * rho * r**3 / m)).T)
#V = G * m * rho / (r * P)
#U = 4 * np.pi * rho * r**3 / m
m_0 = 4*np.pi * integrate(rho * r**2, r)
#V = G * m * rho / (t * P)
#U = 4 * np.pi * rho * t**3 / m
V[0] = 0.
U[0] = 0.
#V = G * m * rho / (r * P)
#U = 4 * np.pi * rho * r**3 / m
V = G * m * rho / (t**2 * P)
U = 4 * np.pi * rho * t**2 / m
V[0] = 0.
U[0] = 0.
    
par1 = K_rho_Gamma1+(Gamma_1rho+Gamma_1P)*K_Gamma1_rho
np.savetxt('PAR', np.hstack((par1, U, V)))
np.savetxt('T', t)
#exit(1)

a = (2/t + drho_dr/rho) #(t*drho_dr + 2*rho)/(t*rho)
a[0]=0
b = 4*np.pi*G*rho**2/P
F = (Gamma_1P + Gamma_1rho) * K_Gamma1_rho + K_rho_Gamma1
c = t**2*rho*derivative(F/(t**2*rho), t)
c[0]=0
np.savetxt('PAR', np.hstack((a, b, c)))
np.savetxt('T', t)

#asdf1, asdf2, w_1, w_2, asdf3, asdf4 = np.loadtxt('out').T
asdf1, asdf2, w_1, w_2, asdf3, asdf4 = np.loadtxt('out').T
    
    #def bvp(a, w):
    #    dw1dr = -(K_rho_Gamma1 + (Gamma_1P+Gamma_1rho)*K_Gamma1_rho) - U * w[1]
    #    dw2dr = U * w[1] - V * w[0]
    #    return np.array((dw1dr, dw2dr))
    
    #def bcs(ra, rb):
    #    return (np.array([ra[0]]), np.array([rb[1]]))
    
    #problem = scikits.bvp_solver.ProblemDefinition(
    #    num_ODE=2, 
    #    num_parameters=0, 
    #    num_left_boundary_conditions=1, 
    #    boundary_points=(0, R), 
    #    function=bvp, 
    #    boundary_conditions=bcs)
    #solution = scikits.bvp_solver.solve(problem,
    #    solution_guess=np.array((len(r), 2)))
    
    # Kosovichev 1999 eqs 43, 44, 45
#   A_1 = np.array( [ [  V, -V ], [          0,        -U ] ] )
#   B_1 = np.array( [ [ -V,  0 ], [          U,         0 ] ] )
#   B_2 = np.array( [ [  0,  0 ], [          0,         0 ] ] )
#   C_1 = np.array( [ [  1,  0 ], [ Gamma_1P  ,         0 ] ] )
#   D_1 = np.array( [ [  1,  0 ], [ Gamma_1rho,         1 ] ] )
#   D_2 = np.array( [ [ -1,  0 ], [          0,  Gamma_1Y ] ] ) 

# Kosovichev 1999 eq 48
#   A = A_1 + np.matmul(np.matmul(B_1, inv(D_1)), C_1)
#   B = np.matmul(np.matmul(B_1, inv(D_1)), D_2) + B_2
#   C = np.matmul(inv(D_1), C_1)
#   D = np.matmul(inv(D_1), D_2)
    
    # now we want to solve for K^(2) = (K_{u,Y}, K_{Y,u}) 
    # we have that (eqs 26, 27)
    #   dy/dx = Ay + Bz_2
    #   z_1 = Cy + Dz_2
    # with (underneath eq 42)
    #   z_2 = (drho/rho, dGamma1/Gamma1)
    # which means (eqs 33, 34)
    #   dw/dx = -A^T w - C^T K^(1)
    #   K^(2) = D^T K^(1) + B^T w
    # with boundary conditions (eq 28)
    #   w * y = 0   at both  r=0 and r=R
    # where (underneath eq 42)
    #   y = (dp/p, dm/m)
    # so w.l.o.g. we have
    #   w_1 = 0 at r=0
    #   w_2 = 0 at r=R
    
    #z_1 = np.array([ K_rho_Gamma1, K_Gamma1_rho ])
    #dw_dx = -A.T * w - np.matmul(C.T, K_1)
    #K_2 = np.matmul(D.T, K_1) + B.T * w 
    
    # set up call to FDM
    #dw_dx = 
    #N = len(r) # number of mesh points
    #M = 2 # number of first order ODEs
    #ML = 1 # # number of boundary conditions at the first boundary
    #PAR = # real array passed to EQN and BCS for calculating info
    #X = np.zeros(( M, N )) # initial guess and eventually computed solution
    #XC = np.ones(( M, N )) # solution after applying deferred correction - not used
    #T = r
    #EQN = 

U_spl = UnivariateSpline(t, U, k=1)
V_spl = UnivariateSpline(t, V, k=1)
K_spl = UnivariateSpline(t, 
    (K_rho_Gamma1+(Gamma_1rho+Gamma_1P)*K_Gamma1_rho), k=1)

KrG_spl = UnivariateSpline(t, K_rho_Gamma1, k=1)
Gr_spl = UnivariateSpline(t, Gamma_1rho, k=1)
GP_spl = UnivariateSpline(t, Gamma_1P, k=1)
KGr_spl = UnivariateSpline(t, K_Gamma1_rho, k=1)

def bvp(t, w):
    U_t = U_spl(t)
    V_t = V_spl(t)
    K_t = K_spl(t)
    
    dw1dr = -U_t*w[1] - KrG_spl(t) * ( Gr_spl(t) + GP_spl(t) ) * KGr_spl(t) / t
    dw2dr = V_t*w[0] + U_t*w[1]
    return np.array((dw1dr, dw2dr))

result = solve_bvp(bvp, 
    lambda ya, yb: np.array([ ya[0], yb[1] ]), 
    t, np.array([ np.zeros(len(t)), np.zeros(len(t)) ]), 
    max_nodes=10000000, 
    tol=1e-5)

w_1 = result.sol(t)[0]
w_2 = result.sol(t)[1]






def bvp(t, w):
    rho_t = rho_spl(t)
    P_t = P_spl(t)
    drho_t = drho_spl(t)
    #dF_t = dF_spl(t)
    GP = GP_spl(t)
    Gr = Gr_spl(t)
    KGr = KGr_spl(t)
    KrG = KrG_spl(t)
    
    F = (GP + Gr) * KGr + KrG
    
    dPsi_dr = F - 4*np.pi*G*t**2*rho_t * complement(rho_t/(t**2*P_t) * w[0], t)
    return np.array([dPsi_dr])

result = solve_bvp(bvp, 
    lambda ya, yb: np.array([ yb[0]-1 ]), 
    t, np.array([ np.ones(len(t)) ]), 
    max_nodes=100000000, 
    tol=1e-5)

w_1 = result.sol(t)[0]
w_2 = result.sol(t)[1]






def bvp(t, w):
    rho_t = rho_spl(t)
    P_t = P_spl(t)
    #dF_t = dF_spl(t)
    GP = GP_spl(t)
    Gr = Gr_spl(t)
    KGr = KGr_spl(t)
    KrG = KrG_spl(t)
    
    F = (GP + Gr) * KGr + KrG
    
    dPsi_dr = t**2 * rho_t * w[1] + F
    dphi_dr = -4*np.pi*G*rho_t / (t**2 * P_t) * w[0]
    return np.array([dPsi_dr, dphi_dr])

result = solve_bvp(bvp, 
    lambda ya, yb: np.array([ yb[0]-1 ]), 
    t, np.array([ np.ones(len(t)) ]), 
    max_nodes=100000000, 
    tol=1e-5)

w_1 = result.sol(t)[0]
w_2 = result.sol(t)[1]


rho_spl = UnivariateSpline(t, rho, k=1)
P_spl = UnivariateSpline(t, P, k=1)
drho_spl = UnivariateSpline(t, drho_dr, k=1)
#dF_spl = UnivariateSpline(t, derivative(((Gamma_1P+Gamma_1rho)*K_Gamma1_rho + K_rho_Gamma1)/(t**2*rho), r), k=1)
KrG_spl = UnivariateSpline(t, K_rho_Gamma1, k=1)
Gr_spl = UnivariateSpline(t, Gamma_1rho, k=1)
GP_spl = UnivariateSpline(t, Gamma_1P, k=1)
KGr_spl = UnivariateSpline(t, K_Gamma1_rho, k=1)

def bvp(t, w):
    rho_t = rho_spl(t)
    P_t = P_spl(t)
    drho_t = derivative(rho_t, t)#drho_spl(t)
    #dF_t = dF_spl(t)
    GP = GP_spl(t)
    Gr = Gr_spl(t)
    KGr = KGr_spl(t)
    KrG = KrG_spl(t)
    
    dF = derivative( ((GP+Gr)*KGr+KrG) / (t**2*rho_t) , t)
    #a = (t*drho_t + 2*rho_t) / (t*rho_t)
    #b = 4*np.pi*G*rho_t**2/P_t
    #f = t**2*rho_t*dF
    
    dPsidr = w[1]
    #dzdr = a*w[1] - b*w[0] + f
    dzdr = (2/t+drho_t/rho_t)*w[1] - 4*np.pi*G*rho_t**2/P_t*w[0] + t**2*rho_t*dF
    return np.array((dPsidr, dzdr))

result = solve_bvp(bvp, 
    lambda ya, yb: np.array([ ya[0], yb[0] ]), 
    t, np.array([ np.zeros(len(t)), np.zeros(len(t)) ]), 
    max_nodes=100000000, 
    tol=1e-5)

w_1 = result.sol(t)[0]
w_2 = result.sol(t)[1]

#K_u_Y = Gamma_1rho * K_Gamma1_rho - K_rho_Gamma1 + V * w_1 - U * w_2
#asdf1, asdf2, w_1, w_2, asdf3, asdf4 = np.loadtxt('out').T
#K_u_Y = -K_rho_Gamma1 - Gamma_1rho*K_Gamma1_rho + V*w_1 - U*w_2
#K_u_Y = P * derivative(w_1 / P, r/R) + Gamma_1P * K_Gamma1_rho
K_u_Y = Gamma_1P*K_Gamma1_rho - P * derivative(w_1/P, t)
print(np.vstack(( r, w_1, w_2, R*K_u_Y )).T)
save(x, K_u_Y, K_Y_u, R, 'u', 'Y', ell, n, False, output_dir, True, True)



save(x, (-K_rho_Gamma1) / max(-K_rho_Gamma1), (-Gamma_1rho*K_Gamma1_rho)/max(-Gamma_1rho*K_Gamma1_rho), 1, 'psi', 'psi', ell, n, False, output_dir, True, True)
save(x, (-K_rho_Gamma1 - Gamma_1rho*K_Gamma1_rho) / max(-K_rho_Gamma1 - Gamma_1rho*K_Gamma1_rho), (V*w_1 - U*w_2)/max(V*w_1 - U*w_2), 1, 'psi', 'psi', ell, n, False, output_dir, True, True)
save(x, V*w_1/max(V*w_1), -U*w_2/max(-U*w_2), 1, 'psi', 'psi', ell, n, False, output_dir, True, True)


"""
rho_spl = UnivariateSpline(t, rho, k=1)
P_spl = UnivariateSpline(t, P, k=1)
def bvp(t, w):
    rho_t = rho_spl(t)
    P_t = P_spl(t)
    K_t = K_spl(t)
    dw1dr = 4*np.pi*G*rho_t*t**2*w[1] - K_t
    dw2dr = - rho_t / (t**2 * P_t) * w[0]
    return np.array((dw1dr, dw2dr))

result = solve_bvp(bvp, 
    lambda ya, yb: np.array([ ya[0]-1, ya[1]-1 ]), 
    t, np.array([ np.zeros(len(t)), np.zeros(len(t)) ]), 
    max_nodes=100000000, 
    tol=1e-3)

#w_1 = np.hstack((0, result.sol(t[1:])[0]))
#w_2 = np.hstack((0, result.sol(t[1:])[1]))
w_1 = result.sol(t)[0]
w_2 = result.sol(t)[1]




rho_spl = UnivariateSpline(t, rho, k=3)
P_spl = UnivariateSpline(t, P, k=3)
K_spl = UnivariateSpline(t, (K_rho_Gamma1-(Gamma_1rho+Gamma_1P)*K_Gamma1_rho), k=3)
def bvp(t, w):
    rho_t = rho_spl(t)
    P_t = P_spl(t)
    K_t = K_spl(t)
    return np.array([4*np.pi*G*rho_t*R*t**2* complement(rho_t / (t**2 * P_t) * w[0], t) - K_t])

def bc(ya, yb):
    return np.array([ yb[0]-1 ])

result = solve_bvp(bvp, bc, t, np.ones((1, len(t))),
    max_nodes=100000000, 
    tol=1e-3)

"""

    # K_Y_u = Gamma_1,Y * K_Gamma1_rho (c.f. Basu & Chaplin 2016, eq 10.40)
    # K_u_Y = P * d(P^-1 psi)/dr + Gamma_1,p * K_Gamma1_rho
    # where psi is obtained by solving
    # d[psi(r)]/dr - 4 pi G rho r^2 int_r^R rho/(r^2 P) psi dr = -z(r)
    # z(r) = K_rho_Gamma1 + (Gamma_1,p + Gamma_1,rho) K_Gamma1_rho
#    K_Y_u = Gamma_1Y * K_Gamma1_rho
#
#    z = K_rho_Gamma1 + ( Gamma_1P + Gamma_1rho ) * K_Gamma1_rho
#    
#    rho_spl = UnivariateSpline(r, rho)
#    P_spl = UnivariateSpline(r, P)
#    z_spl = UnivariateSpline(r, z)
#    
#    drho = rho_spl.derivative()
#    dz = z_spl.derivative()
#    
#    def bvp(t, psi):
#        # psi'' + a(r) psi' + b(r) psi = f(r)
#        # where
#        # a(r) =  P/rho (r rho' - 2 rho)
#        # b(r) = -P/rho [ 1 / ( 4 pi G r^2 rho ) ]
#        # f(r) =  P/rho ( z r rho' - rho r z' + 2 z rho )
#        # rearrange into y_1' = y_2
#        #                y_2' = f(t) - a(t) y_2 - b(t) psi 
#        
#        # remesh onto variable x from r
#        rhot = rho_spl(t) 
#        Pt = P_spl(t) 
#        zt = z_spl(t) 
#        
#        # differentiate with respect to x
#        drhot = drho(t) 
#        dzt = dz(t) 
#        
#        # calculate variables on mesh x
#        #Prhot = Pt / rhot
#        #at = Prhot * ( t * rhot - 2 * rhot )
#        #bt = - Prhot * ( 1 / ( 4 * np.pi * G * t**2 * rhot ) )
#        #ft = Prhot * ( zt * t * drhot - \
#        #               rhot * x * dzt + \
#        #               2 * zt * rhot )
#        
#        drhorho = drhot / rhot
#        
#        at = drhorho - 2 / t
#        bt = 4 * np.pi * G * t * rhot**2 / Pt
#        ft = zt * (2/t + drhorho) - dzt
#        
#        return np.vstack(( psi[1], 
#                           ft - at * psi[1] - bt * psi[0] ))
#    
#    result = solve_bvp(bvp, 
#        lambda ya, yb: np.array([ ya[0]-1, yb[0]-1 ]), 
#        r, 
#        np.array([ np.ones(len(r)), np.ones(len(r)) ]), 
#        max_nodes=1000000)#, tol=1e-5)
    
    
    
    
    
    #def dpsi_dr(psi, s):
    #    return 4*np.pi*G * np.interp(s, r, rho) * s**2 \
    #        * trapz(rho[r>=s] / (r[r>=s]**2 * P[r>=s]) * psi, r[r>=s]) \
    #        - ( np.interp(s, r, K_rho_Gamma1) \
    #            + ( np.interp(s, r, Gamma_1P) + np.interp(s, r, Gamma_1rho) ) \
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
    #        + ( K_rho_Gamma1 + ( Gamma_1P + Gamma_1rho ) * K_Gamma1_rho ))).T
    
    #def dpsi_dr(r2, psi):
    #    psi = psi[0]
    #    return np.vstack((4*np.pi*G*np.interp(r2, r, rho)*r2**2 \
    #        * complement(np.interp(r2, r, rho)/(r2**2*np.interp(r2,r,P))*psi, r2) \
    #        + ( np.interp(r2, r, K_rho_Gamma1) + ( np.interp(r2, r, Gamma_1P) \
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
    
    #print(result)
    #psi = result.sol(r)[0] #interp1d(result['x'], result['y'][0])(r)
    #print('psi', psi)
    
    #result = least_squares(resid, np.ones(len(r)))
    #print(result)
    #psi = result['x']
    #psi = minimize(sqdiff, np.ones(len(r)), method='Nelder-Mead')['x']
    
#    K_u_Y = P * UnivariateSpline(P**-1 * psi, r).derivative()(r) \
#        + Gamma_1P * K_Gamma1_rho
    
    
    return { 
        ('c',      'rho'): (K_c_rho,      K_rho_c2),
        ('c2',     'rho'): (K_c2_rho,     K_rho_c2),
        ('Gamma1', 'rho'): (K_Gamma1_rho, K_rho_Gamma1),
        ('u',      'Y'):   (K_u_Y,        K_Y_u)#,
        #('psi', 'psi'):    (psi,          psi)
    }

