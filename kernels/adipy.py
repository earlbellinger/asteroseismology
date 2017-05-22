import struct
import numpy as np
from scipy.integrate import cumtrapz

"""
cs[23] = mode inertia
cs[26] = frequency without Richardson correction
cs[36] = frequency with Richardson correction
"""

def read_agsm(agsm):
    """Reads an ADIPLS grand summary file."

    Parameters
    ----------
    agsm: str
        Name of the grand summary file, usually starting or ending with agsm.

        
    Returns
    -------
    css: list of arrays
        The cs arrays for each mode.
    """
    
    with open(agsm, 'rb') as f:
        bin_file = f.read()

    css = []


def read_amde(amde):
    """Reads an ADIPLS eigenfunction file, written with nmode=1.

    Parameters
    ----------
    amde: str
        Name of the eigenfunction file, usually starting or ending with amde

        
    Returns
    -------
    css: list of arrays
        The cs arrays for each mode.
    eigs: list of arrays
        The eigenfunction arrays for each mode.
    """

    with open(amde, 'rb') as f:
        bin_file = f.read()

    css = []
    eigs = []

    while len(bin_file)>4:
        bin_file = bin_file[4:] # cut the leading fortran write statement

        fmt = '<' + 50*'d'
        size = struct.calcsize(fmt)
        css.append(np.array(struct.unpack(fmt, bin_file[:size])))
        bin_file = bin_file[size:]
        
        fmt = '<i'
        size = struct.calcsize(fmt)
        nnw = int(struct.unpack(fmt, bin_file[:size])[0])
        bin_file = bin_file[size:]
        
        fmt = '<' + 7*nnw*'d'
        size = struct.calcsize(fmt)
        eigs.append(np.array(struct.unpack(fmt, bin_file[:size])).reshape((-1,7)))
        bin_file = bin_file[size:]
        
        bin_file = bin_file[4:]

    return css, eigs

def read_amdl(amdl):
    """Reads an ADIPLS model file.

    Parameters
    ----------
    amdl: str
        Name of the model file, usually starting or ending with amdl

        
    Returns
    -------
    """

    with open(amdl, 'rb') as f:
        bin_file = f.read()

    bin_file = bin_file[4:]

    fmt = '<i'
    size = struct.calcsize(fmt)
    nmod = int(struct.unpack(fmt, bin_file[:size])[0])
    bin_file = bin_file[size:]

    nn = int(struct.unpack(fmt, bin_file[:size])[0])
    bin_file = bin_file[size:]

    fmt = '<8d'
    size = struct.calcsize(fmt)
    D = np.array(struct.unpack(fmt, bin_file[:size]))
    bin_file = bin_file[size:]

    fmt = '<' + 6*nn*'d'
    size = struct.calcsize(fmt)
    A = np.array(struct.unpack(fmt, bin_file[:size])).reshape((-1,6))

    return nmod, nn, D, A


def read_rkr(rkr_file):
    """Reads an ADIPLS rotational kernel file.

    Parameters
    ----------
    rkr: str
        Name of the eigenfunction file, usually starting or ending with amde

        
    Returns
    -------
    css: list of arrays
        The cs arrays for each mode.
    rkrs: list of arrays
        The kernel arrays for each mode.
    """

    with open(rkr_file, 'rb') as f:
        bin_file = f.read()

    css = []
    rkrs = []

    while len(bin_file)>4:
        bin_file = bin_file[4:] # cut the leading fortran write statement

        fmt = '<' + 50*'d'
        size = struct.calcsize(fmt)
        css.append(np.array(struct.unpack(fmt, bin_file[:size])))
        bin_file = bin_file[size:]
        
        fmt = '<i'
        size = struct.calcsize(fmt)
        nnw = int(struct.unpack(fmt, bin_file[:size])[0])
        bin_file = bin_file[size:]
        
        fmt = '<' + 2*nnw*'d'
        size = struct.calcsize(fmt)
        rkrs.append(np.array(struct.unpack(fmt, bin_file[:size])).reshape((-1,2)))
        bin_file = bin_file[size:]
        
        bin_file = bin_file[4:]

    return css, rkrs


def kernel(variable, ell, nu, eig, fgong):
    """Returns a structural kernel.  I have tried to make this as
    notationally similar to Gough & Thompson (1991) as possible.

    Parameters
    ----------
    variable: str
        Either 'rho' or 'c' to select whether a density or sound speed
        kernel is returned.
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

    def integrate(y, x):
        return np.hstack((0., cumtrapz(y, x)))

    def complement(y, x):
        return np.trapz(y, x) - integrate(y, x)

    omega = 2.*np.pi*nu                              # convert to cyclic frequency
    G = 6.672e-8                                     # gravitational constant
    L2 = ell*(ell+1)
    L = np.sqrt(L2)
    M, R = fgong['glob'][:2]                         # mass and radius from FGONG
    sigma = np.sqrt(R**3/G/M)*omega                  # dimensionless frequency

    r = fgong['var'][::-1,0]                         # radial co-ordinate
    m = M*np.exp(fgong['var'][::-1,1])               # mass co-ordinate
    P = fgong['var'][::-1,3]                         # pressure
    rho = fgong['var'][::-1,4]                       # density
    Gamma1 = fgong['var'][::-1,9]                    # first adiabatic index
    cs2 = Gamma1*P/rho                               # square of the sound speed
    A1 = (m/M)/(r/R)**3
    A2 = (G*m*rho)/(Gamma1*P*r)
    Vg = A2[:]
    drho_dr = -(Vg+fgong['var'][::-1,14])*rho/r      # density gradient

    x = eig[:,0]                                     # dimensionless radius (i.e. r/R)
    i0 = np.argmin((x-1.)**2)                        # do I actually use this!?
    xi_r = eig[:,1]*R                                # radial component of eigenfunction

    dxi_r_dr = np.hstack((0., np.diff(xi_r)/np.diff(r)))

    # chi is the "dilatation"
    if ell == 0:
        xi_h = 0.*xi_r # radial modes have zero horizontal component
        eta = G*m/r**3/omega**2
        chi = Vg/x*(eig[:,1]-sigma**2/A1/x*eig[:,2])
    elif ell > 0:
        xi_h = eig[:,2]*R/L2
        eta = L2*G*m/r**3/omega**2
        chi = Vg/x*(eig[:,1]-sigma**2/A1/L2*eig[:,2]-eig[:,3])
    else:
        raise ValueError('ell must be non-negative')

    # most stellar models include the central point, at which the
    # numerical value of chi and drho_dr might be buggy.
    chi[0] = 0.
    drho_dr[0] = 0.

    S = np.trapz((xi_r**2+L2*xi_h**2)*rho*r**2, r)

    if variable == 'c':
        K = rho*cs2*chi**2*r**2 # c.f. equation (60)
    elif variable == 'rho':
        # first compute the huge bracketed terms in last two lines of equation (61)
        K = (ell+1.)/r**ell*(xi_r-ell*xi_h)*integrate((rho*chi+xi_r*drho_dr)*r**(ell+2.), r) \
            - ell*r**(ell+1.)*(xi_r+(ell+1.)*xi_h)*complement((rho*chi+xi_r*drho_dr)*r**(1.-ell), r)

        # then combine it with the rest
        K = -0.5*(xi_r**2+L2*xi_h**2)*rho*omega**2*r**2 \
            +0.5*rho*cs2*chi**2*r**2 \
            -G*m*(chi*rho+0.5*xi_r*drho_dr)*xi_r \
            -4*np.pi*G*rho*r**2*complement((chi*rho+0.5*xi_r*drho_dr)*xi_r, r) \
            +G*m*rho*xi_r*dxi_r_dr \
            +0.5*G*(m*drho_dr+4.*np.pi*r**2*rho**2)*xi_r**2 \
            -4*np.pi*G*rho/(2.*ell+1.)*K
    else:
        raise ValueError('variable can only be c or rho')

    return K/(S*omega**2)


def kernels(variable, fgong, css, eigs):
    # read_amde returns css, eigs.
    Ks = []

    for cs, eig in zip(css, eigs):
        ell = int(cs[17])
        enn = int(cs[18])
        nu = cs[36]*1e-3
        if nu <= 0.:
            nu = cs[26]*1e-3

        Ks.append(kernel(variable, ell, nu, eig, fgong))
        
    return Ks


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
        
    
