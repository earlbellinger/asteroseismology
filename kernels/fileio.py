import numpy as np

def read_fgong(filename):
    """
    Read in an FGONG file.
    
    This function can read in FONG versions 250 and 300.
    
    The output gives the center first, surface last.
    
    **Example usage**
    
    >>> starg,starl = read_fgong('example.fgong')
    
    *Global stellar properties*
    
    The global stellar properties are stored in a normal dictionary:
    
    >>> print(starg)
    {'initial_x': 0.0, 'lam': 0.0, 'initial_z': 0.02, 'xi': 0.0, 'ddrho_drr_c': -121.4597632, 'na15': 0.0, 'na14': 18596.197489999999, 'photosphere_r': 371909348600.0, 'ddP_drr_c': -172.6751611, 'phi': 0.0, 'beta': 0.0, 'photosphere_L': 1.178675752e+37, 'star_age': 35245914.009999998, 'alpha': 1.8, 'star_mass': 1.39244e+34}

    You can get a list of all the keys in this dictionary, and access individual
    values via:
    
    >>> print(list(starg.keys()))
    ['initial_x', 'lam', 'initial_z', 'xi', 'ddrho_drr_c', 'na15', 'na14', 'photosphere_r', 'ddP_drr_c', 'phi', 'beta', 'photosphere_L', 'star_age', 'alpha', 'star_mass']
    >>> print(starg['star_mass'])
    1.39244e+34
    
    *Local stellar properties*
    
    The local stellar properties are stored in a record array. The benefit of
    a record array is that you can access the columns by there name. You can list
    all the column names with:
    
    >>> print(starl.dtype.names)
    ('radius', 'ln(m/M)', 'temperature', 'pressure', 'density', 'X', 'luminosity', 'opacity', 'eps_nuc', 'gamma1', 'grada', 'delta', 'cp', 'free_e', 'brunt_A', 'rx', 'Z', 'R-r', 'eps_grav', 'Lg', 'xhe3', 'xc12', 'xc13', 'xn14', 'xo16', 'dG1_drho', 'dG1_dp', 'dG1_dY', 'xh2', 'xhe4', 'xli7', 'xbe7', 'xn15', 'xo17', 'xo18', 'xne20', 'xh1', 'na38', 'na39', 'na40')
    
    An entire column is accessed via e.g. ``starl['radius']``. To get only
    the first 5 entries, you can do
    
    >>> print(starl['radius'][:5])
    [  0.00000000e+00   2.63489963e+08   3.31977388e+08   4.18267440e+08
       5.26988298e+08]
    
    You can sub select multiple columns simultaneously, and have a nice string
    representation via
    
    >>> print(plt.mlab.rec2txt(starl[['radius','temperature']][:5]))
              radius    temperature
               0.000   32356218.190
       263489963.200   32355749.110
       331977388.500   32355443.030
       418267440.400   32354966.580
       526988298.100   32354211.870
       
    Record arrays support array indexing, so if you only want to select that
    part of the star with a temperature between 200000K and 5000000K, you can do
    (we only print the last 9 lines but plot everything):
    
    >>> keep = (200000<starl['temperature']) & (starl['temperature']<5000000)
    >>> print(plt.mlab.rec2txt(starl[keep][-9:]))
                 radius   ln(m/M)   temperature       pressure   density       X                                   luminosity   opacity   eps_nuc   gamma1   grada   delta               cp   free_e   brunt_A       rx       Z               R-r   eps_grav      Lg    xhe3    xc12    xc13    xn14    xo16   dG1_drho   dG1_dp   dG1_dY     xh2    xhe4    xli7    xbe7    xn15    xo17    xo18   xne20     xh1    na38    na39    na40
       355098749500.000    -0.000    217776.342   28893242.570     0.000   0.700   11786757520000000148929153677249216512.000     6.093    -0.057    1.512   0.295   1.979   1116313577.000    0.843     4.780   -0.000   0.020   16810599020.000     -0.057   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       355254331000.000    -0.000    215959.487   27989491.870     0.000   0.700   11786757520000000148929153677249216512.000     6.146    -0.055    1.512   0.295   1.977   1114375987.000    0.843     4.440   -0.000   0.020   16655017510.000     -0.055   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       355386669700.000    -0.000    214408.595   27237358.180     0.000   0.700   11786757520000000148929153677249216512.000     6.189    -0.053    1.512   0.295   1.975   1112603487.000    0.843     4.144   -0.000   0.020   16522678890.000     -0.053   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       355603756700.000    -0.000    211839.378   26036219.670     0.000   0.700   11786757520000000148929153677249216512.000     6.258    -0.049    1.512   0.295   1.971   1109132366.000    0.843     3.641   -0.000   0.020   16305591810.000     -0.049   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       355827789500.000    -0.000    209133.425   24838430.420     0.000   0.700   11786757520000000148929153677249216512.000     6.329    -0.045    1.513   0.295   1.966   1104338120.000    0.843     3.088   -0.000   0.020   16081559030.000     -0.045   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       355971535600.000    -0.000    207374.731   24091371.390     0.000   0.700   11786757520000000148929153677249216512.000     6.372    -0.043    1.513   0.295   1.962   1100769436.000    0.843     2.721   -0.000   0.020   15937812930.000     -0.043   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       356088698500.000    -0.000    205931.643   23494524.550     0.000   0.700   11786757520000000148929153677249216512.000     6.406    -0.041    1.513   0.295   1.959   1097649965.000    0.843     2.418   -0.000   0.020   15820650060.000     -0.041   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       356207860400.000    -0.000    204455.624   22898391.190     0.000   0.700   11786757520000000148929153677249216512.000     6.438    -0.038    1.513   0.296   1.955   1094296115.000    0.843     2.111   -0.000   0.020   15701488190.000     -0.038   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000
       356390516500.000    -0.000    202170.683   22005682.050     0.000   0.700   11786757520000000148929153677249216512.000     6.483    -0.035    1.514   0.296   1.949   1088675748.000    0.843     1.638   -0.000   0.020   15518832090.000     -0.035   0.000   0.000   0.003   0.000   0.001   0.010      0.000    0.000    0.000   0.000   0.280   0.000   0.000   0.000   0.000   0.000   0.002   0.000   0.000   0.000   0.000

    You can easily plot this:
    
    >>> p = plt.figure()
    >>> p = plt.plot(starl['radius'],starl['temperature'],'k-')
    >>> p = plt.plot(starl[keep]['radius'],starl[keep]['temperature'],'r-',lw=2)
    
    .. image:: fileio_read_fgong_profile.png
    
    @param filename: name of the FGONG file
    @type filename: str
    @return: global parameters, local parameters
    @rtype: dict, record array
    """
    #   these are the standard definitions and units for FGONG
    glob_pars = [('star_mass','g'),
                 ('photosphere_r','cm'),
                 ('photosphere_L','erg/s'),
                 ('initial_z',None),
                 ('initial_x',None),
                 ('alpha',None),
                 ('phi',None),
                 ('xi',None),
                 ('beta',None),
                 ('lam',None),
                 ('ddP_drr_c',None),
                 ('ddrho_drr_c',None),
                 ('star_age','yr'),
                 ('na14',None),
                 ('na15',None)]
    loc_pars = [('radius','cm'),
                ('ln(m/M)',None),
                ('temperature','K'),
                ('pressure','kg/m/s2'),
                ('density','g/cm3'),
                ('X',None),
                ('luminosity','erg/s'),
                ('opacity','cm2/g'),
                ('eps_nuc',None),
                ('gamma1',None),
                ('grada',None),
                ('delta',None),
                ('cp',None),
                ('free_e',None),
                ('brunt_A',None),
                ('rx',None),
                ('Z',None),
                ('R-r','cm'),
                ('eps_grav',None),
                ('Lg','erg/s'),
                ('xhe3',None),
                ('xc12',None),
                ('xc13',None),
                ('xn14',None),
                ('xo16',None),
                ('dG1_drho',None),
                ('dG1_dp',None),
                ('dG1_dY',None),
                ('xh2',None),
                ('xhe4',None),
                ('xli7',None),
                ('xbe7',None),
                ('xn15',None),
                ('xo17',None),
                ('xo18',None),
                ('xne20',None),
                ('xh1',None),
                ('na38',None),
                ('na39',None),
                ('na40',None)]
    #-- start reading the file
    ff = open(filename,'r')
    lines = ff.readlines()
    #-- skip lines until we are at the definitions line
    while not len(lines[0].strip().split())==4:
        lines = lines[1:]
    #-- take care of data dimensions
    NN,ICONST,IVAR,IVERS = [int(i) for i in lines[0].strip().split()]
    if not ICONST==15:
        raise(ValueError, 'cannot interpret FGONG file: wrong number of global parameters')
    if not IVERS in [300,250]:
        raise(ValueError, 'cannot interpret FGONG file: wrong format version (version = {})'.format(IVERS))
    data = []
    #-- read in all the data
    for line in lines[1:]:
        data.append([line[0*16:1*16],line[1*16:2*16],line[2*16:3*16],
                     line[3*16:4*16],line[4*16:5*16]])
    data = np.ravel(np.array(data,float))
    starg = {glob_pars[i][0]:data[i] for i in range(ICONST)}
    data = data[15:].reshape((NN,IVAR)).T
    #-- reverse the profile to get center ---> surface
    data = data[:,::-1]
    #-- make it into a record array and return the data
    starl = np.rec.fromarrays(data,names=[lp[0] for lp in loc_pars])
    return starg,starl

