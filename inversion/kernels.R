#### Kernel pairs 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

u.name <- 'Sq. iso. sound speed' #'Pressure per density'#

c2_rho <- list(name='(c^2, rho)', 
    short='c2rho', 
    f1='c2', 
    f2='rho', 
    f1.name='Squared sound speed',
    f2.name='Density',
    f1.exp=bquote(c^2), 
    f2.exp=bquote(rho),
    f1.units=bquote(cm^2~s^-2),
    f2.units=bquote(g~cm^-3))

u_Y <- list(name='(u, Y)', 
    short='uY', 
    f1='u', 
    f2='Y',
    f1.name=u.name,
    f2.name='Helium abundance',
    f1.exp=bquote(u), 
    f2.exp=bquote(Y),
    f1.units=bquote(cm^2~s^-2),
    f2.units="")

u_Gamma1 <- list(name='(u, Gamma1)', 
    short='uGamma1', 
    f1='u', 
    f2='Gamma1',
    f1.name=u.name,
    f2.name='First adiabatic exponent',
    f1.exp=bquote(u), 
    f2.exp=bquote(Gamma[1]),
    f1.units=bquote(cm^2~s^-2),
    f2.units=bquote())

rho_Y <- list(name='(rho, Y)', 
    short='rhoY', 
    f1='rho', 
    f2='Y',
    f1.name='Density',
    f2.name='Helium abundance',
    f1.exp=bquote(rho), 
    f2.exp=bquote(Y),
    f1.units=bquote(g~cm^-3),
    f2.units=bquote())

Gamma1_rho <- list(name='(Gamma1, rho)', 
    short='Gamma1rho', 
    f1='Gamma1', 
    f2='rho',
    f1.name='First adiabatic exponent',
    f2.name='Density',
    f1.exp=bquote(Gamma[1]), 
    f2.exp=bquote(rho),
    f1.units=bquote(),
    f2.units=bquote(g~cm^-3))

c2_Gamma1 <- list(name='(c^2, Gamma1)', 
    short='c2Gamma1', 
    f1='c2', 
    f2='Gamma1',
    f1.name='Squared sound speed',
    f2.name='First adiabatic exponent',
    f1.exp=bquote(c^2), 
    f2.exp=bquote(Gamma[1]),
    f1.units=bquote(cm^2~s^-2),
    f2.units=bquote())






rho_c2 <- list(name='(rho, c^2)', 
    short='rhoc2', 
    f1='rho', 
    f2='c2', 
    f1.name='Density',
    f2.name='Squared sound speed',
    f1.exp=bquote(rho), 
    f2.exp=bquote(c^2),
    f1.units=bquote(g~cm^-3),
    f2.units=bquote(cm^2~s^-2))

Y_u <- list(name='(Y, u)', 
    short='Yu', 
    f1='Y', 
    f2='u',
    f1.name='Helium abundance',
    f2.name=u.name,
    f1.exp=bquote(Y), 
    f2.exp=bquote(u),
    f1.units=bquote(),
    f2.units=bquote(cm^2~s^-2))

Gamma1_u <- list(name='(Gamma1, u)', 
    short='Gamma1u', 
    f1='Gamma1', 
    f2='u',
    f1.name='First adiabatic exponent',
    f2.name=u.name,
    f1.exp=bquote(Gamma[1]), 
    f2.exp=bquote(u),
    f1.units=bquote(),
    f2.units=bquote(cm^2~s^-2))

Y_rho <- list(name='(Y, rho)', 
    short='Yrho', 
    f1='Y', 
    f2='rho',
    f1.name='Helium abundance',
    f2.name='Density',
    f1.exp=bquote(Y), 
    f2.exp=bquote(rho),
    f1.units=bquote(),
    f2.units=bquote(g~cm^-3))

rho_Gamma1 <- list(name='(rho, Gamma1)', 
    short='rhoGamma1', 
    f1='rho', 
    f2='Gamma1',
    f1.name='Density',
    f2.name='First adiabatic exponent',
    f1.exp=bquote(rho), 
    f2.exp=bquote(Gamma[1]),
    f1.units=bquote(g~cm^-3),
    f2.units=bquote())

Gamma1_c2 <- list(name='(Gamma1, c^2)', 
    short='Gamma1c2', 
    f1='Gamma1', 
    f2='c2',
    f1.name='First adiabatic exponent',
    f2.name='Squared sound speed',
    f1.exp=bquote(Gamma[1]), 
    f2.exp=bquote(c^2),
    f1.units=bquote(),
    f2.units=bquote(cm^2~s^-2))

k.pairs <- list(c2_rho, u_Y, u_Gamma1, rho_Y, Gamma1_rho, c2_Gamma1,
    rho_c2, Y_u, Gamma1_u, Y_rho, rho_Gamma1, Gamma1_c2)

