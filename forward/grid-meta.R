#### Meta-data for grid.dat 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

get_label <- function(symbol) as.expression(bquote(
        .(seis.names[[symbol]])
      ~ .(seis.labs[[symbol]])
      * .(seis.units[[symbol]])
))

seis.names <- list(
  M              = "Mass", 
  Y              = "Initial helium", 
  Z              = "Initial metallicity",
  alpha          = "Mixing length parameter", 
  age            = "Age", 
  radius         = "Radius", 
  H              = "Fractional hydrogen", 
  He             = "Fractional helium", 
  Hc             = "Fractional core-hydrogen", 
  log_g          = "Surface gravity", 
  L              = "Luminosity", 
  Teff           = "Temperature", 
  Fe.H           = "Metallicity", 
  Dnu_median     = "Large frequency separation", 
  Dnu_slope      = "", 
  Dnu0_median    = "Large frequency separation", 
  Dnu0_slope     = "", 
  dnu02_median   = "Small frequency separation", 
  dnu02_slope    = "", 
  r_sep02_median = "", 
  r_sep02_slope  = "", 
  r_avg01_median = "", 
  r_avg01_slope  = "", 
  dnu13_median   = "Small frequency separation", 
  dnu13_slope    = "", 
  r_sep13_median = "", 
  r_sep13_slope  = "",
  r_avg10_median = "", 
  r_avg10_slope  = ""
)

seis.labs <- list(
  M              = bquote(M), 
  Y              = bquote(Y[0]), 
  Z              = bquote(Z[0]),
  alpha          = bquote(alpha["MLT"]), 
  age            = bquote(tau), 
  radius         = bquote(R), 
  H              = bquote(X(H)), 
  He             = bquote(X(He)), 
  Hc             = bquote(H[c]), 
  log_g          = bquote(log~g), 
  L              = bquote(L), 
  Teff           = bquote(T["eff"]), 
  Fe.H           = bquote("Fe"/"H"), 
  Dnu_median     = bquote("<"*Delta*nu*">"), 
  Dnu_slope      = bquote("<"*d*Delta*nu/d*nu*">"), 
  Dnu0_median    = bquote("<"*Delta*nu[0]*">"), 
  Dnu0_slope     = bquote("<"*d*Delta*nu[0]/d*nu*">"), 
  dnu02_median   = bquote("<"*delta*nu[0*","*2]*">"), 
  dnu02_slope    = bquote("<"*d*delta*nu[0*","*2]/d*nu*">"), 
  r_sep02_median = bquote("<"*r[0*","*2]*">"), 
  r_sep02_slope  = bquote("<"*d*r[0*","*2]/d*nu*">"), 
  r_avg01_median = bquote("<"*r[0*","*1]*">"), 
  r_avg01_slope  = bquote("<"*d*r[0*","*1]/d*nu*">"), 
  dnu13_median   = bquote("<"*delta*nu[1*","*3]*">"), 
  dnu13_slope    = bquote("<"*d*delta*nu[1*","*3]/d*nu*">"), 
  r_sep13_median = bquote("<"*r[1*","*3]*">"), 
  r_sep13_slope  = bquote("<"*d*r[1*","*3]/d*nu*">"),
  r_avg10_median = bquote("<"*r[1*","*0]*">"), 
  r_avg10_slope  = bquote("<"*d*r[1*","*0]/d*nu*">")
)

seis.units <- list(
  M              = bquote("/"*M["\u0298"]), 
  Y              = bquote(), 
  Z              = bquote(),
  alpha          = bquote(), 
  age            = bquote("/Gyr"), 
  radius         = bquote("/"*R["\u0298"]), 
  H              = bquote(), 
  He             = bquote(), 
  Hc             = bquote(), 
  log_g          = bquote(), 
  L              = bquote("/"*L["\u0298"]), 
  Teff           = bquote("/"*K), 
  Fe.H           = bquote(), 
  Dnu_median     = bquote("/"*mu*Hz), 
  Dnu_slope      = bquote(), 
  Dnu0_median    = bquote("/"*mu*Hz), 
  Dnu0_slope     = bquote(), 
  dnu02_median   = bquote("/"*mu*Hz), 
  dnu02_slope    = bquote(), 
  r_sep02_median = bquote("/"*mu*Hz), 
  r_sep02_slope  = bquote(), 
  r_avg01_median = bquote("/"*mu*Hz), 
  r_avg01_slope  = bquote(), 
  dnu13_median   = bquote("/"*mu*Hz), 
  dnu13_slope    = bquote(), 
  r_sep13_median = bquote("/"*mu*Hz), 
  r_sep13_slope  = bquote(),
  r_avg10_median = bquote("/"*mu*Hz), 
  r_avg10_slope  = bquote()
)

seis.latex <- list(
  M              = "M", 
  Y              = "$Y_0$", 
  Z              = "$Z_0$", 
  alpha          = "$\\alpha_{\\text{\"MLT\"}}$", 
  age            = "$\\tau$", 
  radius         = "R", 
  H              = "X(H)", 
  He             = "X(He)", 
  Hc             = "$H_c$", 
  log_g          = "lg $g$", 
  L              = "L", 
  Teff           = "$T_{\text{\"eff\"}}$", 
  Fe.H           = "Fe/H", 
  Dnu_median     = "$\\langle\\Delta\\nu\\rangle$", 
  Dnu_slope      = "$\\langle\\frac{d\\Delta\\nu}{d\nu}\\rangle$", 
  Dnu0_median    = "$\\langle\\Delta\\nu_0\\rangle$", 
  Dnu0_slope     = "$\\langle\\frac{d\\Delta\\nu_0}{d\nu}\\rangle$", 
  dnu02_median   = "$\\langle\\delta\\nu_{02}\\rangle$", 
  dnu02_slope    = "$\\langle\\frac{d\\delta\\nu_{02}}{d\nu}\\rangle$",
  r_sep02_median = "$\\langle r_{02}\\rangle$", 
  r_sep02_slope  = "$\\langle\\frac{dr_{02}}{d\nu}\\rangle$", 
  r_avg01_median = "$\\langle r_{01}\\rangle$",
  r_avg01_slope  = "$\\langle\\frac{dr_{01}}{d\nu}\\rangle$", 
  dnu13_median   = "$\\langle\\delta\\nu_{13}\\rangle$", 
  dnu13_slope    = "$\\langle\\frac{d\\delta\\nu_{13}}{d\nu}\\rangle$", 
  r_sep13_median = "$\\langle r_{13}\\rangle$", 
  r_sep13_slope  = "$\\langle\\frac{dr_{13}}{d\nu}\\rangle$", 
  r_avg10_median = "$\\langle r_{10}\\rangle$", 
  r_avg10_slope  = "$\\langle\\frac{dr_{10}}{d\nu}\\rangle$"
)

Z_levels <- list(
  M      = seq(0.7, 1.3, 0.1),
  Y      = seq(0.22, 0.34, 0.01),
  Z      = log10(seq(10**1e-04, 10**0.04, length=10)),
  alpha  = seq(1.5, 2.5, 0.1),
  age    = 0:14, 
  radius = seq(0.6, 2.1, 0.1),
  H      = seq(0.54, 0.78, 0.02),
  He     = seq(0.22, 0.45, 0.02),
  Hc     = seq(0, 0.78, 0.05)
)

color_levels <- list(
  M      = seq(0.7, 1.3, 0.025),
  Y      = seq(0.22, 0.34, 0.0025),
  Z      = log10(seq(10**1e-04, 10**0.04, length=20)),
  alpha  = seq(1.5, 2.5, 0.025),
  age    = seq(0, 13.8, 0.5),
  radius = seq(0.6, 2.1, 0.025),
  H      = seq(0.54, 0.78, 0.01),
  He     = seq(0.22, 0.45, 0.005),
  Hc     = seq(0, 0.78, 0.015)
)
