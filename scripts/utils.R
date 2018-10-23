#### Utility functions for forward & inverse asteroseismology
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Stellar Ages & Galactic Evolution Group
#### Max-Planck-Institut fur Sonnensystemforschung

invisible(library(magicaxis))
invisible(library(extrafont))

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs9.21/bin/gswin64.exe")

## Plotting values
utils.mar <<- c(3, 4, 1, 1)
utils.mgp <<- c(2, 0.25, 0)
#paper.mgp <<- c(2, 0.15, 0)
utils.font <<- "Palatino Linotype"#"Times" #"Palatino Linotype" #"CM Roman" #"Helvetica" #"Times" #
font <- "Palatino Linotype"#"Times" #"Palatino Linotype" #"CM Roman" #"Helvetica" #"Palatino"
text.cex <- 1
mgp <- c(2, 0.25, 0)
utils.tcl <<- -0.346
#hack.mgp <- c(2, 0.5, 0)


## Constants
solar_Teff = log10(5777)
solar_age = 4.57e9
#solar_radius = 6.955*10**10 # cm 
solar_radius = 6.9598e10 # cm
solar_mass = 1.9892e33 # g
cgrav = 6.67428*10**-8 # cm^3 * g^-1 * s^-2
#solar_mass = 1.9891*10**30 # kg 
solar_scale = sqrt(solar_mass/solar_radius^3)
FeH <- function(Z, H1, Z_div_X_solar=0.02293) log10(Z / H1 / Z_div_X_solar)

# e/Fe_solar for every element according to GS98 
# abundances calculated with log10(Fe/H)+12 
# where Fe=7.5, H=12. 
# solve for Fe, then calculate log10(e/H)+12 and e/Fe 
H_Fe_solar  = 31622.7766017
He_Fe_solar = 2691.53480393
Li_Fe_solar = 3.98107170554e-07
Be_Fe_solar = 7.94328234725e-07
B_Fe_solar  = 1.1220184543e-05
C_Fe_solar  = 10.4712854805
N_Fe_solar  = 2.6302679919
O_Fe_solar  = 21.379620895
F_Fe_solar  = 0.0011481536215
Ne_Fe_solar = 3.80189396321
Na_Fe_solar = 0.0676082975392
Mg_Fe_solar = 1.20226443462
Al_Fe_solar = 0.0933254300797
Si_Fe_solar = 1.1220184543
P_Fe_solar  = 0.00891250938134
S_Fe_solar  = 0.676082975392
Cl_Fe_solar = 0.01
Ar_Fe_solar = 0.0794328234725
K_Fe_solar  = 0.00416869383471
Ca_Fe_solar = 0.0724435960075
Sc_Fe_solar = 4.67735141287e-05
Ti_Fe_solar = 0.00331131121483
V_Fe_solar  = 0.000316227766017
Cr_Fe_solar = 0.0147910838817
Mn_Fe_solar = 0.00776247116629
Fe_Fe_solar = 1.0
Co_Fe_solar = 0.0026302679919
Ni_Fe_solar = 0.0562341325191
Cu_Fe_solar = 0.000512861383992
Zn_Fe_solar = 0.00125892541179
Ga_Fe_solar = 2.39883291902e-05
Ge_Fe_solar = 8.12830516165e-05
Rb_Fe_solar = 1.25892541179e-05
Sr_Fe_solar = 2.95120922667e-05
Y_Fe_solar  = 5.49540873858e-06
Zr_Fe_solar = 1.25892541179e-05
Nb_Fe_solar = 8.31763771103e-07
Mo_Fe_solar = 2.6302679919e-06
Ru_Fe_solar = 2.18776162395e-06
Rh_Fe_solar = 4.16869383471e-07
Pd_Fe_solar = 1.54881661891e-06
Ag_Fe_solar = 2.75422870334e-07
Cd_Fe_solar = 1.86208713666e-06
In_Fe_solar = 1.44543977075e-06
Sn_Fe_solar = 3.16227766017e-06
Sb_Fe_solar = 3.16227766017e-07
Ba_Fe_solar = 4.26579518802e-06
La_Fe_solar = 4.67735141287e-07
Ce_Fe_solar = 1.20226443462e-06
Pr_Fe_solar = 1.62181009736e-07
Nd_Fe_solar = 1e-06
Sm_Fe_solar = 3.2359365693e-07
Eu_Fe_solar = 1.02329299228e-07
Gd_Fe_solar = 4.16869383471e-07
Dy_Fe_solar = 4.3651583224e-07
Ho_Fe_solar = 5.75439937337e-08
Er_Fe_solar = 2.69153480393e-07
Yb_Fe_solar = 3.80189396321e-07
Lu_Fe_solar = 3.6307805477e-08
Hf_Fe_solar = 2.39883291902e-07
W_Fe_solar  = 4.07380277804e-07
Os_Fe_solar = 8.91250938134e-07
Ir_Fe_solar = 7.07945784385e-07
Pt_Fe_solar = 1.99526231497e-06
Au_Fe_solar = 3.2359365693e-07
Tl_Fe_solar = 2.51188643151e-07
Pb_Fe_solar = 2.81838293127e-06

#png_res <- 400
#cex.paper <- 0.8
#cex.slides <- 1.3
#cex.hack <- 1.4
#latex_pt_per_in <- 5 * 72.27

blue <- "#0571b0"
red <- "#CA0020"
orange <- "#F97100"

################################################################################
### PLOTTING ###################################################################
################################################################################
### Make the same plot as a pdf and a png suitable for papers and for slides
# takes a plotting function plot_f that calls `plot`
#
# template for calling:
## my_plotting_function <- function(
##     my_plotting_function_arg1,
##     my_plotting_function_arg2,
##     ...,
##     text.cex=1, mgp=utils.mgp, font=utils.font)
## make_plots(my_plotting_function, 'filename_without_extension',
##     my_plotting_function_arg1,
##     my_plotting_function_arg2,
##     ...)
# note 1: the ...'s are real, include them in the function declaration & call
# note 2: include any hidden arguments that you want in my_plotting_function;
#         they are all defined directly below; they are the things like
#         'filepath', which defaults to 'plots' but can be anything you like, or
#         'wide' and 'make_png', which you can turn off if you don't want those
# note 3: personally I always pass along 'text.cex', 'mgp', and 'font'
#         so that I can do
#         ## magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0),
#         ##         las=1, mgp=mgp, family=font, cex.axis=text.cex)
make_plots <- function(plot_f, filename,
        filepath='plots',
        mar=utils.mar,
        mar.paper=c(2.5, 3, 1, 1),
        mgp.slides=utils.mgp,
        mgp.paper=c(1.5, 0.15, 0),
        oma=c(0,0,0,0),
        hack.mgp=c(2, 0.5, 0), thin.hack=FALSE,
        cex.paper=0.77, cex.slides=1.3, cex.hack=1.4,
        wide=TRUE, thin=TRUE,
        tall=TRUE, short=TRUE,
        make_png=TRUE, make_pdf=TRUE, make_tikz=FALSE,
        paper=TRUE, slides=TRUE,
        paper_pdf_width=6.97522, # inches
        paper_pdf_height=4.17309,
        slides_pdf_width=6.22665,
        slides_pdf_height=4.1511,
        latex_pt_per_in=5*72.27,
        png_res=400,
        font=utils.font, 
        use.cairo=T, ...) {

    paper_png_width <- paper_pdf_width * latex_pt_per_in
    paper_png_height <- paper_pdf_height * latex_pt_per_in
    slides_png_width <- slides_pdf_width * latex_pt_per_in
    slides_png_height <- slides_pdf_height * latex_pt_per_in

    args. <- c(as.list(environment()), list(...))

    if (slides) {
        do.call(widethin, c(list(
                directory=file.path(filepath, 'slides'),
                slides_pdf_width,
                slides_pdf_height,
                slides_png_width,
                slides_png_height,
                text.cex=cex.slides,
                mgp=mgp.slides),
            args.))
    }
    if (paper) {
        if (all(args.$mar == utils.mar)) args.$mar <- mar.paper
        do.call(widethin, c(list(
                directory=file.path(filepath, 'paper'),
                pdf_width=paper_pdf_width,
                pdf_height=paper_pdf_height,
                png_width=paper_png_width,
                png_height=paper_png_height,
                mgp=mgp.paper,
                #mar=mar.paper,
                text.cex=cex.paper),
            args.))
    }
}

widethin <- function(plot_f, filename, directory,
        pdf_width, pdf_height, png_width, png_height, text.cex, use.cairo, ...,
        wide=T, thin=T) {

    if (thin) {
        thin_dir <- file.path(directory, 'thin')
        dir.create(thin_dir, showWarnings=FALSE, recursive=TRUE)
        tallshort(plot_f, filename, thin_dir,
                  pdf_width/2, png_width/2,
                  pdf_height, png_height, text.cex, use.cairo, ...,
                  thin=T, wide=F)
    }

    if (wide) {
        wide_dir <- file.path(directory, 'wide')
        dir.create(wide_dir, showWarnings=FALSE, recursive=TRUE)
        tallshort(plot_f, filename, wide_dir,
                  pdf_width, png_width,
                  pdf_height, png_height, text.cex, use.cairo, ...,
                  thin=F, wide=T)
    }
}

tallshort <- function(plot_f, filename, directory,
        pdf_width, png_width, pdf_height, png_height, text.cex, use.cairo, ...,
        tall=T, short=T, thin.hack=F) {

    if (short) {
        short_dir <- file.path(directory, 'short')
        dir.create(short_dir, showWarnings=FALSE, recursive=TRUE)
        pdfpng(plot_f, filename, short_dir,
               pdf_width, pdf_height/2,
               png_width, png_height/2, text.cex, use.cairo, ...,
               tall=F, short=T)
    }
    if (tall) {
        tall_dir <- file.path(directory, 'tall')
        dir.create(tall_dir, showWarnings=FALSE, recursive=TRUE)
        pdfpng(plot_f, filename, tall_dir,
               pdf_width, pdf_height,
               png_width, png_height, text.cex, use.cairo, ...,
               tall=T, short=F)
    }
}

pdfpng <- function(plot_f, filename, directory,
        pdf_width, pdf_height, png_width, png_height,
        text.cex, png_res, use.cairo, ..., 
        font=utils.font, mar=utils.mar, thin.hack=F,
        make_png=T, make_pdf=T, make_tikz=T, mgp=utils.mgp, oma=c(0,0,0,0)) {
    if (make_png) {
        png_font <- if (font == "CM Roman") "Times" else font 
        png(file.path(directory, paste0(filename, '.png')),
            width=png_width, height=png_height,
            family=png_font, res=png_res, type='cairo')
        par(mar=mar, mgp=mgp, cex.lab=text.cex, family=png_font, oma=oma)
        plot_f(text.cex=if (thin.hack) cex.hack else text.cex,
               font=png_font,
               mgp=if (thin.hack) hack.mgp else mgp, ...)
        dev.off()
    }
    if (make_pdf) {
        cairornot <- if (use.cairo) cairo_pdf else pdf
        cairornot(file.path(directory, paste0(filename, '.pdf')),
            width=pdf_width, height=pdf_height, family=font)
        par(mar=mar, mgp=mgp, cex.lab=text.cex, family=font, oma=oma)
        plot_f(text.cex=if (thin.hack) cex.hack else text.cex,
               font=font, mgp=if (thin.hack) hack.mgp else mgp, ...)
        dev.off()
        if (font == "CM Roman") {
            embed_fonts(file.path(directory, paste0(filename, '.pdf')), 
                outfile=file.path(directory, paste0(filename, '.pdf')))
        }
    }
    if (make_tikz) {
        tikz(file.path(directory, paste0(filename, '.tex')),
            width=pdf_width, height=pdf_height)
        par(mar=mar, mgp=mgp, cex.lab=text.cex, family=font, oma=oma)
        plot_f(text.cex=if (thin.hack) cex.hack else text.cex,
               font=font, mgp=if (thin.hack) hack.mgp else mgp, ...)
        dev.off()
    }
}

################################################################################
### FGONG ######################################################################
################################################################################
read.fgong <- function(filename) {
    ## Create connection
    #con <- file(description=file, open="r")
    
    ## Hopefully you know the number of lines from some other source or
    #com <- paste("wc -l ", file, " | awk '{ print $1 }'", sep="")
    #n <- system(command=com, intern=TRUE)
    
    ## Loop over a file connection
    #for(i in 1:n) {
    #    tmp <- scan(file=con, nlines=1, quiet=TRUE)
        ## do something on a line of data 
    #}
}

################################################################################
### Nearly-equal spacing #######################################################
################################################################################
find_closest <- function(x, y) {
    x.dest <- x
    y.dest <- y
    x.ind <- c()
    y.ind <- c()
    min.cost <- Inf
    ii <- 0
    while (length(x.dest) > 0 && length(y.dest) > 0) {
        ii <- ii + 1
        cost.mat <- outer(x.dest, y.dest, function(x, y) abs(x-y))
        cost <- min(cost.mat)
        indices <- which(cost.mat==cost, arr.ind=T)
        if (cost >= 10*min.cost && ii > 5) break
        if ((cost < min.cost || ii < 5) && cost > 0) 
            min.cost <- min(cost.mat)
        x.ind <- c(x.ind, which(x==x.dest[indices[1]])[1])
        y.ind <- c(y.ind, which(y==y.dest[indices[2]])[1])
        x.dest <- x.dest[-indices[1]]
        y.dest <- y.dest[-indices[2]]
    }
    list(x=x.ind, y=y.ind)
}

# solve linear transport problem to get equally-spaced points
find_closest2 <- function(x, y=NULL, num_points=0) {
    invisible(library(lpSolve))

    x.len <- length(x)

    if (x.len < num_points) {
        print(paste("Too few points"))
        return(NULL)
    }

    if (is.null(y) && num_points > 0)
        y <- seq(max(x), min(x), length=num_points)
    if (!is.null(y) && num_points == 0)
        num_points <- length(y)

    cost.mat  <- outer(y, x, function(x, y) abs(x-y))
    row.signs <- rep("==", num_points)
    row.rhs   <- rep(1, num_points)
    col.signs <- rep("<=", x.len)
    col.rhs   <- rep(1, x.len)
    sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
        col.signs, col.rhs)$solution
    apply(sol, 1, which.max)
}

################################################################################
### Uncertainty calculations ###################################################
################################################################################

gumr <- function(x.n,x.u) {
  z2 <- trunc(log10(x.u))+1
  z1 <- round(x.u/(10^z2), 3)
  y1 <- round(x.n*10^(-z2), 2)
  list(value=y1*10^z2,uncert=z1*10^z2)
}

################################################################################
### Metadata for the grid ######################################################
################################################################################

get_label <- function(symbol) bquote(
        .(seis.names[[symbol]])
      ~ .(seis.labs[[symbol]])
      * .(seis.units[[symbol]])
)

get_label_nameless <- function(symbol) bquote(
        .(seis.labs[[symbol]])
      * .(seis.units[[symbol]])
)

seis.names <- list(
  M               = "Mass",
  Y               = "Initial helium",
  Z               = "Initial metallicity",
  alpha           = "Mixing length parameter",
  diffusion       = "Diffusion multiplication factor",
  overshoot       = "Overshoot",
  age             = "Age",
  radius          = "Radius",
  mass_X          = "Hydrogen mass fraction",
  mass_Y          = "Helium mass fraction",
  mass_cc         = "Convective core mass",
  X_surf          = "Surface hydrogen",
  Y_surf          = "Surface helium",
  X_c             = "Core hydrogen",
  log_g           = "Surface gravity",
  L               = "Luminosity",
  Teff            = "Temperature",
  Fe.H            = "Metallicity",
  Dnu_median      = "Large frequency separation",
  Dnu             = "Large frequency separation",
  Dnu_slope       = "",
  Dnu0_median     = "Large frequency separation",
  Dnu0            = "Large frequency separation",
  Dnu0_slope      = "",
  dnu02_median    = "Small separation",
  dnu02           = "Small separation",
  dnu02_slope     = "",
  r_sep02_median  = "",
  r02             = "",
  r_sep02_slope   = "",
  r_avg01_median  = "",
  r01             = "",
  r_avg01_slope   = "",
  dnu13_median    = "Small separation",
  dnu13           = "Small separation",
  dnu13_slope     = "",
  r_sep13_median  = "",
  r13             = "",
  r_sep13_slope   = "",
  r_avg10_median  = "",
  r10             = "",
  r_avg10_slope   = "",
  radial_velocity = "Radial velocity",
  nu_max          = ""
)

seis.labs <- list(
  M               = bquote(M),
  Y               = bquote(Y[0]),
  Z               = bquote(Z[0]),
  alpha           = bquote(alpha["MLT"]),
  diffusion       = bquote(D),
  overshoot       = bquote(alpha["ov"]),
  age             = bquote(tau),
  radius          = bquote(R),
  mass_X          = bquote(X),
  mass_Y          = bquote(Y),
  mass_cc         = bquote(M["cc"]),
  X_surf          = bquote(X["surf"]),
  Y_surf          = bquote(Y["surf"]),
  X_c             = bquote(X[c]),
  log_g           = bquote(log~g),
  L               = bquote(L),
  Teff            = bquote(T["eff"]),
  Fe.H            = bquote("[Fe"/"H]"),
  Dnu_median      = bquote("<"*Delta*nu*">"),
  Dnu             = bquote("<"*Delta*nu*">"),
  Dnu_slope       = bquote("<"*d*Delta*nu/d*nu*">"),
  Dnu0_median     = bquote("<"*Delta*nu[0]*">"),
  Dnu0            = bquote("<"*Delta*nu[0]*">"),
  Dnu0_slope      = bquote("<"*d*Delta*nu[0]/d*nu*">"),
  dnu02_median    = bquote("<"*delta*nu[0*","*2]*">"),
  dnu02           = bquote("<"*delta*nu[0*","*2]*">"),
  dnu02_slope     = bquote("<"*d*delta*nu[0*","*2]/d*nu*">"),
  r_sep02_median  = bquote("<"*r[0*","*2]*">"),
  r02             = bquote("<"*r[0*","*2]*">"),
  r_sep02_slope   = bquote("<"*d*r[0*","*2]/d*nu*">"),
  r_avg01_median  = bquote("<"*r[0*","*1]*">"),
  r01             = bquote("<"*r[0*","*1]*">"),
  r_avg01_slope   = bquote("<"*d*r[0*","*1]/d*nu*">"),
  dnu13_median    = bquote("<"*delta*nu[1*","*3]*">"),
  dnu13           = bquote("<"*delta*nu[1*","*3]*">"),
  dnu13_slope     = bquote("<"*d*delta*nu[1*","*3]/d*nu*">"),
  r_sep13_median  = bquote("<"*r[1*","*3]*">"),
  r13             = bquote("<"*r[1*","*3]*">"),
  r_sep13_slope   = bquote("<"*d*r[1*","*3]/d*nu*">"),
  r_avg10_median  = bquote("<"*r[1*","*0]*">"),
  r10             = bquote("<"*r[1*","*0]*">"),
  r_avg10_slope   = bquote("<"*d*r[1*","*0]/d*nu*">"),
  radial_velocity = bquote(r[v]),
  nu_max          = bquote(nu["max"]),
  PC1             = bquote(PC[1]),
  PC2             = bquote(PC[2])
)

seis.units <- list(
  Y               = bquote(),
  Z               = bquote(),
  alpha           = bquote(),
  diffusion       = bquote(),
  overshoot       = bquote(),
  age             = bquote("/Gyr"),
  mass_X          = bquote("/"*M["*"]),
  mass_Y          = bquote("/"*M["*"]),
  mass_Y          = bquote("/"*M["*"]),
  X_surf          = bquote(),
  Y_surf          = bquote(),
  X_c             = bquote(),
  log_g           = bquote(" (cgs)"),
  Teff            = bquote("/"*K),
  Fe.H            = bquote(),
  Dnu_median      = bquote("/"*mu*Hz),
  Dnu             = bquote("/"*mu*Hz),
  Dnu_slope       = bquote(),
  Dnu0_median     = bquote("/"*mu*Hz),
  Dnu0            = bquote("/"*mu*Hz),
  Dnu0_slope      = bquote(),
  dnu02_median    = bquote("/"*mu*Hz),
  dnu02           = bquote("/"*mu*Hz),
  dnu02_slope     = bquote(),
  r_sep02_median  = bquote(""),
  r02             = bquote(""),
  r_sep02_slope   = bquote(),
  r_avg01_median  = bquote(""),
  r01             = bquote(""),
  r_avg01_slope   = bquote(),
  dnu13_median    = bquote("/"*mu*Hz),
  dnu13           = bquote("/"*mu*Hz),
  dnu13_slope     = bquote(),
  r_sep13_median  = bquote(""),
  r13             = bquote(""),
  r_sep13_slope   = bquote(),
  r_avg10_median  = bquote(""),
  r10             = bquote(""),
  r_avg10_slope   = bquote(),
  radial_velocity = bquote(m/s),
  nu_max          = bquote(mu*Hz)
)

seis.latex <- list(
  M               = "M",
  Y               = "$Y_0$",
  Z               = "$Z_0$",
  alpha           = "$\\alpha_{\\text{\"MLT\"}}$",
  age             = "$\\tau$",
  radius          = "R",
  mass_X          = "X",
  mass_Y          = "Y",
  X_surf          = "X_{\text{\"surf\"}}",
  Y_surf          = "Y_{\text{\"surf\"}}",
  X_c             = "$X_c$",
  log_g           = "lg $g$",
  L               = "L",
  Teff            = "$T_{\text{\"eff\"}}$",
  Fe.H            = "Fe/H",
  Dnu_median      = "$\\langle\\Delta\\nu\\rangle$",
  Dnu             = "$\\langle\\Delta\\nu\\rangle$",
  Dnu_slope       = "$\\langle\\frac{d\\Delta\\nu}{d\nu}\\rangle$",
  Dnu0_median     = "$\\langle\\Delta\\nu_0\\rangle$",
  Dnu0            = "$\\langle\\Delta\\nu_0\\rangle$",
  Dnu0_slope      = "$\\langle\\frac{d\\Delta\\nu_0}{d\nu}\\rangle$",
  dnu02_median    = "$\\langle\\delta\\nu_{02}\\rangle$",
  dnu02           = "$\\langle\\delta\\nu_{02}\\rangle$",
  dnu02_slope     = "$\\langle\\frac{d\\delta\\nu_{02}}{d\nu}\\rangle$",
  r_sep02_median  = "$\\langle r_{02}\\rangle$",
  r02             = "$\\langle r_{02}\\rangle$",
  r_sep02_slope   = "$\\langle\\frac{dr_{02}}{d\nu}\\rangle$",
  r_avg01_median  = "$\\langle r_{01}\\rangle$",
  r01             = "$\\langle r_{01}\\rangle$",
  r_avg01_slope   = "$\\langle\\frac{dr_{01}}{d\nu}\\rangle$",
  dnu13_median    = "$\\langle\\delta\\nu_{13}\\rangle$",
  dnu13           = "$\\langle\\delta\\nu_{13}\\rangle$",
  dnu13_slope     = "$\\langle\\frac{d\\delta\\nu_{13}}{d\nu}\\rangle$",
  r_sep13_median  = "$\\langle r_{13}\\rangle$",
  r13             = "$\\langle r_{13}\\rangle$",
  r_sep13_slope   = "$\\langle\\frac{dr_{13}}{d\nu}\\rangle$",
  r_avg10_median  = "$\\langle r_{10}\\rangle$",
  r10             = "$\\langle r_{10}\\rangle$",
  r_avg10_slope   = "$\\langle\\frac{dr_{10}}{d\nu}\\rangle$",
  radial_velocity = "$r_v$",
  nu_max          = "$\\nu_{\\max}$"
)


Z_levels <- list(
  M      = seq(0.7, 1.3, 0.1),
  Y      = seq(0.22, 0.34, 0.01),
  Z      = log10(seq(10**1e-04, 10**0.04, length=10)),
  alpha  = seq(1.5, 2.5, 0.1),
  age    = 0:14,
  radius = seq(0.6, 2.1, 0.1),
  mass_X = seq(0.54, 0.78, 0.02),
  mass_Y = seq(0.22, 0.45, 0.02),
  X_c    = seq(0, 0.78, 0.05)
)

color_levels <- list(
  M      = seq(0.7, 1.3, 0.025),
  Y      = seq(0.22, 0.34, 0.0025),
  Z      = log10(seq(10**1e-04, 10**0.04, length=20)),
  alpha  = seq(1.5, 2.5, 0.025),
  age    = seq(0, 13.8, 0.5),
  radius = seq(0.6, 2.1, 0.025),
  mass_X = seq(0.54, 0.78, 0.01),
  mass_Y = seq(0.22, 0.45, 0.005),
  X_c    = seq(0, 0.78, 0.015)
)

## Scatter plot function
scatter_plot <- function(seis.DF, X, Y, Z, combos, col.pal,
        xlim, ylim, xlab, ylab, log='', quartiles=F, ..., thin=F, short=F,
        solar_x=NA, solar_y=NA, text.cex=cex.paper, mgp=utils.mgp) {
    if (thin) par(mar=c(3, 4, 1, 4))
    Z_max <- if (quartiles) quantile(seis.DF[[Z]], 0.95) else max(seis.DF[[Z]])
    Z_min <- if (quartiles) quantile(seis.DF[[Z]], 0.05) else min(seis.DF[[Z]])
    for (simulation_i in 1:nrow(combos)) {
        track <- merge(seis.DF, combos[simulation_i,])
        idx <- floor((track[[Z]]-Z_min)/(Z_max-Z_min)*(length(col.pal)-1))+1
        idx[idx>length(col.pal)] <- length(col.pal)
        idx[idx<1] <- 1
        color <- col.pal[idx]
        use_line <- length(unique(seis.DF[[Z]]))==1
        ys <- track[[Y]]
        xs <- track[[X]]
        if (grepl('x', log) && 0 %in% xs) xs <- xs + min(xs[xs>0])
        if (grepl('y', log) && 0 %in% ys) ys <- ys + min(ys[ys>0])
        if (simulation_i == 1) {
            plot(NA, axes=FALSE, tcl=0,
                 xlab=if(thin) get_label_nameless(X) else xlab,
                 ylab=if(short) get_label_nameless(Y) else ylab,
                 log=log, xlim=xlim, ylim=ylim)
            if (!is.na(solar_x) && !is.na(solar_y)) {
                abline(v=solar_x, lty=3, col='black')
                abline(h=solar_y, lty=3, col='black')
            }
            magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=mgp, cex.axis=text.cex, las=1)
        }
        relation <- ys ~ xs
        if (use_line) lines(relation, col=adjustcolor(color, alpha=.5))
        else {
            segments(xs[-length(xs)], ys[-length(xs)], xs[-1], ys[-1],
                col=adjustcolor(color[-length(color)], alpha=.5))
            #points(relation, col=color, pch=20, cex=0.1)
        }
    }

    points(solar_x, solar_y, pch=1, cex=1)
    points(solar_x, solar_y, pch=20, cex=0.1)

    X_range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3],
                 par()$usr[2]+0.10*X_range, par()$usr[4],
                 signif(quantile(seq(Z_min, Z_max, length=1000),
                                 c(0, 0.25, 0.5, 0.75, 1)), 2),
                 col.pal[1:length(col.pal)],
                 cex=text.cex, gradient='y', align='rb')
    mtext(if(short) get_label_nameless(Z) else get_label(Z),
        4, line=ifelse(thin, 3, 5), cex=text.cex)
}

## A basic normalization function
normalize <- function(x) (x-min(x))/(max(x)-min(x))
standardize <- function(x) (x-mean(x))/sd(x)

## Get mesh
get_mesh <- function(X, Y, Z, seis.DF, cygA_stds) {
    library(akima)

    xx <- normalize(seis.DF[[X]])
    yy <- normalize(seis.DF[[Y]])

    #cygx <- cygA_stds[[X]]/(max(seis.DF[[X]])-min(seis.DF[[X]]))
    #cygy <- cygA_stds[[Y]]/(max(seis.DF[[Y]])-min(seis.DF[[Y]]))

    mesh <- interp(xx, yy, seis.DF[[Z]],
                xo=seq(0, 1, 0.02), #, cygx),
                yo=seq(0, 1, 0.02)) # cygy))

    mesh$x <- seq(min(seis.DF[[X]]), max(seis.DF[[X]]),
        (max(seis.DF[[X]])-min(seis.DF[[X]]))/50)
    #cygA_stds[[X]])
    mesh$y <- seq(min(seis.DF[[Y]]), max(seis.DF[[Y]]),
        (max(seis.DF[[Y]])-min(seis.DF[[Y]]))/50)
    #cygA_stds[[Y]])

    mesh
}

## Mesh plot function
mesh_plot <- function(mesh, X, Z, xlab, ylab, log='', ...,
        solar_x=NA, solar_y=NA, text.cex=cex.paper, mgp=utils.mgp) {
    filled.contour(mesh, log=log,
        xlim=if (X=='Teff') rev(range(mesh$x)) else range(mesh$x),
        levels=color_levels[[Z]],
        color=colorRampPalette(brewer.pal(11, "Spectral")),
        xaxs='i', yaxs='i',
        key.axes={
            axis(4, cex.axis=text.cex, tcl=0, line=0, mgp=mgp)
            mtext(get_label(Z), side=4, las=3, line=3, cex=text.cex)
        },
        plot.axes={
            contour(mesh, add=TRUE, labcex=0.5, levels=Z_levels[[Z]])
            if (!is.na(solar_x) && !is.na(solar_y)) {
                points(solar_x, solar_y, pch=1, cex=1)
                points(solar_x, solar_y, pch=20, cex=0.1)
            }
            magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=mgp, cex.axis=text.cex)
        },
        plot.title={
            title(xlab=xlab, cex.lab=text.cex, line=2)
            title(ylab=ylab, cex.lab=text.cex, line=2)
        })
}

## Calls scatter_plot and mesh_plot
# uses globals seis.DF, col.pal, and combos
scatter_mesh <- function(X, Y, Z, filepath=file.path("plots", "mesh"),
        scatter=1, mesh=1, log='', ...) {
    if (Y == 'Teff' && X == 'L') {
        Y = 'L'
        X = 'Teff'
    }

    xlab <- get_label(X)
    ylab <- get_label(Y)

    ylim <- quantile(seis.DF[[Y]], c(0.0001, 0.9999))
    xlim <- quantile(seis.DF[[X]], c(0.0001, 0.9999))
    if (X == 'Teff') xlim <- rev(xlim)

    has_solar <- any(grepl(X, names(solar_vals))) &&
                 any(grepl(Y, names(solar_vals)))
    solar_x <- if (has_solar) solar_vals[[X]] else NA
    solar_y <- if (has_solar) solar_vals[[Y]] else NA

    if (scatter) {
        make_plots(scatter_plot, paste0(paste(Z, Y, X, sep="_"), "-scatter"),
                   filepath=filepath, thin.hack=1,
                   mar=c(3, 4, 1, 6),
                   seis.DF=seis.DF, X=X, Y=Y, Z=Z, combos=combos,
                   col.pal=col.pal, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                   solar_x=solar_x, solar_y=solar_y, log=log, ...)
    }

    if (mesh) {
        mesh <- get_mesh(X, Y, Z, seis.DF, cygA_stds)
        make_plots(mesh_plot, paste0(paste(Z, Y, X, sep="_"), "-mesh"),
                   filepath=filepath, thin.hack=1,
                   mar=c(3, 4, 1, 0),
                   mesh=mesh, X=X, Z=Z, xlab=xlab, ylab=ylab,
                   solar_x=solar_x, solar_y=solar_y, log=log, ...)
    }
}

plot_lamb_brunt <- function(DF, ...,
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    lambs <- NULL
    for (ell in 1:3) {
        S <- DF$csound**2 / (DF$radius/max(DF$radius) * solar_radius)**2
        lamb <- data.frame(10**6/(2*pi) * sqrt(ell*(ell+1) * S))
        names(lamb) <- ell
        lambs <- if (is.null(lambs)) lamb else cbind(lambs, lamb)
    }
    brunt_N2 <- DF$brunt_N2
    stable <- which(brunt_N2 >= 0)
    indices <- if (any(diff(stable)!=1)) {
        stable[which(diff(stable)!=1)[1]+1] : max(stable)
    } else {
        stable
    }
    brunt <- 10**6/(2*pi) * sqrt(brunt_N2[indices])
    plot(DF$radius[indices]/max(DF$radius), brunt,
        axes=0, tcl=0, type='l',
        log='y', #xaxs='i', #cex=0.5,
        ylim=range(brunt, lambs),#, 0.01, 10**6),
        #ylim=c(1, 10**5),
        xlim=range(DF$radius/max(DF$radius)),#range(0.01, max(DF$radius)),
        #xlim=c(0.088, 0.09),
        #xlim=c(0.0815, 0.0825),
        xlab=expression("Radius" ~ r/R["*"]),
        ylab=expression("Frequency" ~ nu/mu*Hz))
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), las=1,
        mgp=mgp, cex.axis=text.cex)
    for (ell in 1:3) lines(DF$radius/max(DF$radius), lambs[[ell]], lty=ell+1)
    legend("topright", lty=c(1,2), #pch=c(1, NA),
        legend=c("Brunt", "Lamb"), bty='n')
}


################################################################################
### Scaling relations ##########################################################
################################################################################

## R/R_sun = nu_max/nu_max_sun (Dnu/Dnu_sun)^-2 (Teff/Teff_sun)^(1/2)
## M/M_sun = (R/R_sun)^3 (Dnu/Dnu_sun)^2
scaling_nu_max <- function(R, M, Teff, Teff_sun=5772, nu_max_sun=3090) {
    M * nu_max_sun / ( R**2 * sqrt(Teff/Teff_sun) )
}

scaling_nu_max_Viani <- function(R, M, Teff, mu,
        Teff_sun=5772, nu_max_sun=3090, mu_sun=1.2546895396570101) {
    M * nu_max_sun * sqrt(mu / mu_sun) / ( R**2 * sqrt(Teff/Teff_sun) )
}

scaling_Delta_nu_Guggenberger <- function(Teff=5777.739, FeH=0) {
    A <- 0.64 * FeH + 1.78
    lambda <- -0.55 * FeH + 1.23
    omega <- 22.21
    phi <- 0.48 * FeH + 0.12
    B <- 0.66 * FeH + 134.92
    A * exp( lambda * Teff / 10**4 ) * cos( omega * Teff / 10**4 + phi ) + B
}

scaling_Delta_nu <- function(R, M, Delta_nu_Sun=135) {
    Delta_nu_Sun * M / R**3
}

## Delta_nu = 1 / ( 2 * integral_0^R [ dr/c ] )
## c.f. Aerts et al. 2010, Asteroseismology, eq. 3.217
asymptotic_Delta_nu <- function(radius, csound) {
    invisible(library(Bolstad))
    1 / (2*sintegral(radius*solar_radius, 1/(csound*10**6))$value)
}


################################################################################
### B-splines ##################################################################
################################################################################
# x is the location
# i is the ith knot 
# j is the degree of the spline 
# ks are the knot locations 
B <- function(x, i, j, ks) {
    if (j == 0) {
        if (ks[i] <= x && x < ks[i+1]) 1 else 0
    } else {
        den.1 <- ks[i+j] - ks[i]
        den.2 <- ks[i+j+1] - ks[i+1]
        exp.1 <- if (den.1 == 0) 0 else (x-ks[i])/den.1 * B(x, i, j-1, ks)
        exp.2 <- if (den.2 == 0) 0 else (ks[i+j+1]-x)/den.2 * B(x, i+1, j-1, ks)
        exp.1 + exp.2
    }
}

get_knots <- function(x, num_knots, degree=4) {
    #x <- x[x>0]
    n <- num_knots-(2*degree)
    c(rep(min(x), degree), 
      seq(min(x), max(x), length=n),
      rep(max(x), degree))
}

dB_dx <- function(x, i, j, ks) {
    den.1 <- ks[i+j] - ks[i+1]
    den.2 <- ks[i+j-1] - ks[i]
    exp.1 <- if (den.1 == 0) 0 else -B(x, i+1, j-1, ks) / den.1
    exp.2 <- if (den.2 == 0) 0 else B(x, i, j-1, ks) / den.2
    (j-1) * ( exp.1 + exp.2 )
}

d.n_B_dx.n <- function(x, i, j, ks, n=2) {
    if (n==1) return(dB_dx(x, i, j, ks))
    den.1 <- ks[i+j] - ks[i+1]
    den.2 <- ks[i+j-1] - ks[i]
    exp.1 <- if (den.1 == 0) 0 else -d.n_B_dx.n(x, i+1, j-1, ks, n-1) / den.1
    exp.2 <- if (den.2 == 0) 0 else d.n_B_dx.n(x, i, j-1, ks, n-1) / den.2
    (j-1) * ( exp.1 + exp.2 )
}


################################################################################
### MISC #######################################################################
################################################################################

cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

shadowtext <- function(x, y=NULL, labels, col='black', bg='white', 
                       theta= seq(0, 2*pi, length.out=50), r=0.2, ... ) {
    
    xy <- xy.coords(x,y)
    xo <- r*strwidth('A')
    yo <- r*strheight('A')
    
    # draw background text with small shift in x and y in background colour
    for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
    }
    # draw actual text in exact xy position in foreground colour
    text(xy$x, xy$y, labels, col=col, ... )
}
