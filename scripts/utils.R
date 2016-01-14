#### Utility functions for forward & inverse asteroseismology 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Plotting values
utils.mgp <<- c(2, 0.25, 0)
hack.mgp <- c(2, 0.5, 0)

font <- "Palatino"
png_res <- 400
cex.paper <- 1
cex.slides <- 1.3
cex.hack <- 1.4

paper_pdf_width <- 6.97522 # inches
paper_pdf_height <- 4.17309

latex_pt_per_in <- 5 * 72.27
paper_png_width <- paper_pdf_width * latex_pt_per_in
paper_png_height <- paper_pdf_height * latex_pt_per_in

slides_pdf_width <- 6.22665 
slides_pdf_height <- 4.1511

slides_png_width <- slides_pdf_width * latex_pt_per_in
slides_png_height <- slides_pdf_height * latex_pt_per_in

## Make the same plot as a pdf and a png suitable for papers and for slides
# takes a plotting function plot_f that calls `plot` 
make_plots <- function(plot_f, filename, ..., 
        filepath='plots', mar=c(3, 4, 1, 1), mgp=utils.mgp, 
        wide=TRUE, thin=TRUE, 
        make_png=TRUE, make_pdf=TRUE, 
        paper=TRUE, slides=TRUE,
        thin.hack=FALSE) {
    
    if (paper & make_pdf & wide) {
        directory <- file.path(filepath, 'paper', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '.pdf')),
                  width=paper_pdf_width, height=paper_pdf_height, 
                  family=font)
        par(mar=mar, mgp=utils.mgp, cex.lab=cex.paper, family=font)
        plot_f(text.cex=cex.paper, ...)
        dev.off()
    }
    
    if (paper & make_pdf & thin) {
        directory <- file.path(filepath, 'paper', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        if (thin.hack) {
            cairo_pdf(file.path(directory, paste0(filename, '-thin.pdf')),
                      width=paper_pdf_width, height=paper_pdf_height, 
                      family=font)
            par(mar=mar, mgp=hack.mgp, cex.lab=cex.hack, family=font)
            plot_f(text.cex=cex.hack, mgp=hack.mgp, ...)
        } else {
            cairo_pdf(file.path(directory, paste0(filename, '-thin.pdf')),
                      width=paper_pdf_width/2, height=paper_pdf_height, 
                      family=font)
            par(mar=mar, mgp=utils.mgp, cex.lab=cex.paper, family=font)
            plot_f(text.cex=cex.paper, ...)
        }
        dev.off()
    }
    
    if (paper & make_png & wide) {
        directory <- file.path(filepath, 'paper', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '.png')),
                  width=paper_png_width, height=paper_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=utils.mgp, cex.lab=cex.paper, family=font)
        plot_f(text.cex=cex.paper, ...)
        dev.off()
    }
    
    if (paper & make_png & thin) {
        directory <- file.path(filepath, 'paper', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        if (thin.hack) {
            png(file.path(directory, paste0(filename, '-thin.png')),
                      width=paper_png_width, height=paper_png_height, 
                      family=font, res=png_res, type='cairo')
            par(mar=mar, mgp=hack.mgp, cex.lab=cex.hack, family=font)
            plot_f(text.cex=cex.hack, ...)
        } else {
            png(file.path(directory, paste0(filename, '-thin.png')),
                      width=paper_png_width/2, height=paper_png_height, 
                      family=font, res=png_res, type='cairo')
            par(mar=mar, mgp=utils.mgp, cex.lab=cex.paper, family=font)
            plot_f(text.cex=cex.paper, ...)
        }
        dev.off()
    }
    
    if (slides & make_pdf & wide) {
        directory <- file.path(filepath, 'slides', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '-slides.pdf')),
                  width=slides_pdf_width, height=slides_pdf_height, 
                  family=font)
        par(mar=mar, mgp=utils.mgp, cex.lab=cex.slides, family=font)
        plot_f(text.cex=cex.slides, ...)
        dev.off()
    }
    
    if (slides & make_pdf & thin) {
        directory <- file.path(filepath, 'slides', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        if (thin.hack) {
            cairo_pdf(file.path(directory, 
                      paste0(filename, '-slides-thin.pdf')),
                      width=slides_pdf_width, height=slides_pdf_height, 
                      family=font)
            par(mar=mar, mgp=hack.mgp, cex.lab=cex.slides*cex.hack, 
                family=font)
            plot_f(text.cex=cex.slides*cex.hack, ...)
        } else {
            cairo_pdf(file.path(directory, 
                      paste0(filename, '-slides-thin.pdf')),
                      width=slides_pdf_width/2, height=slides_pdf_height, 
                      family=font)
            par(mar=mar, mgp=utils.mgp, cex.lab=cex.slides, family=font)
            plot_f(text.cex=cex.slides, ...)
        }
        dev.off()
    }
    
    if (slides & make_png & wide) {
        directory <- file.path(filepath, 'slides', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '-slides.png')),
                  width=slides_png_width, height=slides_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=utils.mgp, cex.lab=cex.slides, family=font)
        plot_f(text.cex=cex.slides, ...)
        dev.off()
    }
    
    if (slides & make_png & thin) {
        directory <- file.path(filepath, 'slides', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        if (thin.hack) {
            png(file.path(directory, paste0(filename, '-slides-thin.png')), 
                      width=slides_png_width, height=slides_png_height, 
                      family=font, res=png_res, type='cairo')
            par(mar=mar, mgp=hack.mgp, cex.lab=cex.slides*cex.hack, 
                family=font)
            plot_f(text.cex=cex.slides*cex.hack, ...)
        } else {
            png(file.path(directory, paste0(filename, '-slides-thin.png')), 
                      width=slides_png_width/2, height=slides_png_height, 
                      family=font, res=png_res, type='cairo')
            par(mar=mar, mgp=utils.mgp, cex.lab=cex.slides, family=font)
            plot_f(text.cex=cex.slides, ...)
        }
        dev.off()
    }
}

##############################################################################
### Metadata for the grid ####################################################
##############################################################################

get_label <- function(symbol) as.expression(bquote(
        .(seis.names[[symbol]])
      ~ .(seis.labs[[symbol]])
      * .(seis.units[[symbol]])
))

get_label_nameless <- function(symbol) as.expression(bquote(
        .(seis.labs[[symbol]])
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
  dnu02_median   = "Small separation", 
  dnu02_slope    = "", 
  r_sep02_median = "", 
  r_sep02_slope  = "", 
  r_avg01_median = "", 
  r_avg01_slope  = "", 
  dnu13_median   = "Small separation", 
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
  H              = bquote(X), 
  He             = bquote(Y), 
  Hc             = bquote(X[c]), 
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
  log_g          = bquote("/dex (cgs)"), 
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
  H              = "X", 
  He             = "Y", 
  Hc             = "$X_c$", 
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

## Scatter plot function
scatter_plot <- function(seis.DF, X, Y, Z, combos, col.pal, 
        xlim, ylim, xlab, ylab, ..., 
        solar_x=NA, solar_y=NA, text.cex=cex.paper, mgp=utils.mgp) { 
    Z_max <- max(seis.DF[[Z]])
    Z_min <- min(seis.DF[[Z]])
    for (simulation_i in 1:nrow(combos)) {
        track <- merge(seis.DF, combos[simulation_i,])
        color <- col.pal[floor((track[[Z]]-Z_min)/(Z_max-Z_min)*
                                   (length(col.pal)-1))+1]
        use_line <- length(unique(color))==1
        relation <- track[[Y]] ~ track[[X]]
        if (simulation_i == 1) {
            plot(relation, type=ifelse(use_line, 'l', 'p'),
                 pch=20, axes=FALSE, col=color, cex=0.1, tcl=0, 
                 xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab)
            if (!is.na(solar_x) && !is.na(solar_y)) {
                abline(v=solar_x, lty=3, col='black')
                abline(h=solar_y, lty=3, col='black')
            }
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=mgp, cex.axis=text.cex)
        } else {
            if (use_line) lines(relation, col=color[1])
            else points(relation, col=color, pch=20, cex=0.1)
        }
    }
    
    points(solar_x, solar_y, pch=1, cex=1)
    points(solar_x, solar_y, pch=20, cex=0.1)
    X_range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 2), 
                 col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(get_label(Z), 4, line=5, cex=text.cex)
}

## A basic normalization function
normalize <- function(x) (x-min(x))/(max(x)-min(x))

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
mesh_plot <- function(mesh, X, Z, xlab, ylab, ..., 
        solar_x=NA, solar_y=NA, text.cex=cex.paper, mgp=utils.mgp) {
    filled.contour(mesh,
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
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), 
                    mgp=mgp, cex.axis=text.cex)
        },
        plot.title={
            title(xlab=xlab, cex.lab=text.cex, line=2)
            title(ylab=ylab, cex.lab=text.cex, line=2)
        })
}

## Calls scatter_plot and mesh_plot 
# uses globals seis.DF, col.pal, and combos 
scatter_mesh <- function(X, Y, Z, filepath=file.path("plots", "mesh")) {
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
    
    make_plots(scatter_plot, paste0(paste(Z, Y, X, sep="_"), "-scatter"),
               filepath=filepath, thin.hack=1,
               mar=c(3, 4, 1, 6), 
               seis.DF=seis.DF, X=X, Y=Y, Z=Z, combos=combos, 
               col.pal=col.pal, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
               solar_x=solar_x, solar_y=solar_y)
    
    mesh <- get_mesh(X, Y, Z, seis.DF, cygA_stds)
    make_plots(mesh_plot, paste0(paste(Z, Y, X, sep="_"), "-mesh"),
               filepath=filepath, thin.hack=1,
               mar=c(3, 4, 1, 0), 
               mesh=mesh, X=X, Z=Z, xlab=xlab, ylab=ylab, 
               solar_x=solar_x, solar_y=solar_y)
}

