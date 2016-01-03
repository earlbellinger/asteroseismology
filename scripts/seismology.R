#### Seismological calculations for stellar observations and models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path(dirname(sys.frame(1)$ofile), 'utils.R'))

library(matrixStats)
library(magicaxis)
library(RColorBrewer)

dnu.cl <- brewer.pal(4, "BrBG")

## Separation: just the difference between two frequencies 
# nu_{l1,n1} - nu_{l2,n2}
separation <- function(first_l, first_n, second_l, second_n, df) {
    first <- df$l == first_l & df$n == first_n
    second <- df$l == second_l & df$n == second_n
    if (sum(first) == 1 && sum(second) == 1) { # check that it's unique
        difference <- df[first,]$nu - df[second,]$nu
        if (difference < 0) return(NA)
        return(difference)
    }
    return(NA)
}

## Five point averages 
#dd_01= 1/8( nu_[n-1,0] - 4*nu_[n-1,1] + 6*nu_[n,0] - 4*nu[n,  1] + nu_[n+1,0] )
#dd_10=-1/8( nu_[n-1,1] - 4*nu_[n,  0] + 6*nu_[n,1] - 4*nu[n+1,0] + nu_[n+1,1] )
dd <- function(l0, l1, n, df) {
    ell.0 <- df[df$l==0 & df$n>0,]
    ell.1 <- df[df$l==1 & df$n>0,]
    n. <- df[df$n==n,]
    n.minus.one <- df[df$n==n-1,]
    n.plus.one <- df[df$n==n+1,]
    val <- if (l0 == 0 && l1 == 1) { ## dd_01
        ( merge(n.minus.one, ell.0)$nu -
        4*merge(n.minus.one, ell.1)$nu +
        6*merge(n., ell.0)$nu -
        4*merge(n., ell.1)$nu +
          merge(n.plus.one, ell.0)$nu )/8
    } else if (l1 == 0 && l0 == 1) { ## dd_10
        -( merge(n.minus.one, ell.1)$nu -
         4*merge(n., ell.0)$nu +
         6*merge(n., ell.1)$nu -
         4*merge(n.plus.one, ell.0)$nu +
           merge(n.plus.one, ell.1)$nu )/8
    } else NA
    if (length(val) == 0) NA
    else val
}

## Separations and ratios
dnu <- function(l, n, df) separation(l, n, l+2, n-1, df)
Dnu <- function(l, n, df) separation(l, n, l, n-1, df)
r_sep <- function(l, n, df) dnu(l, n, df) / Dnu(1-l, n+l, df)
r_avg <- function(l, n, df) dd(l, 1-l, n, df) / Dnu(1-l, n+l, df)

## Plot Dnu, dnu, r02, ... 
seismology_plot <- function(text.cex, 
        a, b, fit, gaussian_env, w.median, nu_max, l_degs, 
        ylab, dnu.cl, pchs, ...) {
    
    plot(a~b, tck=0, ylab=as.expression(ylab), 
         cex=2*gaussian_env/max(gaussian_env), 
         ylim=range(a, w.median, coef(fit)[1], 2*w.median-coef(fit)[1]), 
         col=if (length(l_degs)==1) 1 else dnu.cl[pchs], 
         pch=if (length(l_degs)==1) 1 else pchs, 
         xlab=expression("Frequency" ~ nu / mu*Hz))
    abline(fit, lty=2)
    abline(v=nu_max, lty=3)
    magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
    if (length(l_degs)>1)
        legend("topleft", pch=l_degs+1, col=dnu.cl, cex=text.cex, #horiz=1,
               #ncol=length(l_degs), #bty="n",
               legend=paste0("\u2113=", l_degs))
}

## Calculate averages of things like f = dnu, Dnu, r_sep, r_avg
# df is the where the result will be stored
# freqs are a data frame with columns l, n, nu
# l_degs are the l's for which this calculation should be made 
# nu_max is the center of the gaussian
# make a plot with filename 'outf' if outf != FALSE
avg <- function(f, df, freqs, l_degs, nu_max, outf=FALSE, ...) {
    sep_name <- deparse(substitute(f))
    a <- c() # contains the computed quantity (e.g. large freq separations)
    b <- c() # contains frequencies of the base mode
    pchs <- c() # if there's more than one l, get different symbols for each
    for (l_deg in l_degs) {
        ell <- freqs[freqs$n > 1 & freqs$l==l_deg,]
        vals <- sapply(unique(ell$n), function(n) f(l_deg, n, freqs))
        not.nan <- complete.cases(vals)
        a <- c(a, vals[not.nan])
        b <- c(b, ell$nu[not.nan])
        pchs = c(pchs, rep(l_deg+1, sum(not.nan)))
    }
    
    if (length(a)==0) return(NA)
    
    # build expression for y label of plot
    ylab <- if (sep_name == 'Dnu' && length(l_degs) > 1) bquote(Delta*nu)
       else if (sep_name == 'Dnu')   bquote(Delta*nu[.(l_degs)])
       else if (sep_name == 'dnu')   bquote(delta*nu[.(l_degs)*','*.(l_degs+2)])
       else if (sep_name == 'r_sep') bquote(r[.(l_degs)*','*.(l_degs+2)])
       else if (sep_name == 'r_avg') bquote(r[.(l_degs)*','*.(1-l_degs)])
    ylab <- bquote(.(ylab) / mu*Hz)
    
    sep_name <- if (sep_name == 'Dnu' && length(l_degs) > 1) paste0(sep_name)
       else if (sep_name == 'Dnu')   paste0(sep_name, l_degs)
       else if (sep_name == 'dnu')   paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_sep') paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_avg') paste0(sep_name, l_degs, 1-l_degs)
    
    fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
    gaussian_env <- dnorm(b, nu_max, fwhm)
    w.median <- weightedMedian(a, gaussian_env)
    df[paste0(sep_name, "_median")] <- w.median
    fit <- lm(a~b, weights=gaussian_env)
    df[paste0(sep_name, "_slope")] <- coef(fit)[2]
    
    print(outf)
    if (outf != FALSE) make_plots(seismology_plot, 
        paste0(outf, '-', sep_name), 
        a, b, fit, gaussian_env, w.median, nu_max, l_degs, 
        ylab, dnu.cl, pchs, ...)
    
    df
}

## build a data frame containing seismological calculations
# freqs is a data frame with l, n, and nu
# nu_max is the frequency of maximum oscillation power
# acoustic_cutoff is the truncation frequency 
# outf is the filename that plots should have (None for no plot)
seismology <- function(freqs, nu_max, ..., acoustic_cutoff=Inf, outf=FALSE) {
    if (nrow(freqs) == 0) {
        print("No frequencies found")
        return(NULL)
    }
    freqs <- unique(freqs[complete.cases(freqs) & freqs$nu < acoustic_cutoff,])
    
    # all (l,n) combinations should be unique; discard otherwise
    for (l_mode in unique(freqs$l)) {
        if (any(duplicated(freqs[freqs$l==l_mode & freqs$n>0,]$n))) {
            print("Duplicated (l,n) combination")
            return(NULL)
        }
    }
    
    seis.DF <- avg(Dnu, NULL, freqs, sort(unique(freqs$l)), nu_max, outf, ...)
    seis.DF <- avg(Dnu, seis.DF, freqs, 0, nu_max, outf, ...)
    for (l_deg in 0:1) {
        if (l_deg==0 || 3 %in% freqs$l) { # some stars only have l=0,1,2
            seis.DF <- avg(dnu, seis.DF, freqs, l_deg, nu_max, outf, ...)
            seis.DF <- avg(r_sep, seis.DF, freqs, l_deg, nu_max, outf, ...)
        }
        seis.DF <- avg(r_avg, seis.DF, freqs, l_deg, nu_max, outf, ...)
    }
    return(seis.DF)
}
