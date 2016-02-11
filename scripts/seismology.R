#### Seismological calculations for stellar observations and models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path(dirname(sys.frame(1)$ofile), 'utils.R'))

library(matrixStats)
library(magicaxis)
library(RColorBrewer)

#dnu.cl <- brewer.pal(4, "BrBG")
dnu.cl <- c("#ca0020", "#f4a582", "#0571b0", "#800080")
#c("#ca0020", "#f4a582", "#92c5de", "#0571b0")

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
        radials <- freqs[freqs$l==l_mode & freqs$n>0,]
        duplicates <- duplicated(radials$n)
        if (any(duplicates)) {
            print("Duplicated (l,n) combination")
            print(radials[duplicates,])
            return(NULL)
        }
    }
    
    # get averages for Dnu, dnu, r01, r02, r10, r13
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

## Calculate averages of things like f = dnu, Dnu, r_sep, r_avg
# DF is the where the result will be stored
# freqs are a data frame with columns l, n, nu
# l_degs are the l's for which this calculation should be made 
# nu_max is the center of the gaussian
# make a plot with filename 'outf' if outf != FALSE
avg <- function(f, DF, freqs, l_degs, nu_max, outf=FALSE, ...) {
    sep_name <- deparse(substitute(f))
    a <- c() # contains the computed quantity (e.g. large freq separations)
    b <- c() # contains frequencies of the base mode
    pchs <- c() # if there's more than one l, get different symbols for each
    p_modes <- freqs[freqs$n >= 1,] 
    for (l_deg in l_degs) {
        ell <- p_modes[p_modes$l==l_deg,]
        if (nrow(ell) == 0) next
        vals <- sapply(unique(ell$n), function(n) f(l_deg, n, freqs))
        not.nan <- complete.cases(vals)
        a <- c(a, vals[not.nan])
        b <- c(b, ell$nu[not.nan])
        if (outf != FALSE) pchs = c(pchs, rep(l_deg+1, sum(not.nan)))
    }
    
    # need 3 points to make something reasonable 
    if (length(a)<=2) {
        #DF[paste0(sep_name, "_median")] <- NA
        #DF[paste0(sep_name, "_slope")] <- NA
        return(DF)
    }
    
    # build expression for y label of plot
    if (outf != FALSE) {
        ylab <- if (sep_name == 'Dnu' && length(l_degs) > 1) 
               bquote(Delta*nu)
           else if (sep_name == 'Dnu')   
               bquote(Delta*nu[.(l_degs)])
           else if (sep_name == 'dnu')   
               bquote(delta*nu[.(l_degs)*','*.(l_degs+2)])
           else if (sep_name == 'r_sep') bquote(r[.(l_degs)*','*.(l_degs+2)])
           else if (sep_name == 'r_avg') bquote(r[.(l_degs)*','*.(1-l_degs)])
        ylab <- bquote(.(ylab) / mu*Hz)
    }
    
    sep_name <- if (sep_name == 'Dnu' && length(l_degs) > 1) paste0(sep_name)
       else if (sep_name == 'Dnu')   paste0(sep_name, l_degs)
       else if (sep_name == 'dnu')   paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_sep') paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_avg') paste0(sep_name, l_degs, 1-l_degs)
    
    fwhm <- (0.66*nu_max**0.88)/(2*sqrt(2*log(2)))
    gaussian_env <- dnorm(b, nu_max, fwhm)
    w.median <- weightedMedian(a, gaussian_env)
    DF[paste0(sep_name, "_median")] <- w.median
    fit <- lm(a~b, weights=gaussian_env)
    DF[paste0(sep_name, "_slope")] <- coef(fit)[2]
    
    if (outf != FALSE) make_plots(seismology_plot, 
        paste0(outf, '-', sep_name), 
        a=a, b=b, fit=fit, gaussian_env=gaussian_env, 
        w.median=w.median, nu_max=nu_max, l_degs=l_degs, 
        ylab=ylab, dnu.cl=dnu.cl, pchs=pchs, sep_name=sep_name, 
        freqs=freqs,
        ...)
    
    DF
}

## Separation: just the difference between two frequencies 
# nu_{l1,n1} - nu_{l2,n2} 
# the difference must be positive, otherwise it returns NA 
separation <- function(first_l, first_n, second_l, second_n, DF) {
    first <- DF$l == first_l & DF$n == first_n
    second <- DF$l == second_l & DF$n == second_n
    if (sum(first) == 1 && sum(second) == 1) { # check that it's unique
        difference <- DF[first,]$nu - DF[second,]$nu
        if (difference >= 0) return(difference)
    }
    return(NA)
}

## Five point averages 
#dd_01= 1/8( nu_[n-1,0] - 4*nu_[n-1,1] + 6*nu_[n,0] - 4*nu[n,  1] + nu_[n+1,0] )
#dd_10=-1/8( nu_[n-1,1] - 4*nu_[n,  0] + 6*nu_[n,1] - 4*nu[n+1,0] + nu_[n+1,1] )
dd <- function(l0, l1, n, DF) {
    p_modes <- DF$n>0
    ell.0 <- DF[DF$l==0 & p_modes,]
    ell.1 <- DF[DF$l==1 & p_modes,]
    ns <- DF[DF$n>=(n-1) || DF$n<=(n+1)]
    n. <- ns[ns$n==n,]#DF[DF$n==n,]
    n.minus.one <- ns[ns$n==(n-1),]#DF[DF$n==n-1,]
    n.plus.one <- ns[ns$n==(n+1),]#DF[DF$n==n+1,]
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
dnu <- function(l, n, DF) separation(l, n, l+2, n-1, DF)
Dnu <- function(l, n, DF) separation(l, n, l, n-1, DF)
r_sep <- function(l, n, DF) dnu(l, n, DF) / Dnu(1-l, n+l, DF)
r_avg <- function(l, n, DF) dd(l, 1-l, n, DF) / Dnu(1-l, n+l, DF)

## Plot Dnu, dnu, r02, ... 
seismology_plot <- function(a, b, fit, gaussian_env, w.median, 
        nu_max, l_degs, ylab, dnu.cl, pchs, freqs, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    if (length(l_degs)==1)
        col.pal <- colorRampPalette(c(dnu.cl[1], dnu.cl[3]))(1001)[1+1000*
            normalize(gaussian_env)]
    plot(a~b, axes=FALSE, tck=0, xaxs='i',
         cex=1.5 * gaussian_env/max(gaussian_env), 
         ylab=as.expression(ylab), 
         xlab=expression("Frequency" ~ nu / mu*Hz), 
         xlim=range(freqs$nu), 
         #ylim=quantile(a, c(0.001, 0.999)), 
         ylim=range(w.median, coef(fit)[1], 2*w.median-coef(fit)[1], a), 
         col=if (length(l_degs)==1) col.pal else dnu.cl[pchs], 
         pch=if (length(l_degs)==1) 1 else pchs)
    abline(fit, lty=2)
    abline(v=nu_max, lty=3)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    if (length(l_degs)>1)
        legend("bottomright", pch=l_degs+1, col=dnu.cl, cex=text.cex, bty="n",
               legend=paste0("\u2113=", l_degs), horiz=1)
}

