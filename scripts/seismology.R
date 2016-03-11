#### Seismological calculations for stellar observations and models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path(dirname(sys.frame(1)$ofile), 'utils.R'))

library(matrixStats)
library(magicaxis)
library(RColorBrewer)

fwhm_conversion <- (2*sqrt(2*log(2)))

#dnu.cl <- brewer.pal(4, "BrBG")
dnu.cl <- c("#ca0020", "#f4a582", "#0571b0", "#800080")
#c("#ca0020", "#f4a582", "#92c5de", "#0571b0")

## build a data frame containing seismological calculations
# freqs is a data frame with l, n, and nu
# nu_max is the frequency of maximum oscillation power
# acoustic_cutoff is the truncation frequency 
# outf is the filename that plots should have (None for no plot)
seismology <- function(freqs, nu_max, 
        ..., min_points_per_deg=15, acoustic_cutoff=Inf, outf=FALSE) {
    if (nrow(freqs) == 0) {
        print("No frequencies found")
        return(NULL)
    }
    #freqs <- freqs[,1:3]
    freqs <- unique(freqs[complete.cases(freqs) & freqs$nu < acoustic_cutoff,])
    
    # get l=0
    seis.DF <- avg(Dnu, freqs, nu_max=nu_max, outf=outf, 
        acoustic_cutoff=acoustic_cutoff, ...)
    
    # make echelle
    if (outf != FALSE) make_plots(echelle_plot, paste0(outf, '-echelle'), 
        freqs=freqs, large_sep=seis.DF[[1]], ...)
    
    # all (l,n) combinations should be unique with no mixing; discard otherwise
    #for (l_mode in unique(freqs$l)) {
    #    radials <- freqs[freqs$l==l_mode & freqs$n>=0,]
    #    duplicates <- duplicated(radials$n)
    #    if (any(duplicates)) {
    #        print("Duplicated (l,n) combination")
    #        print(radials[duplicates,])
    #        #while (duplicated(radials$n)) 
    #        #    radials <- radials[-which(duplicated(radials$n)),]
    #        #freqs <- rbind(freqs[freqs$l!=l_mode | freqs$n<0,], radials)
    #        return(NULL)
    #    }
        #if ('inertia' %in% names(radials) && any(diff(radials$inertia) > 0)) {
         #   #||
         #       #abs(diff(radials$inertia))>100*radials$inertia[-1]) 
         #       #{
         #   print("Mixed modes")
         #   #print(radials)
         #   #while (any(diff(radials$inertia) > 0)) 
         #   #    radials <- radials[c(1, 
         #   #        1+which(diff(radials$inertia) <= 0)),]
         #   #freqs <- rbind(freqs[freqs$l!=l_mode | freqs$n<0,], radials)
         #   return(NULL)
        #}
        #if (nrow(radials) <= min_points_per_deg) {
        #    print("Too few points")
        #    print(radials)
        #    return(NULL)
        #}
    #}
    
    # get averages for Dnu, dnu, r01, r02, r10, r13
    #seis.DF <- avg(Dnu, seis.DF, freqs, sort(unique(freqs$l)), 
    #    nu_max, outf, ...)
    for (l_deg in 0:1) {
        if (l_deg==0 || 3 %in% freqs$l) { # some stars only have l=0,1,2
            seis.DF <- avg(dnu, freqs, l_deg, seis.DF, nu_max, outf, ...)
            seis.DF <- avg(r_sep, freqs, l_deg, seis.DF, nu_max, outf, ...)
        }
        seis.DF <- avg(r_avg, freqs, l_deg, seis.DF, nu_max, outf, ...)
    }
    return(seis.DF)
}

## Calculate averages of things like separator = dnu, Dnu, r_sep, r_avg
# DF is the where the result will be stored
# freqs are a data frame with columns l, n, nu
# l_degs are the l's for which this calculation should be made 
# nu_max is the center of the gaussian
# make a plot with filename 'outf' if outf != FALSE
avg <- function(separator, freqs, 
        l_degs=0, DF=NULL, nu_max=NA, acoustic_cutoff=Inf, outf=FALSE, ...) {
    sep_name <- deparse(substitute(separator))
    #print(sep_name)
    seps <- c() # contains the computed quantity (e.g. large freq separations)
    nus <- c() # contains frequencies of the base mode
    pchs <- c() # if there's more than one l, get different symbols for each
    p_modes <- freqs[freqs$n >= 1,] 
    for (l_deg in l_degs) {
        ell <- p_modes[p_modes$l==l_deg,]
        ell <- ell[!duplicated(ell$n) & !duplicated(ell$n, fromLast=T),]
        if (nrow(ell) == 0) {
            print(paste0("No ell=", ell))
            #print(ell)
            next
        }
        vals <- sapply(unique(ell$n), function(n) separator(l_deg, n, freqs))
        not.nan <- complete.cases(vals)
        seps <- c(seps, vals[not.nan])
        nus <- c(nus, ell$nu[not.nan])
        if (outf != FALSE) pchs = c(pchs, rep(l_deg+1, sum(not.nan)))
    }
    
    # need 3 points to make something reasonable 
    if (length(seps)<=2) {
        print(paste("Too few points for", sep_name))
        print(seps)
        #DF[paste0(sep_name, "_median")] <- NA
        #DF[paste0(sep_name, "_slope")] <- NA
        return(DF)
    }
    
    # build expression for y label of plot
    if (outf != FALSE) {
        ylab <- if (sep_name == 'Dnu' && length(l_degs) > 1) 
               bquote(Delta*nu / mu*Hz)
           else if (sep_name == 'Dnu') 
               bquote(Delta*nu[.(l_degs)] / mu*Hz)
           else if (sep_name == 'dnu')
               bquote(delta*nu[.(l_degs)*','*.(l_degs+2)] / mu*Hz)
           else if (sep_name == 'r_sep') 
               bquote(r[.(l_degs)*','*.(l_degs+2)])
           else if (sep_name == 'r_avg') 
               bquote(r[.(l_degs)*','*.(1-l_degs)])
       ylab <- as.expression(ylab)
    }
    
    sep_name <- if (sep_name == 'Dnu' && length(l_degs) > 1) paste0(sep_name)
       else if (sep_name == 'Dnu')   paste0(sep_name, l_degs)
       else if (sep_name == 'dnu')   paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_sep') paste0(sep_name, l_degs, l_degs+2)
       else if (sep_name == 'r_avg') paste0(sep_name, l_degs, 1-l_degs)
    
    fwhm <- (0.66*nu_max**0.88)/fwhm_conversion
    if (!is.na(nu_max)) {
        gaussian_env <- dnorm(nus, nu_max, fwhm)
        w.median <- weightedMedian(seps, gaussian_env)
        DF[paste0(sep_name, "_median")] <- w.median
        #fit <- lm(seps~nus, weights=gaussian_env)
    } else {
        DF[paste0(sep_name, "_median")] <- median(seps)
        #fit <- lm(seps~nus)
    }
    #DF[paste0(sep_name, "_slope")] <- coef(fit)[2]
    
    if (outf != FALSE) make_plots(seismology_plot, 
        paste0(outf, '-', sep_name), 
        seps=seps, nus=nus, #fit=fit, 
        gaussian_env=gaussian_env, 
        w.median=w.median, nu_max=nu_max, l_degs=l_degs, 
        ylab=ylab, dnu.cl=dnu.cl, pchs=pchs, sep_name=sep_name, 
        freqs=freqs,
        ...)
    
    DF
}

## Separation: just the difference between two frequencies 
# nu_{l1,n1} - nu_{l2,n2} 
# the difference must uniquely exist and be positive, otherwise it returns NA 
separation <- function(first_l, first_n, second_l, second_n, DF, use_n=T) {
    first <- DF[DF$l == first_l & DF$n == first_n,]
    if (nrow(first) != 1) return(NA)
    second <- if (use_n) {
        DF[DF$l == second_l & DF$n == second_n,]
    } else {
        n.diff <- first_n - second_n
        second.modes <- DF[DF$l == second_l & DF$nu < first$nu,]
        if (all(is.na(second.modes))) return (NA)
        second.modes[order(second.modes$nu, decreasing=T),][n.diff,]
    }
    if (nrow(second) != 1) return(NA)
    difference <- first$nu - second$nu
    if (difference <= 0) return(NA)
    return(difference)
}

## Five point averages 
#dd_01= 1/8( nu_[n-1,0] - 4*nu_[n-1,1] + 6*nu_[n,0] - 4*nu[n,  1] + nu_[n+1,0] )
#dd_10=-1/8( nu_[n-1,1] - 4*nu_[n,  0] + 6*nu_[n,1] - 4*nu[n+1,0] + nu_[n+1,1] )
# with the subscript of nu being (n, l)
dd <- function(l0, l1, n, DF, use_n=T) {
    if (l0 == 0 && l1 != 1 || l0 == 1 && l1 != 0) return(NA)
    
    p_modes <- DF$n>=0
    ell.0 <- DF[DF$l==0 & p_modes,]
    ell.1 <- DF[DF$l==1 & p_modes,]
    if (use_n) {
        ns <- DF[DF$n>=(n-1) & DF$n<=(n+1),]
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
        
        if (length(val) != 1) {
            NA
        } else {
           val
        }
    
    } else {
        ref.l <- if (l0 == 0 && l1 == 1) ell.0 else ell.1
        ref.nu <- ref.l[ref.l$n == n,]$nu
        if (length(ref.nu) != 1) return(NA)
        lesser <- ref.l$nu < ref.nu
        greater <- ref.l$nu > ref.nu
        if (!any(lesser) || !any(greater)) return(NA)
        same.minus1 <- max(ref.l[lesser,]$nu)
        same.plus1 <- min(ref.l[greater,]$nu)
        
        other.l <- if (l0 == 0 && l1 == 1) ell.1 else ell.0
        lesser <- other.l$nu < ref.nu
        greater <- other.l$nu > ref.nu
        if (!any(lesser) || !any(greater)) return(NA)
        other.minus1 <- max(other.l[lesser,]$nu)
        other.plus1 <- min(other.l[greater,]$nu)
        
        const <- if (l0 == 0 && l1 == 1) 1/8 else -1/8
        const * ( same.minus1 - 4 * other.minus1 + 6 * ref.nu - 
                                4 * other.plus1 + same.plus1 )
    }
}

## Separations and ratios
dnu <- function(l, n, DF) separation(l, n, l+2, n-1, DF, use_n=F)
Dnu <- function(l, n, DF) separation(l, n, l, n-1, DF, use_n=F)
r_sep <- function(l, n, DF) dnu(l, n, DF) / Dnu(1-l, n+l, DF)
r_avg <- function(l, n, DF) dd(l, 1-l, n, DF, use_n=F) / Dnu(1-l, n+l, DF)

## Plot Dnu, dnu, r02, ... 
seismology_plot <- function(seps, nus, #fit, 
        gaussian_env, w.median, 
        nu_max, l_degs, ylab, dnu.cl, pchs, freqs, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    fit <- lm(seps~nus, weights=gaussian_env)
    if (length(l_degs)==1)
        col.pal <- colorRampPalette(c(dnu.cl[1], dnu.cl[3]))(1001)[1+1000*
            normalize(gaussian_env)]
    plot(seps~nus, axes=FALSE, tck=0, #xaxs='i',
         cex=1.5 * gaussian_env/max(gaussian_env), 
         ylab=ylab, 
         xlab=expression("Frequency" ~ nu / mu*Hz), 
         #xlim=c(1000, max(freqs$nu)), 
         ylim=range(w.median, 
                    coef(fit)[1], 
                    2*w.median-coef(fit)[1]), 
                    #seps), 
         col=if (length(l_degs)==1) col.pal else dnu.cl[pchs], 
         pch=if (length(l_degs)==1) 1 else pchs)
    abline(fit, lty=2)
    abline(v=nu_max, lty=3)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    if (length(l_degs)>1)
        legend("topleft", pch=l_degs+1, col=dnu.cl, 
               #cex=text.cex, #bty="n",
               cex=0.8*text.cex,
               legend=paste0("\u2113=", l_degs))#, horiz=1)
}

echelle_plot <- function(freqs, large_sep=NA, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    if (is.na(large_sep)) large_sep <- avg(Dnu, NULL, freqs, 0)[[1]] 
    plot(freqs$nu %% large_sep, freqs$nu, 
         tck=0, axes=FALSE, 
         pch=freqs$l+1, 
         col=dnu.cl[freqs$l+1], 
         xlab=expression((nu ~ mod ~ Delta * nu) / mu * Hz), 
         ylab=expression(nu/mu*Hz))
    for (l_deg in 0:3) {
        lines(freqs[freqs$l==l_deg,]$nu %% large_sep, 
              freqs[freqs$l==l_deg,]$nu,
              lty=2)
    }
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    l_degs <- sort(unique(freqs$l))
    legend("top", pch=l_degs+1, col=dnu.cl[l_degs+1], cex=0.8*text.cex, 
           legend=c(paste0("\u2113=", l_degs)))
}

