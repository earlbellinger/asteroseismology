#### Seismological calculations for stellar observations and models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path(dirname(sys.frame(1)$ofile), 'utils.R'))

invisible(library(matrixStats))
invisible(library(magicaxis))
invisible(library(RColorBrewer))
invisible(library(mblm))

fwhm_conversion <- (2*sqrt(2*log(2)))

dnu.cl <- c("#ca0020", "#f4a582", "#0571b0", "#800080")

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
    ells <- unique(freqs$l)
    
    # get l=0
    seis.DF <- get_average("Dnu", freqs, 0, nu_max, outf, ...)
    
    # make echelle
    if (outf != FALSE) make_plots(echelle_plot, paste0(outf, '-echelle'), 
        freqs=freqs, large_sep=seis.DF[[1]], ...)
    
    # get averages for dnu, r01, r02, r10, r13
    for (l_deg in 0:1) {
        if (l_deg %in% ells && (l_deg+2) %in% ells) { 
            dnu_median <- get_average("dnu", freqs, l_deg, nu_max, outf, ...)
            if (!is.null(dnu_median)) seis.DF <- cbind(seis.DF, dnu_median)
            
            if ((1-l_deg) %in% ells) {
                rsep_median <- get_average("r_sep", freqs, l_deg, nu_max, 
                    outf, ...)
                if (!is.null(rsep_median)) 
                    seis.DF <- cbind(seis.DF, rsep_median)
            }
        }
        if (0 %in% ells && 1 %in% ells) {
            ravg_median <- get_average("r_avg", freqs, l_deg, nu_max, outf, ...)
            if (!is.null(ravg_median)) seis.DF <- cbind(seis.DF, ravg_median)
        }
    }
    return(seis.DF)
}

get_average <- function(sep_name, freqs, l_deg, nu_max=NA, outf=F, ...) {
    result <- NULL
    
    vals <- get_separations(sep_name, freqs, l_deg)
    
    nus <- vals$nus
    seps <- vals$separations
    
    not.na <- complete.cases(nus) & complete.cases(seps)
    nus <- nus[not.na]
    seps <- seps[not.na]
    
    if (length(seps)<2) {
        print(paste("Too few points for", sep_name))
        return(NULL)
    }
    
    sep_name <- if (sep_name == 'Dnu')   paste0(sep_name, l_deg)
           else if (sep_name == 'dnu')   paste0(sep_name, l_deg, l_deg+2)
           else if (sep_name == 'r_sep') paste0('r', l_deg, l_deg+2)
           else if (sep_name == 'r_avg') paste0('r', l_deg, 1-l_deg)
    
    if (!is.na(nu_max)) {
        fwhm <- (0.66*nu_max**0.88)/fwhm_conversion
        gaussian_env <- dnorm(nus, nu_max, fwhm)
        w.median <- weightedMedian(seps, gaussian_env)
        #int.weights <- floor(gaussian_env*10000)
        #new.nus <- rep(nus, int.weights)
        #new.seps <- rep(seps, int.weights)
        #fit <- mblm(new.seps~new.nus)
        #fit <- lm(seps~nus, weights=gaussian_env)
        result <- data.frame(w.median)#, coef(fit)[[2]])
    } else {
        #fit <- lm(seps~nus)
        result <- data.frame(median(seps))#, coef(fit)[[2]])
    }
    colnames(result) <- sep_name #c(sep_name, paste0(sep_name, "_slope"))
    
    if (outf != FALSE && !is.na(result)) make_plots(seismology_plot, 
        paste0(outf, '-', sep_name), 
        seps=seps, nus=nus, 
        gaussian_env=gaussian_env, 
        w.median=result[sep_name], nu_max=nu_max, 
        ylab=as.expression(get_label_nameless(sep_name)), 
        dnu.cl=dnu.cl, sep_name=sep_name, 
        freqs=freqs,
        ...)
    
    result
}

get_separations <- function(sep_name, freqs, l_deg) {
    if (sep_name == 'Dnu') {
        # Dnu(n, l) = nu_[n, l] - nu_[n-1, l]
        nus <- sort(freqs[freqs$l==l_deg,]$nu)
        list(nus=nus[-1], separations=diff(nus))
    
    } else if (sep_name == 'dnu') {
        # dnu(n, l) = nu_[n, l] - nu_[n-1, l+2]
        x <- freqs[freqs$l==l_deg,]$nu
        y <- freqs[freqs$l==l_deg+2,]$nu
        indices <- find_closest(x, y)
        list(nus=x[indices$x], separations=x[indices$x]-y[indices$y])
    
    } else if (sep_name == 'r_sep') { 
        # r(l, n) = dnu(n, l) / Dnu(n+l, 1-l)
        # dnu(n, l) = nu_[n, l] - nu_[n-1, l+2]
        # Dnu(n, l) = nu_[n, l] - nu_[n-1, l]
        dnus <- get_separations("dnu", freqs, l_deg)
        Dnus <- get_separations("Dnu", freqs, 1-l_deg)
        x <- dnus$nu
        y <- Dnus$nu
        indices <- find_closest(x, y)
        dnus <- dnus$separations[indices$x]
        Dnus <- Dnus$separations[indices$y]
        list(nus=x[indices$x], separations=dnus/Dnus)
    
    } else if (sep_name == 'r_avg') {
        # r_[l, 1-l](n) = dd_[l, 1-l](n) / Dnu(n+l, 1-l)
        # dd_01(n) =  1/8 ( nu_[n-1,0] - 4*nu_[n-1,1] + 
        #                 6*nu_[n,  0] - 4*nu_[n,  1] + nu_[n+1, 0] )
        # dd_10(n) = -1/8 ( nu_[n-1,1] - 4*nu_[n,  0] + 
        #                 6*nu_[n,  1] - 4*nu_[n+1,0] + nu_[n+1, 1] )
        nus.0 <- sort(freqs[freqs$l==0,]$nu)
        nus.1 <- sort(freqs[freqs$l==1,]$nu)
        x <- if (l_deg == 0) nus.0 else nus.1 # x are freqs of the same l 
        y <- if (l_deg == 0) nus.1 else nus.0 # y are freqs of the other l 
        Dnus.x <- if (l_deg == 0) diff(nus.0) else diff(nus.1)
        Dnus.factor <- if (l_deg == 0) 8 else -8
        large_sep <- median(Dnus.x, na.rm=1)
        
        separations <- c()
        nus <- c()
        Dnus <- c()
        for (n in 1:length(x)) {
            # nu_[n-1] of the same l
            xs <- x[x<x[n]]
            x.smaller <- find_closest(x[n], xs)
            if (is.null(x.smaller$x)) next
            x.smaller <- xs[x.smaller$y]
            if (x.smaller < (x[n]-1.5*large_sep)) next
            
            # nu_[n-1] or nu_[n] of the other l, 
            # depending on whether l=0 or l=1
            ys <- y[y<x[n] & y>x.smaller]
            y.smaller <- find_closest(x[n], ys)
            if (is.null(y.smaller$x)) next
            y.smaller <- ys[y.smaller$y]
            
            # nu_[n] or nu_[n+1] of the other l
            eligible <- y>x[n]
            ys <- y[eligible]
            y.bigger <- find_closest(x[n], ys)
            if (is.null(y.bigger$x)) next
            y.bigger <- ys[y.bigger$y]
            
            # nu_[n+1] of the same l 
            xs <- x[x>x[n] & x>y.bigger]
            x.bigger <- find_closest(x[n], xs)
            if (is.null(x.bigger$x)) next
            x.bigger <- xs[x.bigger$y]
            if (x.bigger > (x[n]+1.5*large_sep)) next
            
            separations <- c(separations,
                x.smaller - 4*y.smaller + 6*x[n] - 4*y.bigger + x.bigger)
            nus <- c(nus, x[n])
            Dnus <- c(Dnus, Dnus.factor * (y.bigger - y.smaller))
            
        }
        list(nus=nus, separations=separations/Dnus)
    }
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
        d.freqs <<- freqs # destructively gets modified by separation
        vals <- sapply(unique(ell$n), function(n) separator(l_deg, n, d.freqs))
        #print(d.freqs)
        not.nan <- complete.cases(vals)
        is.paired <- vals < 2*median(vals, na.rm=T)
        seps <- c(seps, vals[not.nan & is.paired])
        nus <- c(nus, ell$nu[not.nan & is.paired])
        if (outf != FALSE) pchs <- c(pchs, rep(l_deg+1, sum(not.nan)))
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
separation <- function(first_l, first_n, second_l, second_n, DF, use_n=F) {
    first <- DF[DF$l == first_l & DF$n == first_n,]
    if (nrow(first) != 1) return(NA)
    
    if (use_n) { 
        second <- DF[DF$l == second_l & DF$n == second_n,]
    } else {
        n.diff <- first_n - second_n
        second.modes <- DF[DF$l == second_l & DF$nu < first$nu,]
        if (all(is.na(second.modes))) return (NA)
        second <- second.modes[order(second.modes$nu, decreasing=T),][n.diff,]
        if (first_l != second_l) {
            same.ell <- DF[DF$l == first_l & 
                           DF$nu < first$nu & 
                           DF$nu > second$nu,]
            if (nrow(same.ell) > 0) return(NA)
        }
    }
    if (nrow(second) != 1) return(NA)
    difference <- first$nu - second$nu
    if (difference <= 0) return(NA)
    d.freqs <<- DF[-which(DF$nu == second$nu),]
    return(difference)
}

## Five point averages 
#dd_01= 1/8( nu_[n-1,0] - 4*nu_[n-1,1] + 6*nu_[n,0] - 4*nu[n,  1] + nu_[n+1,0] )
#dd_10=-1/8( nu_[n-1,1] - 4*nu_[n,  0] + 6*nu_[n,1] - 4*nu[n+1,0] + nu_[n+1,1] )
# with the subscript of nu being (n, l)
dd <- function(l0, l1, n, DF, use_n=F) {
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
        nu_max, ylab, dnu.cl, pchs, freqs, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    fit <- lm(seps~nus, weights=gaussian_env)
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
         col=col.pal, 
         pch=1)
    abline(fit, lty=2)
    abline(v=nu_max, lty=3)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
}

echelle_plot <- function(freqs, large_sep=NA, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    if (is.na(large_sep)) large_sep <- avg(Dnu, freqs)[[1]] 
    nus <- freqs$nu
    nus <- c(nus %% large_sep, (nus %% large_sep) + large_sep)
    plot(nus, c(freqs$nu, freqs$nu), 
         tck=0, axes=FALSE, 
         pch=freqs$l+1, 
         col=dnu.cl[freqs$l+1], 
         xlab=expression((nu ~ mod ~ Delta * nu) / mu * Hz), 
         ylab=expression(nu/mu*Hz),
         xlim=c(0, large_sep*2))
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(0,1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    ticks <- axTicks(1)
    ticks <- ticks[ticks < large_sep]
    ticks <- c(ticks, ticks+large_sep)
    axis(1, tcl=0.25, mgp=mgp, cex.axis=text.cex, tick=F, 
        at=ticks, labels=ticks%%large_sep)
    l_degs <- sort(unique(freqs$l))
    abline(v=large_sep, lty=3)
    legend("bottomleft", pch=l_degs+1, col=dnu.cl[l_degs+1], cex=0.8*text.cex, 
           legend=c(paste0("\u2113=", l_degs)))
}

