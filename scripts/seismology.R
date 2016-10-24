#### Seismological calculations for stellar observations and models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path(dirname(sys.frame(1)$ofile), 'utils.R'))

invisible(library(matrixStats))
invisible(library(magicaxis))
invisible(library(RColorBrewer))
#invisible(library(mblm))

fwhm_conversion <- (2*sqrt(2*log(2)))

dnu.cl <- c("#ca0020", "#f4a582", "#0571b0", "#800080")

labs <- list(
    Dnu             = bquote(Delta*nu/mu*Hz), 
    Dnu0            = bquote(Delta*nu[0]/mu*Hz), 
    dnu02           = bquote(delta*nu[0*","*2]/mu*Hz), 
    r02             = bquote(r[0*","*2]), 
    r01             = bquote(r[0*","*1]), 
    dnu13           = bquote(delta*nu[1*","*3]/mu*Hz), 
    r13             = bquote(r[1*","*3]), 
    r10             = bquote(r[1*","*0])
)

parse_freqs <- function(fname, gyre=F) {
    if (gyre) {
        freqs <- read.table(fname, skip=5, header=1)
        data.frame(l=freqs$l, n=freqs$n_pg, n_p=freqs$n_p, n_g=freqs$n_g,
            nu=freqs$Re.freq., E=freqs$E_norm)
    } else {
        #read.table(fname, col.names=c('l', 'n', 'nu'))
        read.table(fname, col.names=c('l', 'n', 'nu', 'E'))
    }
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
    ells <- unique(freqs$l)
    
    if (!'n_g' %in% names(freqs))
        freqs <- cbind(freqs, data.frame(n_g=0))
    
    # remove mixed modes
    mixed.1 <- freqs[freqs$n_g != 0 & freqs$n_p != 0 & freqs$l == 1,]
    mixed.2 <- freqs[freqs$n_g != 0 & freqs$n_p != 0 & freqs$l == 2,]
    freqs <- freqs[freqs$n_g == 0,]
    #ell.0 <- freqs[freqs$l==0,]
    #e.0 <- splinefun(ell.0$nu, ell.0$E)
    #for (ell in ells) {
    #    l.ell <- freqs[freqs$l == ell,]
    #    mixed <- l.ell$E / e.0(l.ell$nu) > 10
    #    if (any(mixed)) {
    #        print(paste("Removing", length(mixed), "l =", ell, "mixed modes"))
    #        freqs <- freqs[-which(freqs$nu == l.ell[mixed,]$nu),]
    #    }
    #}
    
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
    
    ## interpolate echelle diagram
    ell.0 <- freqs[freqs$l == 0,]
    #central.nu <- ell.0$nu[find_closest(ell.0$nu, nu_max)$x]
    #echelle.center <- seis.DF$Dnu0 / 2
    #modded.nu <- central.nu %% seis.DF$Dnu0
    #central.offset <- (modded.nu - echelle.center)
    #centered.nus <- if (central.offset > 0) ell.0$nu - central.offset
    #else ell.0$nu + abs(central.offset)
    #echelle <- splinefun(centered.nus, centered.nus %% seis.DF$Dnu0)
    
    # find frequencies of mixed modes 
    for (mm.ell in 1:2) {
        mixed <- if (mm.ell == 1) mixed.1 else mixed.2
        if (nrow(mixed) > 0) {
            #print("nrow mixed > 0")
            
            # only accept one mode between adjacent l=0s
            nus <- c() # mixed$nu
            for (ii in 1:nrow(ell.0)) {
                # pick the one with the lowest inertia
                candidates <- mixed[mixed$nu >= ell.0$nu[ii] & 
                                    mixed$nu <= ell.0$nu[ii+1],]
                nus <- c(nus, candidates[which.min(candidates$E),]$nu)
            }
            
            #print(mixed$nu)
            dist_to_nu_max <- (nus - nu_max)**2
            mixed <- data.frame(rbind(nus[order(dist_to_nu_max)]))
            #print(mixed)
            names(mixed) <- paste0('nu_cross.', mm.ell, '.', 1:ncol(mixed))
            #print(mixed)
            seis.DF <- cbind(seis.DF, mixed)
            
            ## interpolate within echelle diagram and save distance
            #mixed.offset <- if (central.offset > 0) mixed - central.offset
            #else mixed + abs(central.offset)
            #echelle.dist <- mixed.offset %% seis.DF$Dnu0 - echelle(mixed.offset)
            #names(echelle.dist) <- paste0('echelle.', mm.ell, '.', 
            #    1:ncol(echelle.dist))
            #seis.DF <- cbind(seis.DF, echelle.dist)
        }
    }
    #print(seis.DF)
    return(seis.DF)
}

get_average <- function(sep_name, freqs, l_deg, nu_max=NA, outf=F, ...) {
    
    result <- NULL
    
    vals <- get_separations(sep_name, freqs, l_deg, nu_max)
    
    nus <- vals$nus
    seps <- vals$separations
    
    #if (sep_name == "dnu") {
    #    if (is.na(Dnu)) Dnu <- get_average("Dnu", freqs, 0)$Dnu0
    #    keep <- abs(seps) < Dnu/2
    #    nus <- nus[keep]
    #    seps <- seps[keep]
    #}
    
    if (is.null(nus) | is.null(freqs)) return(NULL)
    not.na <- complete.cases(nus) & complete.cases(seps)
    nus <- nus[not.na]
    seps <- seps[not.na]
    
    if (length(seps)<5) {
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
        ratios <- dnorm(nus, nu_max, fwhm)/dnorm(nu_max, nu_max, fwhm)
        if (any(ratios < 0.05))
            w.median <- weightedMedian(seps, gaussian_env)
        else {
            print(paste("All points too far from nu_max for", sep_name))
            return(NULL)
        }
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
        ylab=as.expression(labs[sep_name]), 
        dnu.cl=dnu.cl, sep_name=sep_name, 
        freqs=freqs,
        ...)
    
    result
}

get_separations <- function(sep_name, freqs, l_deg, nu_max=NA) {
    if (sep_name == 'Dnu') {
        # Dnu(n, l) = nu_[n, l] - nu_[n-1, l]
        nus <- sort(freqs[freqs$l==l_deg,]$nu)
        Dnus <- diff(nus)
        list(nus=nus[-1], separations=Dnus)
    
    } else if (sep_name == 'dnu') {
        # dnu(n, l) = nu_[n, l] - nu_[n-1, l+2]
        x <- freqs[freqs$l==l_deg,]$nu
        y <- freqs[freqs$l==l_deg+2,]$nu
        indices <- find_closest(x, y)
        nus <- x[indices$x]
        separations <- x[indices$x]-y[indices$y]
        Dnu <- get_average("Dnu", freqs, 0, nu_max)$Dnu0
        rem <- abs(separations) > Dnu / 2
        separations[rem] <- NA
        list(nus=nus, separations=separations)
    
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
seismology_plot <- function(seps, nus, 
        gaussian_env, w.median, 
        nu_max, ylab, dnu.cl, pchs, freqs, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    fit <- lm(seps~nus, weights=gaussian_env)
    col.pal <- colorRampPalette(c(dnu.cl[1], dnu.cl[3]))(1001)[1+1000*
        normalize(gaussian_env)]
    plot(seps~nus, axes=FALSE, tck=0, #xaxs='i',
         cex=1.5 * gaussian_env/max(gaussian_env), 
         ylab=ylab, 
         xlab="", 
         xlim=range(freqs$nu), 
         ylim=range(w.median, 
                    coef(fit)[1], 
                    2*w.median-coef(fit)[1]), 
                    #seps), 
         col=col.pal, 
         pch=1)
    #abline(fit, lty=2)
    abline(v=nu_max, lty=3)
    magaxis(side=2:4, family=font, tcl=0.25, labels=c(1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    new.mgp <- sapply(mgp-c(0.2, 0.2, 0), function(x) max(0,x))
    par(mgp=new.mgp)
    magaxis(side=1, family=font, tcl=0.25, labels=1, mgp=new.mgp, 
        las=1, cex.axis=text.cex)
    title(xlab=expression("Frequency" ~ nu / mu*Hz))
}

echelle_plot <- function(freqs, large_sep=NA, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    if (is.na(large_sep)) large_sep <- avg(Dnu, freqs)[[1]] 
    nus <- freqs$nu
    nus <- c(nus %% large_sep, (nus %% large_sep) + large_sep)
    fwhm <- (0.66*nu_max**0.88)/fwhm_conversion
    gaussian_env <- dnorm(freqs$nu, nu_max, fwhm)
    cexs <- rep(2 * gaussian_env / max(gaussian_env), 2)
    plot(nus, rep(freqs$nu, 2), 
         tck=0, axes=FALSE, 
         pch=freqs$l+1, 
         cex=cexs, 
         col=dnu.cl[freqs$l+1], 
         xlab="", 
         ylab=expression(nu/mu*Hz),
         xlim=c(0, large_sep*2))
    magaxis(side=2:4, family=font, tcl=0.25, labels=c(1,0,0), mgp=mgp, 
        las=1, cex.axis=text.cex)
    magaxis(side=1, family=font, tcl=0.25, labels=F, mgp=mgp-c(0, 0.2, 0), 
        las=1, cex.axis=text.cex)
    par(mgp=mgp-c(0.2, 0, 0))
    title(xlab=expression((nu ~ mod ~ Delta * nu[0]) / mu * Hz))
    ticks <- axTicks(1)
    ticks <- ticks[ticks < large_sep]
    ticks <- c(ticks, ticks+large_sep)
    new.mgp <- sapply(mgp-c(0, 0.2, 0), function(x) max(0,x))
    axis(1, tcl=0.25, mgp=new.mgp, cex.axis=text.cex, tick=F, 
        at=ticks, labels=ticks%%large_sep)
    l_degs <- sort(unique(freqs$l))
    abline(v=large_sep, lty=3)
    legend("left", pch=l_degs+1, col=dnu.cl[l_degs+1], cex=0.8*text.cex, 
           inset=c(0.02, 0), legend=c(paste0("l=", l_degs)))
}

plot_power_spectrum <- function(freqs, nu_max, colors=0, ..., 
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    nus <- freqs$nu
    fwhm <- (0.66*nu_max**0.88)/fwhm_conversion
    power <- dnorm(nus, nu_max, fwhm)
    plot(nus, power/max(power),
         #ylim=c(0, 1.05),
         ylab="",#"Power", 
         xlim=range(nus[power/max(power)>0.001]),
         xlab=expression("Frequency"~nu/mu*Hz),
         type='h', lwd=1.5, 
         xaxs='i', yaxs='i', axes=FALSE, 
         col=adjustcolor('black', alpha=0.5))
    #points(nus, power/max(power), type='h', lwd=2,
    #    col=adjustcolor(col.pal[freqs$l+1], alpha=0.25))
    new.mgp <- sapply(mgp-c(0.2, 0.2, 0), function(x) max(0,x))
    magaxis(1, labels=1, tcl=-0.25, cex.axis=text.cex, mgp=new.mgp, 
        family=font)
    legend("topleft", bty="n", col=col.pal, lty=1, cex=text.cex, legend=c(
        expression("l"==0),
        expression("l"==1),
        expression("l"==2),
        expression("l"==3)
    ))
    
    ## annotate Dnu0
    ell.0 <- freqs$l==0
    ell.pwrs <- power[ell.0]
    max.pwr <- max(ell.pwrs)
    max.pwr.idx <- which(power[ell.0] == max.pwr)
    right <- freqs[ell.0,]$nu[max.pwr.idx-3]
    left <- freqs[ell.0,]$nu[max.pwr.idx-4]
    right.pwr <- ell.pwrs[max.pwr.idx-3]/max(power)
    points(right, right.pwr, 
        type='h', col=col.pal[1], lwd=2)
    left.pwr <- ell.pwrs[max.pwr.idx-4]/max(power)
    points(left, left.pwr, 
        type='h', col=col.pal[1], lwd=2)
    arrows(left+1, left.pwr-0.15, right-1, length=0.1)
    arrows(right-1, left.pwr-0.15, left+1, length=0.1)
    text((left+right)/2 - 35, left.pwr+0.19, cex=text.cex, 
        expression(Delta*nu[0]))
    
    ## annotate dd01
    ell.1 <- freqs$l==1
    ell.0.nus <- freqs[ell.0,]$nu
    ell.1.nus <- freqs[ell.1,]$nu
    center <- freqs[ell.0,]$nu[max.pwr.idx]
    
    without.center <- ell.0.nus[-which(ell.0.nus==center)]
    first.zero <- find_closest(center, without.center)$y
    first.zero.nu <- without.center[first.zero]
    
    second.zero <- find_closest(center, without.center[-first.zero])$y
    second.zero.nu <- without.center[-first.zero][second.zero]
    
    first.one <- find_closest(center, ell.1.nus)$y
    first.one.nu <- ell.1.nus[first.one]
    
    second.one <- find_closest(center, ell.1.nus[-first.one])$y
    second.one.nu <- ell.1.nus[-first.one][second.one]
    
    center.pwr <- power[which(freqs$nu==center)]/max(power)
    points(center, center.pwr, 
        type='h', col=col.pal[1], lwd=2)
    first.zero.nu.pwr <- power[which(freqs$nu==first.zero.nu)]/max(power)
    points(first.zero.nu, first.zero.nu.pwr, 
        type='h', col=col.pal[1], lwd=2)
    second.zero.nu.pwr <- power[which(freqs$nu==second.zero.nu)]/max(power)
    points(second.zero.nu, second.zero.nu.pwr, 
        type='h', col=col.pal[1], lwd=2)
    first.one.nu.pwr <- power[which(freqs$nu==first.one.nu)]/max(power)
    points(first.one.nu, first.one.nu.pwr, 
        type='h', col=col.pal[2], lwd=2)
    second.one.nu.pwr <- power[which(freqs$nu==second.one.nu)]/max(power)
    points(second.one.nu, second.one.nu.pwr, 
        type='h', col=col.pal[2], lwd=2)
    others <- c(first.zero.nu, second.zero.nu, first.one.nu, second.one.nu)
    for (other in 1:4) {
        arrows(center, left.pwr-0.15, others[other], 
        ifelse(other == 1 || other == 2, -0.15, -0.05) + left.pwr, length=0.1)
    }
    text(0.97*second.one.nu, 0.95, expression(dd[0*","*1]), cex=text.cex, 
        adj=c(1, NA))
    
    ## annotate dnu02
    ell.2.nus <- freqs[freqs$l==2,]$nu
    ell.pwrs <- power[ell.0]
    max.pwr <- max(ell.pwrs)
    max.pwr.idx <- which(power[ell.0] == max.pwr)
    right <- freqs[ell.0,]$nu[max.pwr.idx+3]
    left <- ell.2.nus[find_closest(right, ell.2.nus)$y]
    points(right, power[which(freqs$nu == right)]/max(power), 
        type='h', col=col.pal[1], lwd=2)
    points(left, power[which(freqs$nu == left)]/max(power), 
        type='h', col=col.pal[3], lwd=2)
    arrows(left-1, left.pwr-0.15, left, length=0.1)
    arrows(right+1, left.pwr-0.15, right, length=0.1)
    arrows(left, left.pwr-0.15, right, length=0)
    text(left*0.99, left.pwr+0.26, cex=text.cex, 
        expression(delta*nu[0*","*2]), adj=c(0, NA))
    
    ## annotate dnu13
    ell.3.nus <- freqs[freqs$l==3,]$nu
    ell.pwrs <- power[ell.1]
    max.pwr <- max(ell.pwrs)
    max.pwr.idx <- which(power[ell.1] == max.pwr)
    right <- freqs[ell.1,]$nu[max.pwr.idx+5]
    left <- ell.3.nus[find_closest(right, ell.3.nus)$y]
    points(right, power[which(freqs$nu == right)]/max(power), 
        type='h', col=col.pal[2], lwd=2)
    points(left, power[which(freqs$nu == left)]/max(power), 
        type='h', col=col.pal[4], lwd=2)
    arrows(left-1, left.pwr-0.15, left, length=0.1)
    arrows(right+1, left.pwr-0.15, right, length=0.1)
    arrows(left, left.pwr-0.15, right, length=0)
    text(left*0.99, left.pwr, cex=text.cex, 
        expression(delta*nu[1*","*3]), adj=c(0, NA))
}

## nu_max scaling relation
nu_max_scaling <- function(M, R, Teff, Teff_sun=5777, nu_max_sun=3090) {
    M * R**-2 * (Teff/Teff_sun)**(-1/2) * nu_max_sun
}

