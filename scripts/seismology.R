#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

options(error=traceback)

library(matrixStats)
library(magicaxis)
library(RColorBrewer)
library(parallel)
library(parallelMap)

cl <- brewer.pal(4, "BrBG")
#rbPal <- colorRampPalette(c('red','blue'))

Z_div_X_solar = 0.02293

freqs_cols <- c('l', 'n', 'nu', 'inertia')
profile_pattern <- 'profile.+.data$'
freqs_pattern <- 'profile.+-freqs.dat$'

plot_width <- 6.97522
plot_height <- 4.17309
font <- "Palatino"

separation_dir <- 'plots'
dir.create(separation_dir, showWarnings=FALSE)

################################################################################
### Seismological calculations #################################################
################################################################################

## Separation: just the difference between two frequencies 
separation <- function(first_l, first_n, second_l, second_n, df) {
    # nu_{l1,n1} - nu_{l2,n2}
    first <- df$l == first_l & df$n == first_n
    second <- df$l == second_l & df$n == second_n
    if (sum(first) == 1 && sum(second) == 1) { # check that it's unique
        difference <- df[first,]$nu - df[second,]$nu
        if (difference < 0) return(NA)
        return(difference)
    }
    return(NA)
}

# Five point averages defined by
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

get_averages <- function(f, df, freqs, l_degs, nu_max, outf=FALSE) {
    # calcualte averages of things like f = dnu, Dnu, r_sep, r_avg
    # df is the where the result will be stored
    # freqs are a data frame with columns l, n, nu
    # l_degs are the l's for which this calculation should be made 
    # nu_max is the center of the gaussian
    # make a plot with filename 'outf' if outf != FALSE
    sep_name <- deparse(substitute(f))
    a <- c() # contains the computed quantity (e.g. large freq separations)
    b <- c() # contains frequencies of the base mode
    pchs <- c() # if there's more than one l, get different symbols for each
    #err <- c() # uncertainties on frequencies if they are known (not models)
    for (l_deg in l_degs) {
        ell <- freqs[freqs$n > 1 & freqs$l==l_deg,]
        vals <- sapply(unique(ell$n), function(n) f(l_deg, n, freqs))
        not.nan <- complete.cases(vals)
        a <- c(a, vals[not.nan])
        b <- c(b, ell$nu[not.nan])
        pchs = c(pchs, rep(l_deg+1, sum(not.nan)))
        #if ("dnu" %in% names(freqs)) err <- c(err, 1/ell$dnu[not.nan])
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
    #if ("dnu" %in% names(freqs)) err 
    wm <- weightedMedian(a, gaussian_env)
    df[paste0(sep_name, "_median")] <- wm
    fit <- lm(a~b, weights=gaussian_env)
    df[paste0(sep_name, "_slope")] <- coef(fit)[2]
    #df[paste0(sep_name, "_intercept")] <- coef(fit)[1]
    #lower.bound = wm-coef(fit)[1]
    #upper.bound = wm+coef(fit)[1]
    
    if (outf != FALSE) {
        cairo_pdf(file.path(separation_dir, 
                            paste0(sep_name, '-', outf, '.pdf')),
                  width=plot_width, height=plot_height, family=font)
        par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
        plot(a~b, tck=0, ylab=as.expression(ylab), 
             cex=gaussian_env*1.75/max(gaussian_env),
            #ylim=c(lower.bound-lower.bound*0.05, upper.bound+upper.bound*0.05),
             ylim=range(wm, coef(fit)[1], wm+(wm-coef(fit)[1])),
             col=if (length(l_degs)==1) 1 else cl[pchs], 
             pch=if (length(l_degs)==1) 1 else pchs, 
             xlab=expression("Frequency" ~ nu / mu*Hz))
        abline(fit, lty=2)
        abline(v=nu_max, lty=3)
        magaxis(side=1:4, family=font, tcl=0.25, labels=FALSE)
        if (length(l_degs)>1)
            legend("topright", pch=l_degs+1, col=cl, cex=0.75,
                   ncol=length(l_degs), #bty="n",
                   legend=paste0("\u2113=", l_degs))
        dev.off()
    }
    
    df
}

seismology <- function(freqs, nu_max, acoustic_cutoff=Inf, outf=FALSE) {
    if (nrow(freqs) == 0) {
        print("No frequencies found")
        return(NULL)
    }
    freqs <- unique(freqs[complete.cases(freqs) & freqs$nu < acoustic_cutoff,])
    
    # fix radial modes because ADIPLS breaks sometimes
    for (l_mode in unique(freqs$l)) {
        # grab the relevant l's and n's
        selection <- freqs$l==l_mode & freqs$n>0
        ell <- freqs[selection,]
        ns <- ell$n
        # check if any n's are duplicated and if so, shift them down
        if (any(duplicated(ns))) {
            return(NULL) # screw it, just discard this data
            dup <- which(duplicated(ns)) # grab duplicated (hopef. only one)
            if (length(dup) > 1) { # hopeless
                print(paste0("Duplicated l=", l_mode, " mode, exiting"))
                return(NULL) 
            }
            toshift <- 1:(dup-1) # find the ones to shift
            ell$n[toshift] <- ns[toshift] - 1 # calculate new n vals 
            freqs[selection,] <- ell # replace the old with the new
            freqs <- freqs[!(freqs$l==l_mode & freqs$n==0),] # overwrite data
        }
    }
    
    seis.DF <- NULL
    seis.DF <- get_averages(Dnu, seis.DF, freqs, sort(unique(freqs$l)), 
        nu_max, outf)
    seis.DF <- get_averages(Dnu, seis.DF, freqs, 0, nu_max, outf)
    for (l_deg in 0:1) {
        if (l_deg==0 || 3 %in% freqs$l) {
            seis.DF <- get_averages(dnu, seis.DF, freqs, l_deg, nu_max, outf)
            seis.DF <- get_averages(r_sep, seis.DF, freqs, l_deg, nu_max, outf)
        }
        seis.DF <- get_averages(r_avg, seis.DF, freqs, l_deg, nu_max, outf)
    }
    return(seis.DF)
}

################################################################################
### Obtain observable properties from models ###################################
################################################################################

get_obs <- function(profile_file, freqs_file, ev_history) {
    print(profile_file)
    
    profile_header <- read.table(profile_file, header=TRUE, nrows=1, skip=1)
    hstry <- ev_history[ev_history$model_number==profile_header$model_number,]
    if (nrow(hstry) == 0) {#|| hstry$mass_conv_core > 0) 
        print(paste("Model", profile_file, "failed"))
        return(NULL)
    }
    
    obs.DF <- NULL
    ## Things we want to predict
    obs.DF["age"] <- profile_header$star_age/10**9
    #obs.DF["mass"] <- profile_header$star_mass
    obs.DF["radius"] <- profile_header$photosphere_r
    obs.DF["H"] <- profile_header$star_mass_h1/profile_header$star_mass
    obs.DF["Hc"] <- hstry$center_h1
    obs.DF["He"] <- (profile_header$star_mass_he3 + 
            profile_header$star_mass_he4)/profile_header$star_mass
    obs.DF["mass_cc"] <- hstry$mass_conv_core/profile_header$star_mass
    
    ## Things we can observe
    obs.DF["log_g"] <- hstry$log_g
    obs.DF["L"] <- profile_header$photosphere_L
    obs.DF["Teff"] <- profile_header$Teff
    obs.DF["Fe/H"] <- log10(10**hstry$log_surf_z/hstry$surface_h1/Z_div_X_solar)
    
    #if (hstry$mass_conv_core > 0) {
    #    print(paste("ConvectiveCore", profile_file, 
    #                obs.DF["age"], obs.DF["mass"], obs.DF["He"],
    #                hstry$mass_conv_core, 
    #                hstry$mass_conv_core/profile_header$star_mass))
    #}
    
    freqs <- read.table(freqs_file, col.names=freqs_cols, fill=TRUE)
    acoustic_cutoff <- hstry$acoustic_cutoff/(2*pi)
    nu_max <- hstry$nu_max
    seis.DF <- seismology(freqs, nu_max, acoustic_cutoff, 
        outf=ifelse(sample(0:10000, 1)==0, gsub("/", "-", freqs_file), FALSE))
        
    return(merge(rbind(obs.DF), rbind(seis.DF)))
}

### Obtain evolutionary tracks from a MESA directory
parse_dir <- function(directory) {
    #print(directory)
    
    # parse dirname string e.g. "M=1.0_Y=0.28"
    params.DF <- NULL
    for (var in unlist(strsplit(basename(directory), '_'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params.DF[nameval[1]] <- as.numeric(nameval[2])
    }
    
    # obtain history
    log_dir <- file.path(directory, "LOGS")
    logs <- list.files(log_dir)
    if (length(logs) <= 1) {
        print(paste(directory, "No logs found!"))
        return(NA)
    }
    ev_history <- read.table(file.path(log_dir, 'history.data'), 
            header=TRUE, skip=5)
    
    # figure out which profiles & frequency files to use
    profile_candidates <- logs[grep(profile_pattern, logs)]
    freq_file_candidates <- logs[grep(freqs_pattern, logs)]
    profile_files <- c()
    freq_files <- c()
    for (profile_file in profile_candidates) {
        freq_name <- sub(".data", "-freqs.dat", profile_file, fixed=TRUE)
        if (freq_name %in% freq_file_candidates) {
            profile_files <- c(profile_files, profile_file)
            freq_files <- c(freq_files, freq_name)
        }
    }
    if (length(profile_files) <= 2) {
        print("Too few profile files")
        return(NA)
    }
    
    # obtain observable information
    parallelStartMulticore(max(1, detectCores()))
    obs.DF <- do.call(plyr:::rbind.fill, 
        parallelMap(function(profile_file, freqs_file)
                get_obs(profile_file, freqs_file, ev_history), 
            profile_file=file.path(log_dir, profile_files), 
            freqs_file=file.path(log_dir, freq_files)))
    
    return(merge(rbind(params.DF), obs.DF[with(obs.DF, order(age)),]))
}

### Parse command line arguments and process directory
args <- commandArgs(TRUE)
if (length(args)>0) {
    print(args[1])
    DF <- unique(parse_dir(args[1]))
    DF <- DF[complete.cases(DF),]
    DF <- DF[order(DF$age),]
    
    # remove PMS
    # find all models within 2% of starting H
    # if the temperature and luminosity changes have outliers,
    # then remove all up to the last outlier 
    hs <- 1-DF$Y[1]-DF$Z[1]
    dteff <- diff(DF$Teff)
    dl <- diff(DF$L)
    pms <- with(DF, 
        which(100*(hs-H[-1])/(hs) < 2
            & dteff %in% boxplot.stats(dteff)$out
            & dl %in% boxplot.stats(dl)$out
        )
    )
    if (any(pms)) DF <- DF[-1:-(1+max(pms)),]
    
    # set ZAMS age and eliminate models older than the universe
    minage <- min(DF$age)
    DF$age <- DF$age - minage
    DF <- DF[DF$age <= 13.82,]
    
    # set ZAMS composition
    DF$Y <- DF$He[1]
    DF$Z <- 1 - DF$H[1] - DF$He[1]
    
    # make plot
    png(file.path(separation_dir, paste0('HR-', args[1], '.png')),
                  res=400, width=250*plot_width, height=250*plot_height, 
                  family=font)
    #cairo_pdf(file.path(separation_dir, paste0(args[1], '-HR.pdf')),
    #          width=plot_width, height=plot_height, family=font)
    col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(1000)
    par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
    plot(log10(DF$Teff), log10(DF$L), type='l', tcl=0, axes=FALSE,
        xlab=expression(log ~ T[eff]/K ),
        ylab=expression(log ~ L / L['\u0298'] ),
        xlim=rev(log10(range(DF$Teff))))
    abline(v=log10(5771), lty=3, col='lightgray')
    abline(h=0, lty=3, col='lightgray')
    points(log10(5771), 0, pch=1, cex=1)
    points(log10(5771), 0, pch=20, cex=0.1)
    points(log10(DF$Teff), log10(DF$L), pch=1, cex=0.1, 
        col=col.pal[
            floor((DF$Hc-min(DF$Hc)) /
                  (max(DF$Hc)-min(DF$Hc))*length(col.pal)) + 1])
    magaxis(side=1:4, family=font, tcl=0.25, 
            mgp=c(2, 0.25, 0), labels=c(1,1,0,0))
    var1range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        signif(quantile(seq(min(DF$Hc), max(DF$Hc), length=1000), 
            c(0, 0.25, 0.5, 0.75, 1)), 3), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(expression(H[c]), 4, line=4.5, cex=1.3)
    dev.off()
    
    # check if the track should be kept 
    if (nrow(DF) < 5) {
        print("Rejecting track: too few data points")
        print(DF)
    } else if (ncol(DF) < 5) {
        print("Rejecting track: too few columns")
        print(DF)
    } else if (any(DF$Teff>=7000)) {
        print("Rejecting track: too high temperatures")
        print(DF$Teff)
    } else if (any(DF[['Fe/H']] <= -5)) {
        print("Rejecting track: too low metallicities")
        print(fivenum(DF[['Fe/H']]))
    } else if ( 100*(hs-max(DF$H))/hs > 2 ) {
        print("Rejecting track: inaccurate starting H")
        print(c(hs, max(DF$Hc)))
    } else {
        write.table(DF, paste0(args[1], '.dat'), 
                    quote=FALSE, sep='\t', row.names=FALSE)
    }
}
