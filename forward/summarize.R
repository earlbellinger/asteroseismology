#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source('../scripts/seismology.R')
library(parallel)
library(parallelMap)

## Constants
Z_div_X_solar = 0.02293

profile.pattern <- 'profile.+.data$'
freqs.pattern <- 'profile.+-freqs.dat$'
freqs.cols <- c('l', 'n', 'nu', 'inertia')

### Obtain observable properties from models 
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
    
    freqs <- read.table(freqs_file, col.names=freqs.cols, fill=TRUE)
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
    profile_candidates <- logs[grep(profile.pattern, logs)]
    freq_file_candidates <- logs[grep(freqs.pattern, logs)]
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

plot_HR <- function(text.cex, ...) {
    col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(1000)
    plot(log10(DF$Teff), log10(DF$L), 
        type='l', tcl=0, axes=FALSE,
        xlab=expression('lg' ~ T[eff]/K ),
        ylab=expression('lg' ~ L / L['\u0298'] ),
        xlim=rev(log10(range(DF$Teff))))
    abline(v=log10(5771), lty=3, col='lightgray')
    abline(h=0, lty=3, col='lightgray')
    points(log10(5771), 0, pch=1, cex=1)
    points(log10(5771), 0, pch=20, cex=0.1)
    points(log10(DF$Teff), log10(DF$L), pch=1, cex=0.1, 
        col=col.pal[
            floor((DF$Hc-min(DF$Hc)) /
                  (max(DF$Hc)-min(DF$Hc))*length(col.pal)) + 1])
    magaxis(side=1:4, family=font, tcl=0.25, mgp=utils.mgp, labels=c(1,1,0,0)) 
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        signif(quantile(seq(min(DF$Hc), max(DF$Hc), length=1000), 
            c(0, 0.25, 0.5, 0.75, 1)), 3), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(expression(H[c]), 4, line=4.5, cex=text.cex)
    dev.off()
}

### Parse command line arguments and process directory
args <- commandArgs(TRUE)
if (length(args)>0) {
    print(args[1])
    DF <- unique(parse_dir(args[1]))
    DF <- DF[complete.cases(DF),]
    DF <- DF[order(DF$age),]
    
    ## remove PMS
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
    make_plots(plot_HR, mar=c(3, 4, 1, 6))
    
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
    } else { # Make a table of results! 
        write.table(DF, paste0(args[1], '.dat'), 
                    quote=FALSE, sep='\t', row.names=FALSE)
    }
}
