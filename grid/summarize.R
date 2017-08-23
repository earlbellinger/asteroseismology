#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'seismology.R'))
source(file.path('..', 'scripts', 'utils.R'))
library(parallel)
library(parallelMap)
library(lpSolve)
library(magicaxis)
library(mblm)

#options(warn=2) 

Z_div_X_solar = 0.02293
profile.pattern <- 'profile.+.data$'
freqs.pattern <- 'profile.+-freqs.dat$'
freqs.cols <- c('l', 'n', 'nu')#, 'inertia')
num_points <- 64
X_c.lim <- 1e-2
ev.stages <- list(
    'LOGS_MS'=1, 
    'LOGS_SG'=2, 
    'LOGS_RGB'=3, 
    'LOGS_BUMP'=4, 
    'LOGS_HEB'=5)

### Obtain observable properties from models 
summarize <- function(pro_file, freqs_file, ev.DF, dname) {
    print(pro_file)
    
    pro_header <- read.table(pro_file, header=TRUE, nrows=1, skip=1)
    pro_body <- read.table(pro_file, header=TRUE, nrows=1, skip=5)
    hstry <- ev.DF[ev.DF$model_number==pro_header$model_number,]
    if (nrow(hstry) == 0) {#|| hstry$mass_conv_core > 0) 
        print(paste("Model", pro_file, "failed: no history information"))
        return(NULL)
    }
    
    obs.DF <- NULL
    
    obs.DF["age"] <- hstry$star_age / 10**9 #pro_header$star_age/10**9
    obs.DF["M_current"] <- hstry$star_mass
    
    obs.DF["mass_cc"] <- hstry$mass_conv_core/pro_header$star_mass
    
    obs.DF["X_c"] <- hstry$center_h1
    obs.DF["Y_c"] <- hstry$center_he3 + hstry$center_he4 
    
    obs.DF["log_center_T"] <- hstry$log_center_T
    obs.DF["log_center_Rho"] <- hstry$log_center_Rho
    obs.DF["log_center_P"] <- hstry$log_center_P
    obs.DF["center_mu"] <- hstry$center_mu
    obs.DF["center_degeneracy"] <- hstry$center_degeneracy
    
    obs.DF["mass_X"] <- pro_header$star_mass_h1/pro_header$star_mass
    obs.DF["mass_Y"] <- (pro_header$star_mass_he3 + 
            pro_header$star_mass_he4)/pro_header$star_mass
    
    obs.DF["X_surf"] <- hstry$surface_h1
    obs.DF["Y_surf"] <- hstry$surface_he4 + hstry$surface_he3
    obs.DF["C_surf"] <- hstry$surface_c12
    obs.DF["N_surf"] <- hstry$surface_n14
    obs.DF["O_surf"] <- hstry$surface_o16
    
    obs.DF["log_LH"] <- hstry$log_LH
    obs.DF["log_LHe"] <- hstry$log_LHe
    
    obs.DF["radius"] <- pro_header$photosphere_r
    obs.DF["L"] <- pro_header$photosphere_L
    obs.DF["Teff"] <- pro_header$Teff
    obs.DF["log_Teff"] <- hstry$log_Teff
    obs.DF["log_L"] <- hstry$log_L
    obs.DF["log_R"] <- hstry$log_R
    obs.DF["log_g"] <- hstry$log_g
    
    obs.DF["Fe/H"] <- log10(10**hstry$log_surf_cell_z / 
            hstry$surface_h1 / Z_div_X_solar)
    
    obs.DF["cz_base"] <- hstry$cz_bot_radius
    obs.DF["acoustic_cutoff"] <- hstry$acoustic_cutoff
    obs.DF["acoustic_radius"] <- hstry$acoustic_radius
    
    obs.DF["surface_mu"] <- pro_body$mu
    obs.DF["delta_Pg_asym"] <- hstry$delta_Pg
    obs.DF["nu_max_classic"] <- scaling_nu_max(R=obs.DF[["radius"]], 
        M=hstry[["star_mass"]], Teff=obs.DF[["Teff"]])
    obs.DF["nu_max"] <- scaling_nu_max_Viani(R=obs.DF[["radius"]], 
        M=hstry[["star_mass"]], Teff=obs.DF[["Teff"]], 
        mu=obs.DF[["surface_mu"]]) 
    obs.DF["delta_nu_asym"] <- hstry$delta_nu
    
    freqs <- parse_freqs(freqs_file, gyre=T)
    obs.DF["Dnu0_classic"] <- seismology(freqs, 
        nu_max=obs.DF[["nu_max_classic"]])[["Dnu0"]]
    seis.DF <- seismology(freqs, nu_max=obs.DF[["nu_max"]])
    
    #as.data.frame(cbind(obs.DF, seis.DF))
    merge(rbind(obs.DF), rbind(seis.DF))
}

### Obtain evolutionary tracks from a MESA directory
parse_dir <- function(directory, min_num_models=10, dname='simulations') {
    ## parse dirname string e.g. "M=1.0_Y=0.28"
    params.DF <- NULL
    for (var in unlist(strsplit(basename(directory), '_'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params.DF[nameval[1]] <- as.numeric(nameval[2])
    }
    
    ## obtain history
    dir_files <- list.files(directory)
    log_dirs <- dir_files[grep('LOGS_', dir_files)]
    if (length(log_dirs) == 0) {
        #print("No logs found!")
        #return(data.frame())
        stop("No logs found!")
    }
    
    track <- data.frame()
    for (log_dir in file.path(directory, log_dirs)) {
        logs <- list.files(log_dir)
        if (length(logs) <= 1) {
            print(paste(log_dir, "No logs found!"))
            next 
        }
        ev.DF <- read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5)
        
        ## figure out which profiles & frequency files to use
        profile_candidates <- logs[grep(profile.pattern, logs)]
        freq_file_candidates <- logs[grep(freqs.pattern, logs)]
        pro_files <- c()
        freq_files <- c()
        for (pro_file in profile_candidates) {
            freq_name <- sub(".data", "-freqs.dat", pro_file, fixed=TRUE)
            if (freq_name %in% freq_file_candidates) {
                pro_files <- c(pro_files, pro_file)
                freq_files <- c(freq_files, freq_name)
            }
        }
        if (length(pro_files) <= min_num_models) {
            print(paste(directory, "Too few profile files"))
            next 
        }
        
        ## call summarize on all pairs 
        parallelStartMulticore(max(1,
            as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
        obs.DF <- do.call(plyr:::rbind.fill, 
            parallelMap(function(pro_file, freqs_file)
            #Map(function(pro_file, freqs_file)
                    summarize(pro_file, freqs_file, ev.DF, dname=dname), 
                pro_file=file.path(log_dir, pro_files), 
                freqs_file=file.path(log_dir, freq_files)))
        #print(obs.DF)
        DF <- merge(rbind(params.DF), obs.DF[order(obs.DF$age),])
        
        ## if main sequence, crop PMS
        #if (grepl('LOGS_MS', log_dir)) {
        #    decreasing_L <- which(diff(DF$L) < 0 & DF$center_h1[-1] > 0.55)
        #    if (any(decreasing_L)) {
        #        goes_back_up <- diff(decreasing_L) > 1
        #        pms <- max(decreasing_L)
        #        print(paste(fgong.dir, "Clipping", pms, "points"))
        #        DF <- DF[-1:-pms,]
        #    }
        #}
        
        ## solve linear transport problem to get equally-spaced points 
        space_var <- ifelse(grepl('LOGS_MS', log_dir), 'X_c', 'age')
        x <- DF[[space_var]]
        nrow.DF <- length(x)
        if (nrow.DF < num_points) {
            print(paste(filename, "has too few points"))
            return(NULL)
        }
        ideal <- seq(max(x), min(x), length=num_points)
        cost.mat  <- outer(ideal, x, function(x, y) abs(x-y))
        row.signs <- rep("==", num_points)
        row.rhs   <- rep(1, num_points)
        col.signs <- rep("<=", nrow.DF)
        col.rhs   <- rep(1, nrow.DF)
        sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
            col.signs, col.rhs)$solution
        new.DF <- DF[apply(sol, 1, which.max),]
        
        new.DF$ev.stage <- rep(ev.stages[basename(log_dir)][[1]], nrow(new.DF))
        new.DF <- new.DF[order(new.DF$age),]
        new.DF$ev.tau <- new.DF$ev.stage + ((1:nrow(new.DF))-1)/nrow(new.DF)
        track <- plyr:::rbind.fill(track, new.DF)
    }
    
    #track <- track[order(track$age),]
    
    #PMS_age <- min(track['age'])
    #track['age'] <- track['age'] - PMS_age # remove PMS age 
    #tams_age <- with(track[track['X_c']>=(X_c.lim/100),],
    #    splinefun(X_c, age)(X_c.lim))
    
    #tau_MS <- track['age']/tams_age
    #names(tau_MS) <- 'tau_MS'
    #cbind(track, tau_MS)
    
    track[order(track$age),]
}

### Hertzsprung-Russell diagram
plot_HR <- function(DF, ev.DF, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp,
        col.pal=colorRampPalette(c(blue, 'black', red))(21)
        #brewer.pal(11, "Spectral"))(21)
        ) {
    plot(#log10(DF$Teff), log10(DF$L), 
        #ev.DF$log_L ~ ev.DF$log_Teff,
        #type='l', 
        NA,
        tcl=0, axes=FALSE,
        xlab=expression('Temperature' ~ 'lg'*(T[eff]/K)),
        ylab=expression('Luminosity' ~ 'lg'*(L / L['\u0298'])),
        xlim=rev(log10(range(DF$Teff))),
        ylim=log10(range(DF$L)))
    #abline(v=log10(5771), lty=3, col='lightgray')
    abline(v=log10(7000), lty=2, col='gray')
    #abline(h=0, lty=3, col='lightgray')
    points(log10(5777), 0, pch=1, cex=1)
    points(log10(5777), 0, pch=20, cex=0.1)
    lines(ev.DF$log_Teff, ev.DF$log_L, lty=1, col='gray')
    points(log10(DF$Teff), log10(DF$L), cex=0.2, 
        pch=ifelse(DF$ev.stage%%2, 20, 3),
        #col=col.pal[floor((DF$X_c-min(DF$X_c)) / (max(DF$X_c)-min(DF$X_c))
        col=col.pal[floor((DF$age-min(DF$age)) / (max(DF$age)-min(DF$age))
                          * (length(col.pal)-1)) + 1])
    magaxis(side=1:4, family=font, tcl=0.25, mgp=utils.mgp, las=1,
        cex.axis=text.cex, labels=c(1,1,0,0)) 
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        round(quantile(seq(min(DF$age), max(DF$age),#min(DF$X_c), max(DF$X_c), 
            length=length(col.pal)), c(0, 0.25, 0.5, 0.75, 1)), 3), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(expression("Age" ~ tau/"Gyr"), 4, 
        line=5.5, cex=text.cex)
    par(family="Luxi Mono")
    legend("bottomright", bty='n', xjust=1, cex=text.cex/2, col="gray", 
        legend=c(
            as.expression(bquote(M == .(DF$M[1]))), 
            as.expression(bquote(Y == .(DF$Y[1]))), 
            as.expression(bquote(Z == .(DF$Z[1])))
        )
    )
}

###############################################################################
### Parse command line arguments and process directory ########################
###############################################################################
args <- commandArgs(TRUE)
if (length(args)>0) {
    directory <- args[1]
    print(directory)
    dname <- dirname(directory)
    parsed_dir <- parse_dir(directory, dname=dname)
    
    
    
    DF <- unique(parsed_dir)
    DF <- DF[order(DF$age),]
    
    print(paste("Produced", nrow(DF), "summaries"))
    print(sapply(DF, fivenum))
    
    ## check if the track should be kept 
    if (nrow(DF) < 5) {
        print("Rejecting track: too few data points")
        print(DF)
    } else if (ncol(DF) < 5) {
        print("Rejecting track: too few columns")
        print(DF)
    } else { # Make a table of results! 
        print(paste('Writing table to', paste0(directory, '.dat')))
        write.table(DF, paste0(directory, '.dat'), quote=FALSE, sep='\t', 
            row.names=FALSE)
        
        ## make plots
        ev.DF <- data.frame()
        dir_files <- list.files(directory)
        log_dirs <- file.path(directory, dir_files[grep('LOGS_', dir_files)])
        for (log_dir in log_dirs) {
            ev.DF <- plyr:::rbind.fill(ev.DF, 
                read.table(file.path(log_dir, 'history.data'), 
                    header=TRUE, skip=5))
            ev.DF <- ev.DF[order(ev.DF$star_age),]
            make_plots(plot_HR, paste0(basename(directory), "-HR"), 
                filepath=file.path('plots', dirname(directory), 'HR'), 
                mar=c(3, 4, 1, 7), DF=DF, ev.DF=ev.DF)
        }
    }
    
    warnings()
}

