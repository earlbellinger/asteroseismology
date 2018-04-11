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
library(RColorBrewer)

#options(warn=2) 

Z_div_X_solar = 0.02293
profile.pattern <- 'profile.+.data$'
freqs.pattern <- 'profile.+-freqs.dat$'
freqs.cols <- c('l', 'n', 'nu')#, 'inertia')
num_points <- 64
X_c.lim <- 1e-2
ev_stages <- list(
    'LOGS_MS'=1, 
    'LOGS_SG'=2, 
    'LOGS_RGB'=3, 
    'LOGS_BUMP'=4, 
    'LOGS_HEB'=5)

exh <- 1e-5

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
    
    obs.DF["mass_cc"] <- hstry$mass_conv_core/hstry$star_mass
    obs.DF["mass_he_core"] <- hstry$he_core_mass
    
    obs.DF["h_exh_core_mass"] <- hstry$h_exh_core_mass
    obs.DF["h_exh_core_radius"] <- hstry$h_exh_core_radius
    
    obs.DF["radius_cc"] <- if (obs.DF['mass_cc'] <= 0) 0 else {
        if (hstry$conv_mx1_bot_r <= 0.1) {
            hstry$conv_mx1_top_r / 10**(hstry$log_R)
        } else {
            hstry$conv_mx2_top_r / 10**(hstry$log_R)
        }
    }
    
    
    obs.DF["mass_X"] <- pro_header$star_mass_h1/pro_header$star_mass
    obs.DF["mass_Y"] <- (pro_header$star_mass_he3 + 
            pro_header$star_mass_he4)/pro_header$star_mass
    
    obs.DF["log_LH"] <- hstry$log_LH
    obs.DF["log_LHe"] <- hstry$log_LHe
    obs.DF["log_center_T"] <- hstry$log_center_T
    obs.DF["log_center_Rho"] <- hstry$log_center_Rho
    obs.DF["log_center_P"] <- hstry$log_center_P
    obs.DF["center_mu"] <- hstry$center_mu
    obs.DF["center_degeneracy"] <- hstry$center_degeneracy
    
    obs.DF["X_c"]  <- hstry$center_h1 + hstry$center_h2
    obs.DF["Y_c"]  <- hstry$center_he3 + hstry$center_he4 
    obs.DF["Li_c"] <- hstry$center_li7
    obs.DF["Be_c"] <- hstry$center_be7
    obs.DF["B_c"]  <- hstry$center_b8
    obs.DF["C_c"]  <- hstry$center_c12 + hstry$center_c13
    obs.DF["N_c"]  <- hstry$center_n13 + hstry$center_n14 + hstry$center_n15
    obs.DF["O_c"]  <- hstry$center_o14 + hstry$center_o15 + hstry$center_o16 +
        hstry$center_o17 + hstry$center_o18
    obs.DF["F_c"]  <- hstry$center_f17 + hstry$center_f18 + hstry$center_f19
    obs.DF["Ne_c"] <- hstry$center_ne18 + hstry$center_ne19 + 
        hstry$center_ne20 + hstry$center_ne22
    obs.DF["Mg_c"] <- hstry$center_mg22 + hstry$center_mg24
    
    obs.DF["H1_c"]   <- hstry$center_h1 
    obs.DF["H2_c"]   <- hstry$center_h2 
    obs.DF["He3_c"]  <- hstry$center_he3 
    obs.DF["He4_c"]  <- hstry$center_he4 
    obs.DF["Li7_c"]  <- hstry$center_li7 
    obs.DF["Be7_c"]  <- hstry$center_be7 
    obs.DF["B8_c"]   <- hstry$center_b8 
    obs.DF["C12_c"]  <- hstry$center_c12 
    obs.DF["C13_c"]  <- hstry$center_c13 
    obs.DF["N13_c"]  <- hstry$center_n13 
    obs.DF["N14_c"]  <- hstry$center_n14 
    obs.DF["N15_c"]  <- hstry$center_n15 
    obs.DF["O14_c"]  <- hstry$center_o14 
    obs.DF["O15_c"]  <- hstry$center_o15 
    obs.DF["O16_c"]  <- hstry$center_o16 
    obs.DF["O17_c"]  <- hstry$center_o17 
    obs.DF["O18_c"]  <- hstry$center_o18 
    obs.DF["F17_c"]  <- hstry$center_f17 
    obs.DF["F18_c"]  <- hstry$center_f18 
    obs.DF["F19_c"]  <- hstry$center_f19 
    obs.DF["Ne18_c"] <- hstry$center_ne18 
    obs.DF["Ne19_c"] <- hstry$center_ne19 
    obs.DF["Ne20_c"] <- hstry$center_ne20 
    obs.DF["Ne22_c"] <- hstry$center_ne22 
    obs.DF["Mg22_c"] <- hstry$center_mg22 
    obs.DF["Mg24_c"] <- hstry$center_mg24 
    
    obs.DF["X_surf"]  <- hstry$surface_h1 + hstry$surface_h2
    obs.DF["Y_surf"]  <- hstry$surface_he3 + hstry$surface_he4 
    obs.DF["Li_surf"] <- hstry$surface_li7
    obs.DF["Be_surf"] <- hstry$surface_be7
    obs.DF["B_surf"]  <- hstry$surface_b8
    obs.DF["C_surf"]  <- hstry$surface_c12 + hstry$surface_c13
    obs.DF["N_surf"]  <- hstry$surface_n13 + hstry$surface_n14 + 
        hstry$surface_n15
    obs.DF["O_surf"]  <- hstry$surface_o14 + hstry$surface_o15 + 
        hstry$surface_o16 + hstry$surface_o17 + hstry$surface_o18
    obs.DF["F_surf"]  <- hstry$surface_f17 + hstry$surface_f18 + 
        hstry$surface_f19
    obs.DF["Ne_surf"] <- hstry$surface_ne18 + hstry$surface_ne19 + 
        hstry$surface_ne20 + hstry$surface_ne22
    obs.DF["Mg_surf"] <- hstry$surface_mg22 + hstry$surface_mg24
    
    obs.DF["H1_surf"]   <- hstry$surface_h1 
    obs.DF["H2_surf"]   <- hstry$surface_h2 
    obs.DF["He3_surf"]  <- hstry$surface_he3 
    obs.DF["He4_surf"]  <- hstry$surface_he4 
    obs.DF["Li7_surf"]  <- hstry$surface_li7 
    obs.DF["Be7_surf"]  <- hstry$surface_be7 
    obs.DF["B8_surf"]   <- hstry$surface_b8 
    obs.DF["C12_surf"]  <- hstry$surface_c12 
    obs.DF["C13_surf"]  <- hstry$surface_c13 
    obs.DF["N13_surf"]  <- hstry$surface_n13 
    obs.DF["N14_surf"]  <- hstry$surface_n14 
    obs.DF["N15_surf"]  <- hstry$surface_n15 
    obs.DF["O14_surf"]  <- hstry$surface_o14 
    obs.DF["O15_surf"]  <- hstry$surface_o15 
    obs.DF["O16_surf"]  <- hstry$surface_o16 
    obs.DF["O17_surf"]  <- hstry$surface_o17 
    obs.DF["O18_surf"]  <- hstry$surface_o18 
    obs.DF["F17_surf"]  <- hstry$surface_f17 
    obs.DF["F18_surf"]  <- hstry$surface_f18 
    obs.DF["F19_surf"]  <- hstry$surface_f19 
    obs.DF["Ne18_surf"] <- hstry$surface_ne18 
    obs.DF["Ne19_surf"] <- hstry$surface_ne19 
    obs.DF["Ne20_surf"] <- hstry$surface_ne20 
    obs.DF["Ne22_surf"] <- hstry$surface_ne22 
    obs.DF["Mg22_surf"] <- hstry$surface_mg22 
    obs.DF["Mg24_surf"] <- hstry$surface_mg24 
    
    obs.DF["radius"] <- pro_header$photosphere_r
    obs.DF["L"] <- pro_header$photosphere_L
    obs.DF["Teff"] <- pro_header$Teff
    obs.DF["log_Teff"] <- hstry$log_Teff
    obs.DF["log_L"] <- hstry$log_L
    obs.DF["log_R"] <- hstry$log_R
    obs.DF["log_g"] <- hstry$log_g
    
    obs.DF["Fe_H"] <- log10(10**hstry$log_surf_cell_z / 
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
    seis <- seismology(freqs, nu_max=obs.DF[["nu_max_classic"]],
        min_points=3, check_nu_max=T)
    if ("Dnu0" %in% names(seis)) obs.DF["Dnu0_classic"] <- seis[["Dnu0"]]
    seis.DF <- seismology(freqs, nu_max=obs.DF[["nu_max"]],
        min_points=3, check_nu_max=T)
    
    #as.data.frame(cbind(obs.DF, seis.DF))
    merge(rbind(obs.DF), rbind(seis.DF))
}

### Obtain evolutionary tracks from a MESA directory
parse_dir <- function(directory, min_num_models=10, dname='simulations',
        num_points=num_points) {
    ## parse dirname string e.g. "M=1.0_Y=0.28"
    params.DF <- NULL
    trackfile <- file.path(directory, 'track')
    if (file.exists(trackfile)) {
        params.DF <- read.table(trackfile, header=1)
    } else {
        for (var in unlist(strsplit(basename(directory), '_'))) { 
            nameval <- unlist(strsplit(var, "=")) 
            params.DF[nameval[1]] <- as.numeric(nameval[2])
        }
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
        
        DF$lg_X_c <- log10(DF$X_c) #center_h1)
        DF$lg_Y_c <- log10(DF$Y_c) #center_he3 + DF$center_he4)
        ## solve linear transport problem to get equally-spaced points 
        space_var <- ifelse(grepl('LOGS_MS', log_dir), 'lg_X_c', 
                     ifelse(grepl('LOGS_HEB', log_dir), 'lg_Y_c', 
                     #ifelse(grepl('LOGS_BUMP', log_dir), 'center_degeneracy', 
                        'age'))
        x <- DF[[space_var]]
        nrow.DF <- length(x)
        if (nrow.DF < num_points) {
            print(paste(log_dir, "has too few points"))
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
        
        if (grepl('LOGS_BUMP', log_dir)) {
            new.DF <- rbind(DF[which.max(DF$L),], new.DF)#[-1,])
        }
        
        new.DF$ev_stage <- rep(ev_stages[basename(log_dir)][[1]], nrow(new.DF))
        new.DF <- new.DF[order(new.DF$age),]
        new.DF$ev_tau <- new.DF$ev_stage + ((1:nrow(new.DF))-1)/nrow(new.DF)
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
plot_HR <- function(DF, ev.DF, plot_legend=T, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp,
        col.pal=brewer.pal(5, "Spectral")
        #colorRampPalette(c(blue, 'black', red))(21)
        #brewer.pal(11, "Spectral"))(21)
        ) {
    plot(#log10(DF$Teff), log10(DF$L), 
        #ev.DF$log_L ~ ev.DF$log_Teff,
        #type='l', 
        NA,
        tcl=0, axes=FALSE,
        xlab="",
        ylab="",
        xlim=rev(range(ev.DF$log_Teff)),
        ylim=range(ev.DF$log_L))
    #abline(v=log10(5771), lty=3, col='lightgray')
    abline(v=log10(7000), lty=2, lwd=2, col='gray')
    #abline(h=0, lty=3, col='lightgray')
    lines(ev.DF$log_Teff, ev.DF$log_L, lty=1, lwd=2, col='gray')
    points(log10(DF$Teff), log10(DF$L), cex=0.1, 
        pch=20, lwd=1, #ifelse(DF$ev_stage%%2, 20, 3),
        col=col.pal[DF$ev_stage])
        #col=col.pal[floor((DF$X_c-min(DF$X_c)) / (max(DF$X_c)-min(DF$X_c))
        #col=col.pal[floor((DF$age-min(DF$age)) / (max(DF$age)-min(DF$age))
        #                  * (length(col.pal)-1)) + 1])
    magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1,
        cex.axis=text.cex, labels=T) 
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), las=1,
        cex.axis=text.cex, labels=T) 
    magaxis(side=3:4, family=font, tcl=0, mgp=mgp, las=1,
        cex.axis=text.cex, labels=F)
    
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4],
        c("MS", "SG", "RGB", "BUMP", "CLUMP"),
        col.pal, gradient='y', align='rb', cex=text.cex)
    
    par(mgp=mgp+c(-0.5, 0, 0))
    title(xlab=expression('Temperature' ~ 'lg'*(T['eff']/K)))
    par(mgp=mgp+c(-0.7, 0, 0))
    title(ylab=expression('Luminosity' ~ 'lg'*(L / L['solar'])))
    
    if (plot_legend) {
        par(family="Luxi Mono")
        legend("bottom", xjust=1, cex=text.cex/2, col="gray", inset=c(0.03, 0.03),
            legend=c(
                as.expression(bquote(M == .(signif(DF$M[1],3)))), 
                as.expression(bquote(Y == .(signif(DF$Y[1],3)))), 
                as.expression(bquote(Z == .(signif(DF$Z[1],3)))), 
                as.expression(bquote(alpha["MLT"] == .(signif(DF$alpha[1],3)))), 
                as.expression(bquote(alpha["ov"] == .(signif(DF$ov[1],3)))), 
                as.expression(bquote(D == .(signif(DF$diffusion[1],3)))),
                as.expression(bquote(G == .(signif(DF$settling[1],3)))),
                as.expression(bquote(eta == .(signif(DF$eta[1],3))))
            )
        )
    }
    
    points(log10(5777), 0, pch=20, cex=0.1, lwd=1.5)
    points(log10(5777), 0, pch=1,  cex=0.8, lwd=1.5)
}



plot_chem_c <- function(DF, ev.DF, plot_legend=T, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    
    col.pal <- brewer.pal(7, "Spectral")
    col.pal[4] <- '#7851A9'
    col.pal[5] <- col.pal[6]
    col.pal[6] <- '#CA1F7B'
    tmp <- col.pal[7]
    col.pal[7] <- col.pal[5]
    col.pal[5] <- tmp
    
    base <- 2
    ages <- base**(ev.DF$star_age/10**9)
    ages2 <- base**DF$age
    
    labs <- pretty(log(ages, base))
    labs[1] <- 0
    yticks <- base**labs
    yticks[1] <- 0
    
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(0, max(ages)), #base**max(ceiling(ev.DF$star_age/(2*10**8))/5) ),
        ylim=log10(c(10**-6, 1.5)))
    
    lines(ages, log10(ev.DF$center_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ages, log10(ev.DF$center_ne20 + ev.DF$center_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ages, log10(ev.DF$center_o16 + ev.DF$center_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ages, log10(ev.DF$center_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ages, log10(ev.DF$center_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ages, log10(ev.DF$center_he3 + ev.DF$center_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ages, log10(ev.DF$center_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    if (F) {
    points(ages2, log10(DF$Mg_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[7])
    points(ages2, log10(DF$Ne_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[6])
    points(ages2, log10(DF$O_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[5])
    points(ages2, log10(DF$N_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[4])
    points(ages2, log10(DF$C_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[3])
    points(ages2, log10(DF$Y_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[2])
    points(ages2, log10(DF$X_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[1])
    }
    
    stages <- base**DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages, lwd=1, lty=3)
    abline(v=base**0, lwd=1, lty=3)
    abline(v=base**max(DF$age), lwd=1, lty=3)
    abline(h=log10(exh), lwd=1, lty=2)
    off <- max(yticks)/200
    text(base**0-off*1.5,   log10(9e-6), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[1]-off*1.8, log10(0.2),  labels="TAMS", pos=4,          cex=text.cex/1.66)
    text(stages[2]-off*1.8, log10(0.2),  labels="RGB",  pos=4,          cex=text.cex/1.66)
    text(stages[3]-off*1.5, log10(0.3),  labels="BUMP", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[4]-off*1.5, log10(5e-6), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4], log10(0.5), labels="RC", pos=4)
    
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
        labels=labs, tcl=0)
    axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 1), 
        tcl=-0.25, labels=F)
    axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 0.2), 
        tcl=-0.125, labels=F)
    
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(1)) 
    #magaxis(side=3:4, tcl=0, cex.axis=text.cex, labels=c(F, F), mgp=mgp)
    
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4],
        rev(c("H", "He", "C", "N", "O", "Ne", "Mg")),
        rev(col.pal), gradient='y', align='rb', cex=text.cex)
    
    par(mgp=mgp+c(-0.5, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    par(mgp=mgp+c(0.3, 0, 0))
    title(ylab=expression('Central Abundance'~X["c"]))
    
    if (plot_legend) {
        par(family="Luxi Mono")
        legend("bottom", xjust=1, cex=text.cex/2, col="gray", inset=c(0.03, 0.03),
            legend=c(
                as.expression(bquote(M == .(signif(DF$M[1],3)))), 
                as.expression(bquote(Y == .(signif(DF$Y[1],3)))), 
                as.expression(bquote(Z == .(signif(DF$Z[1],3)))), 
                as.expression(bquote(alpha["MLT"] == .(signif(DF$alpha[1],3)))), 
                as.expression(bquote(alpha["ov"] == .(signif(DF$ov[1],3)))), 
                as.expression(bquote(D == .(signif(DF$diffusion[1],3)))),
                as.expression(bquote(G == .(signif(DF$settling[1],3)))),
                as.expression(bquote(eta == .(signif(DF$eta[1],3))))
            )
        )
    }
}




plot_chem_c2 <- function(DF, ev.DF, plot_legend=T, text.cex=1, ..., 
        font=utils.font, mgp=utils.mgp) {
    
    
    col.pal <- brewer.pal(7, "Spectral")
    col.pal[4] <- '#7851A9'
    col.pal[5] <- col.pal[6]
    col.pal[6] <- '#CA1F7B'
    tmp <- col.pal[7]
    col.pal[7] <- col.pal[5]
    col.pal[5] <- tmp
    #col.pal[3] <- 
    
    base <- 2
    ages <- (ev.DF$star_age/10**9)#base**(ev.DF$star_age/10**9)
    ages2 <- DF$age#base**DF$age
    
    change.age <- DF$age[as.logical(diff(DF$ev_stage))][2]
    if (is.na(change.age)) change.age <- max(DF$age)
    
    labs <- pretty(log(ages, base))
    labs[1] <- 0
    yticks <- labs#base**labs
    yticks[1] <- 0
    
    par(mfrow=c(1,2), mar=c(2.5, 3, 0.5, 0.25))
    
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(0, change.age), 
        ylim=log10(c(exh, 1.5)))
    
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_ne20 + ev.DF$center_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_o16 + ev.DF$center_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_he3 + ev.DF$center_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ev.DF$star_age/10**9, log10(ev.DF$center_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    stages <- DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages[1], lwd=1, lty=3)
    abline(v=0, lwd=1, lty=3)
    #abline(v=change.age, lwd=1, lty=3)
    #abline(h=log10(exh), lwd=1, lty=2)
    off <- change.age/18
    text(0-off,         log10(0.1), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.33)
    text(stages[1]-off, log10(0.1), labels="TAMS", pos=4, srt=-90, cex=text.cex/1.33)
    #text(stages[2]-off, log10(0.2),  labels="RGB",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[3]-off, log10(0.3),  labels="BUMP", pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4]-off, log10(5e-6), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4], log10(0.5), labels="RC", pos=4)
    
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    magaxis(side=1, family=font, tcl=-0.25, mgp=mgp+c(0, -0.1, 0), 
        las=1, cex.axis=text.cex, labels=c(1)) 
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(1))  
    magaxis(side=4, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(F)) 
    par(mgp=mgp+c(-0.3, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    #magaxis(side=3:4, tcl=0, cex.axis=text.cex, labels=c(F, F), mgp=mgp)
    par(mgp=mgp+c(0.3, 0, 0))
    title(ylab=expression('Central Abundance'~X["c"]))
    
    
    
    change.age <- DF$age[as.logical(diff(DF$ev_stage))][3]
    par(mar=c(2.5, 0.75, 0.5, 2.5))
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(change.age, max(ages)), 
        ylim=log10(c(exh, 1.5)))
    
    lines(ages, log10(ev.DF$center_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ages, log10(ev.DF$center_ne20 + ev.DF$center_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ages, log10(ev.DF$center_o16 + ev.DF$center_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ages, log10(ev.DF$center_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ages, log10(ev.DF$center_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ages, log10(ev.DF$center_he3 + ev.DF$center_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ages, log10(ev.DF$center_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    if (F) {
    points(ages2, log10(DF$Mg_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[7])
    points(ages2, log10(DF$Ne_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[6])
    points(ages2, log10(DF$O_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[5])
    points(ages2, log10(DF$N_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[4])
    points(ages2, log10(DF$C_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[3])
    points(ages2, log10(DF$Y_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[2])
    points(ages2, log10(DF$X_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[1])
    }
    
    stages <- DF$age[as.logical(diff(DF$ev_stage))]#base**DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages[-1], lwd=1, lty=3)
    #abline(v=change.age, lwd=1, lty=3)
    abline(v=max(DF$age), lwd=1, lty=3)
    #abline(h=log10(exh), lwd=1, lty=2)
    off <- max(ages)/900
    #text(0-off,   log10(9e-6), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[1]-off, log10(0.2),  labels="TAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[2]-off, log10(0.1),    labels="RGB",  pos=4, srt=-90, cex=text.cex/1.33)
    text(stages[3]-off, log10(0.1),    labels="BUMP", pos=4, srt=-90, cex=text.cex/1.33)
    text(stages[4]-off, log10(3.5e-5), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.33)
    #text(stages[4], log10(0.5), labels="RC", pos=4)
    
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    magaxis(side=1, family=font, tcl=-0.25, mgp=mgp+c(0, -0.1, 0), majorn=3, 
        las=1, cex.axis=text.cex, labels=c(T)) 
    #axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
    #    labels=labs, tcl=0)
    #axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 1), 
    #    tcl=-0.25, labels=F)
    #axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 0.2), 
    #    tcl=-0.125, labels=F)
    
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(F)) 
    #magaxis(side=3:4, tcl=0, cex.axis=text.cex, labels=c(F, F), mgp=mgp)
    
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.15*var1range, par()$usr[4],
        rev(c("H", "He", "C", "N", "O", "Ne", "Mg")),
        rev(col.pal), gradient='y', align='rb', cex=text.cex)
    
    par(mgp=mgp+c(-0.3, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    
    if (plot_legend) {
        par(family="Luxi Mono")
        legend("bottom", xjust=1, cex=text.cex/2, col="gray", inset=c(0.03, 0.03),
            legend=c(
                as.expression(bquote(M == .(signif(DF$M[1],3)))), 
                as.expression(bquote(Y == .(signif(DF$Y[1],3)))), 
                as.expression(bquote(Z == .(signif(DF$Z[1],3)))), 
                as.expression(bquote(alpha["MLT"] == .(signif(DF$alpha[1],3)))), 
                as.expression(bquote(alpha["ov"] == .(signif(DF$ov[1],3)))), 
                as.expression(bquote(D == .(signif(DF$diffusion[1],3)))),
                as.expression(bquote(G == .(signif(DF$settling[1],3)))),
                as.expression(bquote(eta == .(signif(DF$eta[1],3))))
            )
        )
    }
}


plot_chem_surf <- function(DF, ev.DF, plot_legend=T, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    
    col.pal <- brewer.pal(7, "Spectral")
    col.pal[4] <- '#7851A9'
    col.pal[5] <- col.pal[6]
    col.pal[6] <- '#CA1F7B'
    tmp <- col.pal[7]
    col.pal[7] <- col.pal[5]
    col.pal[5] <- tmp
    
    base <- 2
    ages <- base**(ev.DF$star_age/10**9)
    ages2 <- base**DF$age
    
    labs <- pretty(log(ages, base))
    labs[1] <- 0
    yticks <- base**labs
    yticks[1] <- 0
    
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(0, max(ages)),
        ylim=log10(c(10**-4, 1.2)))
    
    lines(ages, log10(ev.DF$surface_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ages, log10(ev.DF$surface_ne20 + ev.DF$surface_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ages, log10(ev.DF$surface_o16 + ev.DF$surface_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ages, log10(ev.DF$surface_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ages, log10(ev.DF$surface_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ages, log10(ev.DF$surface_he3 + ev.DF$surface_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ages, log10(ev.DF$surface_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    if (F) {
    points(ages2, log10(DF$Mg_surf), pch=20, cex=0.5, lwd=0.1, col=col.pal[7])
    points(ages2, log10(DF$Ne_surf), pch=20, cex=0.5, lwd=0.1, col=col.pal[6])
    points(ages2, log10(DF$O_surf),  pch=20, cex=0.5, lwd=0.1, col=col.pal[5])
    points(ages2, log10(DF$N_surf),  pch=20, cex=0.5, lwd=0.1, col=col.pal[4])
    points(ages2, log10(DF$C_surf),  pch=20, cex=0.5, lwd=0.1, col=col.pal[3])
    points(ages2, log10(DF$Y_surf),  pch=20, cex=0.5, lwd=0.1, col=col.pal[2])
    points(ages2, log10(DF$X_surf),  pch=20, cex=0.5, lwd=0.1, col=col.pal[1])
    }
    
    stages <- base**DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages, lwd=1, lty=3)
    abline(v=base**0, lwd=1, lty=3)
    abline(v=base**max(DF$age), lwd=1, lty=3)
    off <- max(yticks)/200
    text(base**0-off*1.5,   log10(0.06), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[1]-off*1.5, log10(0.06), labels="TAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[2]-off*1.5, log10(0.05), labels="RGB",  pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[3]-off*1.5, log10(0.06), labels="BUMP", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[4]-off*1.5, log10(0.04), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.66)
    
    #text(stages[1]-off, log10(0.03), labels="SG", pos=4, cex=text.cex/1.3)
    #text(stages[2]-off, log10(0.03), labels="RGB", pos=4, cex=text.cex/1.3)
    #text(stages[3]-off*1.5, log10(0.08), labels="BUMP", pos=4, srt=-90, cex=text.cex/2)
    #text(stages[4]-off*2, log10(0.04), labels="HeB", pos=4, cex=text.cex/2)
    ##text(stages[4], log10(0.5), labels="RC", pos=4)
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
        labels=labs, tcl=0)
    axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 1), 
        tcl=-0.25, labels=F)
    axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 0.2), 
        tcl=-0.125, labels=F)
    #axis(1, tick=T, at=base**seq(0, max(log(ages2, base)), 0.2), 
    #    tcl=-0.125, labels=F)
    
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(1)) 
    
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4],
        rev(c("H", "He", "C", "N", "O", "Ne", "Mg")),
        rev(col.pal), gradient='y', align='rb', cex=text.cex)
    
    par(mgp=mgp+c(-0.5, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    par(mgp=mgp+c(0.3, 0, 0))
    title(ylab=expression('Surface Abundance'~X["surf"]))
    
    if (plot_legend) {
        par(family="Luxi Mono")
        legend("bottom", xjust=1, cex=text.cex/2, col="gray", inset=c(0.03, 0.03),
            legend=c(
                as.expression(bquote(M == .(signif(DF$M[1],3)))), 
                as.expression(bquote(Y == .(signif(DF$Y[1],3)))), 
                as.expression(bquote(Z == .(signif(DF$Z[1],3)))), 
                as.expression(bquote(alpha["MLT"] == .(signif(DF$alpha[1],3)))), 
                as.expression(bquote(alpha["ov"] == .(signif(DF$ov[1],3)))), 
                as.expression(bquote(D == .(signif(DF$diffusion[1],3)))),
                as.expression(bquote(G == .(signif(DF$settling[1],3)))),
                as.expression(bquote(eta == .(signif(DF$eta[1],3))))
            )
        )
    }
}




plot_chem_surf2 <- function(DF, ev.DF, plot_legend=T, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    
    col.pal <- brewer.pal(7, "Spectral")
    col.pal[4] <- '#7851A9'
    col.pal[5] <- col.pal[6]
    col.pal[6] <- '#CA1F7B'
    tmp <- col.pal[7]
    col.pal[7] <- col.pal[5]
    col.pal[5] <- tmp
    
    base <- 2
    ages <- (ev.DF$star_age/10**9)#base**(ev.DF$star_age/10**9)
    ages2 <- DF$age#base**DF$age
    
    change.age <- DF$age[as.logical(diff(DF$ev_stage))][2]
    
    labs <- pretty(log(ages, base))
    labs[1] <- 0
    yticks <- labs#base**labs
    yticks[1] <- 0
    
    par(mfrow=c(1,2), mar=c(2, 3, 0.5, 0.25))
    
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(0, change.age), 
        ylim=log10(c(exh, 1.5)))
    
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_ne20 + ev.DF$surface_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_o16 + ev.DF$surface_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_he3 + ev.DF$surface_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ev.DF$star_age/10**9, log10(ev.DF$surface_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    stages <- DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages[1], lwd=1, lty=3)
    abline(v=0, lwd=1, lty=3)
    #abline(v=change.age, lwd=1, lty=3)
    #abline(h=log10(exh), lwd=1, lty=2)
    off <- change.age/18
    text(0-off,         log10(0.1), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[1]-off, log10(0.1),  labels="TAMS", pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[2]-off, log10(0.2),  labels="RGB",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[3]-off, log10(0.3),  labels="BUMP", pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4]-off, log10(5e-6), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4], log10(0.5), labels="RC", pos=4)
    
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, 
        las=1, cex.axis=text.cex, labels=c(1)) 
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(1))  
    magaxis(side=4, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(F)) 
    par(mgp=mgp+c(-0.5, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    #magaxis(side=3:4, tcl=0, cex.axis=text.cex, labels=c(F, F), mgp=mgp)
    par(mgp=mgp+c(0.3, 0, 0))
    title(ylab=expression('Surface Abundance'~X["c"]))
    
    
    change.age <- DF$age[as.logical(diff(DF$ev_stage))][3]
    par(mar=c(2, 0.75, 0.5, 2.5))
    plot(NA,
        tcl=0, axes=FALSE, yaxs='i', #xaxs='i', 
        xlab="",
        ylab="",
        xlim=c(change.age, max(ages)), 
        ylim=log10(c(exh, 1.5)))
    
    lines(ages, log10(ev.DF$surface_mg24), lty=1, lwd=3, 
        col=col.pal[7])
    lines(ages, log10(ev.DF$surface_ne20 + ev.DF$surface_ne22), lty=1, lwd=3, 
        col=col.pal[6])
    lines(ages, log10(ev.DF$surface_o16 + ev.DF$surface_o18), lty=1, lwd=3, 
        col=col.pal[5])
    lines(ages, log10(ev.DF$surface_n14), lty=1, lwd=3, 
        col=col.pal[4])
    lines(ages, log10(ev.DF$surface_c12), lty=1, lwd=3, 
        col=col.pal[3])
    lines(ages, log10(ev.DF$surface_he3 + ev.DF$surface_he4), lty=1, lwd=3, 
        col=col.pal[2])
    lines(ages, log10(ev.DF$surface_h1), lty=1, lwd=3, 
        col=col.pal[1])
    
    if (F) {
    points(ages2, log10(DF$Mg_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[7])
    points(ages2, log10(DF$Ne_c), pch=20, cex=0.5, lwd=0.1, col=col.pal[6])
    points(ages2, log10(DF$O_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[5])
    points(ages2, log10(DF$N_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[4])
    points(ages2, log10(DF$C_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[3])
    points(ages2, log10(DF$Y_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[2])
    points(ages2, log10(DF$X_c),  pch=20, cex=0.5, lwd=0.1, col=col.pal[1])
    }
    
    stages <- DF$age[as.logical(diff(DF$ev_stage))]#base**DF$age[as.logical(diff(DF$ev_stage))]
    abline(v=stages[-1], lwd=1, lty=3)
    #abline(v=change.age, lwd=1, lty=3)
    abline(v=max(DF$age), lwd=1, lty=3)
    #abline(h=log10(exh), lwd=1, lty=2)
    off <- max(ages)/900
    #text(0-off,   log10(9e-6), labels="ZAMS", pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[1]-off, log10(0.2),  labels="TAMS", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[2]-off, log10(0.1),  labels="RGB",  pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[3]-off, log10(0.1),  labels="BUMP", pos=4, srt=-90, cex=text.cex/1.66)
    text(stages[4]-off, log10(3.5e-5), labels="TIP",  pos=4, srt=-90, cex=text.cex/1.66)
    #text(stages[4], log10(0.5), labels="RC", pos=4)
    
    #magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, las=1, unlog='y', 
    #    cex.axis=text.cex, labels=c(1)) 
    
    magaxis(side=1:4, tcl=0, cex.axis=text.cex, labels=F, mgp=mgp)
    magaxis(side=1, family=font, tcl=-0.25, mgp=mgp, majorn=3, 
        las=1, cex.axis=text.cex, labels=c(T)) 
    #axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
    #    labels=labs, tcl=0)
    #axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 1), 
    #    tcl=-0.25, labels=F)
    #axis(1, tick=T, at=base**seq(0, max(ev.DF$star_age)/10**9, 0.2), 
    #    tcl=-0.125, labels=F)
    
    magaxis(side=2, family=font, tcl=-0.25, mgp=mgp+c(0, 0.2, 0), 
        las=1, unlog='y', 
        cex.axis=text.cex, labels=c(F)) 
    #magaxis(side=3:4, tcl=0, cex.axis=text.cex, labels=c(F, F), mgp=mgp)
    
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4],
        rev(c("H", "He", "C", "N", "O", "Ne", "Mg")),
        rev(col.pal), gradient='y', align='rb', cex=text.cex)
    
    par(mgp=mgp+c(-0.5, 0, 0))
    title(xlab=expression('Age'~tau/Gyr))
    
    if (plot_legend) {
        par(family="Luxi Mono")
        legend("bottom", xjust=1, cex=text.cex/2, col="gray", inset=c(0.03, 0.03),
            legend=c(
                as.expression(bquote(M == .(signif(DF$M[1],3)))), 
                as.expression(bquote(Y == .(signif(DF$Y[1],3)))), 
                as.expression(bquote(Z == .(signif(DF$Z[1],3)))), 
                as.expression(bquote(alpha["MLT"] == .(signif(DF$alpha[1],3)))), 
                as.expression(bquote(alpha["ov"] == .(signif(DF$ov[1],3)))), 
                as.expression(bquote(D == .(signif(DF$diffusion[1],3)))),
                as.expression(bquote(G == .(signif(DF$settling[1],3)))),
                as.expression(bquote(eta == .(signif(DF$eta[1],3))))
            )
        )
    }
}



###############################################################################
### Parse command line arguments and process directory ########################
###############################################################################
args <- commandArgs(TRUE)
if (length(args)>0) {
    directory <- args[1]
    print(directory)
    dname <- dirname(directory)
    if (length(args)>1) num_points <- as.numeric(args[2])
    parsed_dir <- parse_dir(directory, dname=dname, num_points=num_points)
    
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
        
        if (F) {
        ## make plots
        ev.DF <- data.frame()
        dir_files <- list.files(directory)
        log_dirs <- file.path(directory, dir_files[grep('LOGS_', dir_files)])
        for (log_dir in log_dirs) {
            ev.DF <- plyr:::rbind.fill(ev.DF, 
                read.table(file.path(log_dir, 'history.data'), 
                    header=TRUE, skip=5))
        }
        ev.DF <- ev.DF[order(ev.DF$star_age),]
        
        make_plots(plot_HR, paste0(basename(directory), "-HR"), 
            filepath=file.path('plots', dirname(directory), 'HR'), 
            mar=c(2, 2, 0.5, 3.5), DF=DF, ev.DF=ev.DF, use.cairo=T,
            #thin=F, short=F, make_pdf=F, slides=F)
            wide=F, tall=F, slides=F)
        
        #make_plots(plot_chem_c, paste0(basename(directory), "-chem_c"), 
        #    filepath=file.path('plots', dirname(directory), 'chem_c'), 
        #    mar=c(2, 3, 0.5, 2.5), 
        #    DF=DF, ev.DF=ev.DF, use.cairo=T, 
        #    plot_legend=F, 
        #    #make_png=F, 
        #    wide=F, tall=F, slides=F)#,
        #    #thin=F, short=F, make_pdf=F, slides=F)
        
        make_plots(plot_chem_c2, paste0(basename(directory), "-chem_c2"), 
            filepath=file.path('plots', dirname(directory), 'chem_c2'), 
            mar=c(2.5, 3, 0.5, 2.5), 
            DF=DF, ev.DF=ev.DF, use.cairo=T, 
            #plot_legend=F, 
            cex.paper=0.66, 
            #make_png=F, 
            wide=F, tall=F, slides=F)#,
            #thin=F, short=F, make_pdf=F, slides=F)
        
        #make_plots(plot_chem_surf2, paste0(basename(directory), "-chem_surf2"), 
        #    filepath=file.path('plots', dirname(directory), 'chem_surf2'), 
        #    mar=c(2, 3, 0.5, 2.5), 
        #    DF=DF, ev.DF=ev.DF, use.cairo=T, 
        #    plot_legend=F, 
        #    #make_png=F, 
        #    wide=F, tall=F, slides=F)#,
        #    #thin=F, short=F, make_pdf=F, slides=F)
        
        make_plots(plot_chem_surf, paste0(basename(directory), "-chem_surf"), 
            filepath=file.path('plots', dirname(directory), 'chem_surf'), 
            mar=c(2, 3, 0.5, 2.5), 
            DF=DF, ev.DF=ev.DF, use.cairo=T,
            #plot_legend=F, 
            #make_png=F, 
            wide=F, tall=F, slides=F)#,
            #thin=F, short=F, make_pdf=F, slides=F)
        }
    }
    
    warnings()
}

