#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'seismology.R'))
library(parallel)
library(parallelMap)

Z_div_X_solar = 0.02293
profile.pattern <- 'profile.+.data$'
freqs.pattern <- 'profile.+-freqs.dat$'
freqs.cols <- c('l', 'n', 'nu', 'inertia')

### Obtain observable properties from models 
summarize <- function(pro_file, freqs_file, ev.DF) {
    print(pro_file)
    
    pro_header <- read.table(pro_file, header=TRUE, nrows=1, skip=1)
    hstry <- ev.DF[ev.DF$model_number==pro_header$model_number,]
    if (nrow(hstry) == 0) {#|| hstry$mass_conv_core > 0) 
        print(paste("Model", pro_file, "failed"))
        return(NULL)
    }
    
    obs.DF <- NULL
    
    ## Model properties 
    obs.DF["age"] <- pro_header$star_age/10**9
    obs.DF["mass_cc"] <- hstry$mass_conv_core/pro_header$star_mass
    obs.DF["mass_X"] <- pro_header$star_mass_h1/pro_header$star_mass
    obs.DF["mass_Y"] <- (pro_header$star_mass_he3 + 
            pro_header$star_mass_he4)/pro_header$star_mass
    obs.DF["X_c"] <- hstry$center_h1
    obs.DF["X_surf"] <- hstry$surface_h1
    obs.DF["Y_surf"] <- hstry$surface_he3 + hstry$surface_he4
    
    ## Observable properties 
    obs.DF["radius"] <- pro_header$photosphere_r
    obs.DF["L"] <- pro_header$photosphere_L
    obs.DF["log_g"] <- hstry$log_g
    obs.DF["Teff"] <- pro_header$Teff
    obs.DF["Fe/H"] <- log10(10**hstry$log_surf_cell_z / 
            hstry$surface_h1 / Z_div_X_solar)
    
    ## Seismology 
    freqs <- read.table(freqs_file, col.names=freqs.cols, fill=TRUE)
    #acoustic_cutoff <- hstry$acoustic_cutoff/(2*pi)
    nu_max <- hstry$nu_max
    seis.DF <- seismology(freqs, nu_max, #acoustic_cutoff=acoustic_cutoff, 
        outf=ifelse(sample(0:1000, 1) == 0, gsub("/", "-", freqs_file), FALSE),
        filepath=file.path('plots', 'separation'))
    
    if (all(is.na(seis.DF))) return(NULL)
    merge(rbind(obs.DF), rbind(seis.DF))
}

### Obtain evolutionary tracks from a MESA directory
parse_dir <- function(directory, min_num_models=15) {
    ## parse dirname string e.g. "M=1.0_Y=0.28"
    params.DF <- NULL
    for (var in unlist(strsplit(basename(directory), '_'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params.DF[nameval[1]] <- as.numeric(nameval[2])
    }
    
    ## obtain history
    log_dir <- file.path(directory, "LOGS")
    logs <- list.files(log_dir)
    if (length(logs) <= 1) {
        print(paste(directory, "No logs found!"))
        return(NA)
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
        return(NA)
    }
    
    ## call summarize on all pairs 
    parallelStartMulticore(max(1, detectCores()))
    obs.DF <- do.call(plyr:::rbind.fill, 
        parallelMap(function(pro_file, freqs_file)
                summarize(pro_file, freqs_file, ev.DF), 
            pro_file=file.path(log_dir, pro_files), 
            freqs_file=file.path(log_dir, freq_files)))
    
    merge(rbind(params.DF), obs.DF[with(obs.DF, order(age)),])
}

### Hertzsprung-Russell diagram
plot_HR <- function(DF, ev.DF, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp,
        col.pal=colorRampPalette(brewer.pal(11, "Spectral"))(21)) {
    plot(#log10(DF$Teff), log10(DF$L), 
        ev.DF$log_L ~ ev.DF$log_Teff,
        type='l', tcl=0, axes=FALSE,
        xlab=expression('Temperature' ~ 'lg'*(T[eff]/K)),
        ylab=expression('Luminosity' ~ 'lg'*(L / L['\u0298'])),
        xlim=rev(log10(range(DF$Teff))),
        ylim=log10(range(DF$L)))
    abline(v=log10(5771), lty=3, col='lightgray')
    abline(v=log10(7000), lty=4, col='gray')
    abline(h=0, lty=3, col='lightgray')
    points(log10(5771), 0, pch=1, cex=1)
    points(log10(5771), 0, pch=20, cex=0.1)
    points(log10(DF$Teff), log10(DF$L), pch=1, cex=0.1, 
        col=col.pal[floor((DF$X_c-min(DF$X_c)) / (max(DF$X_c)-min(DF$X_c))
                          * (length(col.pal)-1)) + 1])
    magaxis(side=1:4, family=font, tcl=0.25, mgp=utils.mgp, las=1,
        cex.axis=text.cex, labels=c(1,1,0,0)) 
    var1range <- diff(par()$usr)[1] # Add colorbar
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        round(quantile(seq(min(DF$X_c), max(DF$X_c), length=length(col.pal)), 
            c(0, 0.25, 0.5, 0.75, 1)), 3), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(expression("Core-hydrogen" ~ X[c]), 4, line=5.5, cex=text.cex)
    par(family="Luxi Mono")
    legend("bottom", bty='n', xjust=1, cex=text.cex/2, col="gray", 
        legend=c(
            as.expression(bquote(M == .(DF$M[1]))), 
            as.expression(bquote(Y == .(DF$Y[1]))), 
            as.expression(bquote(Z == .(DF$Z[1]))), 
            as.expression(bquote(alpha == .(DF$alpha[1]))),
            as.expression(bquote(f == .(DF$overshoot[1]))),
            as.expression(bquote(D == .(DF$diffusion[1])))
        )
    )
}

## Kippenhahn diagram of mass vs age showing convective regions 
hatches <- function(conv_bot, conv_top, age) {
    disconts <- c(0, which(abs(diff(conv_bot))>=0.1 | abs(diff(conv_top))>=0.1),
                  length(conv_top)+1)
    for (ii in 1:(length(disconts)-1)) {
        selection <- (disconts[ii]+1):(disconts[ii+1]-1)
        polygon(x =  c(age[selection], rev(age[selection])), 
            y = c(conv_bot[selection], rev(conv_top[selection])), 
            density=10, angle=60, border='black',
            col=ifelse(all(conv_bot[selection] < 0.1), "purple", "blue"))
    }
}

plot_Kippenhahn <- function(DF, ev.DF, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    age <- ev.DF$star_age/10**9
    plot(NA, 
        type='l', axes=FALSE, 
        xlab=expression("Age"~tau/"Gyr"), 
        ylab=expression("Mass"~m/M['\u0298']), 
        xaxs='i', yaxs='i',
        ylim=c(0, max(DF$M)), 
        xlim=c(min(age), max(DF$age)))
    magaxis(1:4, labels=c(1,1,0,0), mgp=c(2, 0.5, 0), family=font, las=1, 
        tcl=-0.25, cex.axis=text.cex)
    hatches(ev.DF$conv_mx1_bot * ev.DF$star_mass, 
            ev.DF$conv_mx1_top * ev.DF$star_mass, age)
    hatches(ev.DF$conv_mx2_bot * ev.DF$star_mass, 
            ev.DF$conv_mx2_top * ev.DF$star_mass, age)
    abline(v=min(DF$age), lty=2, col='gray')
    #lines(age, ev.DF$he_core_mass * ev.DF$star_mass, 
    #    lty=3, col="darkred")
    par(family="Luxi Mono")
    #legend("right", bty='n', cex=text.cex/2, lty=c(2, 3, 4),
    #    col=c("gray", "darkred"),
    #    legend=c("ZAMS", "Helium core"))
    legend("left", bty='n', xjust=1, cex=text.cex/2, col="gray", 
        legend=c(
            as.expression(bquote(M == .(DF$M[1]))), 
            as.expression(bquote(Y == .(DF$Y[1]))), 
            as.expression(bquote(Z == .(DF$Z[1]))), 
            as.expression(bquote(alpha == .(DF$alpha[1]))),
            as.expression(bquote(f == .(DF$overshoot[1]))),
            as.expression(bquote(D == .(DF$diffusion[1])))
        )
    )
}

## Plot frequency separations
plot_separations <- function(DF, ..., 
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    attach(DF)
    plot(NA, axes=0, xaxs='i', yaxs='i',
        xlim=range(age),
        ylim=c(0, log10(200)),
        xlab=expression("Age"~tau/"Gyr"),
        ylab=expression("Frequency"~nu/mu*Hz))
    #magaxis(1:3, labels=c(1,1,0), mgp=c(2, 0.5, 0), family=font, las=1, 
    #    tcl=-0.25, unlog='y')
    magaxis(1:2, labels=c(1,1), mgp=c(2, 0.5, 0), family=font, las=1, 
        tcl=-0.25, unlog='y')
    lines(log10(Dnu0_median) ~ age, lty=2)
    lines(log10(dnu02_median) ~ age, lty=3)
    lines(log10(dnu13_median) ~ age, lty=4)
    legend("bottomleft", lty=c(2,3,4), bty='n', 
       legend=c(expression(Delta*nu), 
                expression(delta*nu[0*","*2]),
                expression(delta*nu[1*","*3])))
    
    ## Build top x-axis
    #tick.locs <- pretty(age)
    #locations <- Map(function(x) 
    #    which(abs(age-x)==min(abs(age-x))), x=tick.locs)
    #xc_vals <- round(X_c[unlist(locations)], 3)
    #axis(3, at=tick.locs, tcl=-0.25, labels=xc_vals)
    
    tick.locs <- as.numeric(round(quantile(
        seq(min(X_c), max(X_c), length.out=1000), 
        c(1, 0.8, 0.6, 0.4, 0.2, 0)), 2))
    if (tick.locs[length(tick.locs)] < 10**-2) {
        tick.locs[length(tick.locs)] <- 10**-2
        tick.locs <- c(tick.locs, 0)
    } 
    locations <- Map(function(x) 
        which(abs(X_c-x)==min(abs(X_c-x))), x=tick.locs)
    age_vals <- round(age[unlist(locations)], 3)
    axis(3, at=c(min(age), age_vals[-1]), tcl=-0.25, 
        labels=c(tick.locs))
    
    xc_minors <- c()
    for (xc_ii in 1:(length(tick.locs)-1)) {
        xc <- tick.locs[xc_ii]
        mins <- as.numeric(quantile(
                seq(tick.locs[xc_ii], tick.locs[xc_ii+1], length.out=100),
            c(0.8, 0.6, 0.4, 0.2)))
        xc_minors <- c(xc_minors, mins)
    }
    minor.locations <- Map(function(x) 
        which(abs(X_c-x)==min(abs(X_c-x))), x=xc_minors)
    minor.locs <- age[unlist(minor.locations)]
    axis(3, at=minor.locs, labels=F, tcl=-0.125)
    
    mtext(expression("Fractional core-hydrogen abundance"~X[c]), line=1.25)
    
    par(new=T)
    plot(NA, axes=0, xaxs='i', yaxs='i', xlab='', ylab='',
        xlim=range(age),
        ylim=c(0, max(r_sep02_median, r_sep13_median, 
                      r_avg01_median, r_avg10_median, 0.2)))
    lines(r_sep02_median ~ age, lty=2, col='darkred')
    lines(r_sep13_median ~ age, lty=3, col='darkred')
    lines(r_avg01_median ~ age, lty=4, col='darkred')
    magaxis(4, labels=1, mgp=c(2, 0.5, 0), family=font, cex.axis=text.cex,
        las=1, tcl=-0.25)
    mtext(expression("Frequency ratio"~r), side=4, line=2.5)
    legend("bottomright", lty=c(2,3,4,5), col='darkred', bty='n', 
       legend=c(expression(r[0*","*2]), 
                expression(r[1*","*3]), 
                expression(r[0*","*1])))
    
    legend("bottom", bty='n', xjust=1, cex=text.cex/2, col="gray", 
        legend=c(
            as.expression(bquote(M == .(M[1]))), 
            as.expression(bquote(Y == .(Y[1]))), 
            as.expression(bquote(Z == .(Z[1]))), 
            as.expression(bquote(alpha == .(alpha[1]))),
            as.expression(bquote(f == .(overshoot[1]))),
            as.expression(bquote(D == .(diffusion[1])))
        )
    )
    detach(DF)
}

###############################################################################
### Parse command line arguments and process directory ########################
###############################################################################
args <- commandArgs(TRUE)
if (length(args)>0) {
    directory <- args[1]
    print(directory)
    DF <- unique(parse_dir(directory))
    DF <- DF[complete.cases(DF),]
    DF <- DF[order(DF$age),]
    
    print(sapply(DF, fivenum))
    
    ## check if the track should be kept 
    hs <- 1-DF$Y[1]-DF$Z[1]
    rejection <- ""
    if (nrow(DF) < 5) {
        print("Rejecting track: too few data points")
        print(DF)
        rejection <- "-r"
    } else if (ncol(DF) < 5) {
        print("Rejecting track: too few columns")
        print(DF)
        rejection <- "-r"
    #} else if ( 100*(hs-max(DF$H))/hs > 1.5 ) {
    #    print("Rejecting track: inaccurate starting H")
    #    print(c(hs, max(DF$X_c)))
    #    rejection <- "-pms"
    } else { # Make a table of results! 
        write.table(DF, paste0(directory, '.dat'), quote=FALSE, sep='\t', 
            row.names=FALSE)
    }
    
    ## make plots
    ev.DF <- read.table(file.path(directory, 'LOGS', 'history.data'), 
            header=TRUE, skip=5)
    make_plots(plot_HR, paste0(basename(directory), "-HR", rejection), 
        filepath=file.path('plots', dirname(directory), 'HR'), 
        mar=c(3, 4, 1, 7), DF=DF, ev.DF=ev.DF)
    make_plots(plot_Kippenhahn, 
        paste0(basename(directory), "-Kippenhahn", rejection), 
        filepath=file.path('plots', dirname(directory), 'Kippenhahn'), 
        DF=DF, ev.DF=ev.DF)
    make_plots(plot_separations, 
        paste0(basename(directory), "-separations", rejection), 
        filepath=file.path('plots', dirname(directory), 'separations'), 
        DF=DF, mar=c(3, 4, 3, 4))
}

