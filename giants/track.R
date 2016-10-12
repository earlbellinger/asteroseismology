library(parallel)
library(parallelMap)

source(file.path('..', 'scripts', 'seismology.R'))

log_dir <- file.path('simulations', 
    'M=1_Y=0.268_Z=0.0198_alpha=1.86_overshoot=0.05_diffusion=1', 'LOGS')
freqs_dir <- file.path('simulations', 
    'M=1_Y=0.268_Z=0.0198_alpha=1.86_overshoot=0.05_diffusion=1', 'LOGS')

ev.DF <- read.table(file.path(log_dir, 'history.data'), 
            header=TRUE, skip=5)

logs <- list.files(log_dir)

freq_files <- logs[grep('profile.+-freqs.dat$', logs)]
prof_files <- sub('-freqs.dat', '.data', freq_files)
#mdl_nums <- as.numeric(sub('.data', '', sub('^[^0-9]+', '', prof_files)))

parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
mdl_nums <- as.numeric(parallelMap(function(pro_file)
        read.table(pro_file, header=TRUE, nrows=1, skip=1)$model_number,
    pro_file=file.path(log_dir, prof_files)))

low.mdl <- mdl_nums < 2000+min(ev.DF$model_number)
#low.mdl <- mdl_nums %in% mdl_nums[
#    find_closest2(ev.DF$star_age, num_points=256)]

#with(ev.DF[ev.DF$model_number %in% mdl_nums[low.mdl],], 
#    magplot(log_Teff, log_L, xlim=rev(range(log_Teff)), type='l'))

get_freqs <- function(model_number) {
    idx <- which(mdl_nums == model_number)
    freq_file <- freq_files[idx]
    #prof_file <- prof_files[idx]
    #parse_freqs(file.path(freqs_dir, freq_file), gyre=1)
    read.table(file.path(freqs_dir, freq_file), 
        col.names=c('l', 'n', 'nu', 'E'))
}

collapse_freqs <- function(model_number) {
    print(model_number)
    DF <- get_freqs(model_number)
    if (nrow(DF) <= 0) return(data.frame(NA))
    #DF <- DF[DF$l == 0,]
    freqs <- rbind(DF$nu)
    colnames(freqs) <- paste0('l', DF$l, '.', 'n', DF$n)
    as.data.frame(freqs)
}

track <- ev.DF[ev.DF$model_number %in% mdl_nums[low.mdl],]

track_freqs <- do.call(plyr:::rbind.fill, 
    parallelMap(collapse_freqs, sort(mdl_nums[low.mdl])))


hr_echelle <- function(model_number) {
    par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,5,1,1), cex=1.5, pch=20)
    
    selector <- track$model_number == model_number
    model <- track[selector,]
    
    idx <- which(mdl_nums == model_number)
    DF <- read.table(file.path(log_dir, prof_files[idx]), 
        skip=5, header=1)
    radius <- read.table(file.path(log_dir, prof_files[idx]), 
        skip=1, nrow=1, header=1)$photosphere_r
    freqs <- get_freqs(model_number)
    large_sep <- avg(Dnu, freqs)[[1]]
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freq_range <- c(nu_max-5*large_sep, nu_max+5*large_sep)
    #freqs <- freqs[freqs$nu >= freq_range[1] & freqs$nu <= freq_range[2],]
    
    ##gyre <- read.table(file.path('simulations', 
    ##    'M=1_Y=0.266_Z=0.018_alpha=1.81_overshoot=0.07_diffusion=1', 'LOGS4',
    ##    freq_files[idx]), 
    ##  skip=5, header=1)
    #gyre <- get_freqs(model_number)
    ##gyre <- gyre[gyre$Re.freq. >= freq_range[1] & 
    ##             gyre$Re.freq. <= freq_range[2],]
    #gyre <- gyre[gyre$nu >= freq_range[1] &
    #             gyre$nu <= freq_range[2],]
    
    ages <- track$star_age/10**9
    age <- ages[selector]
    
    ## hr
    with(ev.DF, 
        plot(log_Teff, log_L, axes=F, tcl=0, type='l', 
            xlim=rev(range(log_Teff)), 
            xlab=expression("Temperature" ~ log(T["eff"]/K)),
            ylab=expression("Luminosity" ~ log(L/L["☉"]))))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    with(model, points(log_Teff, log_L))
    legend("topleft", bty='n', inset=c(-0.1, -0.05),
        legend=c(paste("Age:", signif(age, 4), "/", 
                       signif(max(ages), 4), "Gyr"),
                 paste("Model:", 
                       model$model_number-min(track$model_number)+1, "/",
                       max(track$model_number)-min(track$model_number)+1)))
    
    # propagation
    ell <- 1
    S <- DF$csound**2 / (DF$radius * solar_radius)**2
    lambs <- 10**6/(2*pi) * sqrt(ell*(ell+1) * S)
    
    #brunt_N2 <- DF$brunt_N2
    #stable <- which(brunt_N2 >= 0)
    #indices <- if (any(diff(stable)!=1)) {
    #    stable[which(diff(stable)!=1)[1]+1] : max(stable)
    #} else {
    #    stable
    #}
    #brunt <- 10**6/(2*pi) * sqrt(brunt_N2[indices])
    
    plot(NA, axes=F, tcl=0, type='l', 
        log='xy', #xaxs='i', #cex=0.5,
        #ylim=freq_range,
        ylim=range(lambs),#, 0.01, 10**6), 
        #ylim=c(1, 10**5),
        xlim=range(DF$radius/max(DF$radius)),#range(0.01, max(DF$radius)), 
        xlab=expression("Radius" ~ r/R["*"]),
        ylab=expression("Frequency" ~ nu/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    abline(h=freq_range[1], lty=3)
    abline(h=freq_range[2], lty=3)
    #for (ell in 1:3) lines(DF$radius/max(DF$radius), lambs[[ell]], lty=ell+1)
    boundaries <- which(diff(DF$gradr - DF$grada < 0) != 0)
    abline(v=DF$radius[boundaries]/max(DF$radius), lty=3, col='darkgray')
    #lines(DF$radius[indices]/max(DF$radius), brunt)
    lines(DF$radius/max(DF$radius), lambs, lty=2, col=blue)
    ##legend("topright", lty=c(1,2), #pch=c(1, NA),
    ##    legend=c("Brunt", "Lamb"), bty='n')
    
    ## freq_evolv_plot
    radials <- track_freqs[grepl('l0', colnames(track_freqs))]
    dipoles <- track_freqs[grepl('l1', colnames(track_freqs))]
    plot(NA, axes=F, xaxs='i',
        xlim=c(max(min(ages), age-.15),
               min(max(ages), age+.15)), 
        #ylim=range(radials[selector,], na.rm=1), 
        ylim=freq_range,
        xlab=expression("Star age"~tau/Gyr),
        ylab=expression("Frequency"~nu/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    for (col_i in 1:ncol(dipoles)) lines(ages, dipoles[,col_i], lty=1)
    for (col_i in 1:ncol(radials)) lines(ages, radials[,col_i], lty=2, col=red)
    for (row_i in 1:ncol(dipoles)) 
        points(age, dipoles[selector, row_i], cex=0.5)
    for (row_i in 1:ncol(dipoles)) 
        points(age, radials[selector, row_i], cex=0.5, col=red)
    legend("bottomleft", lty=c(2,1), pch=20, cex=0.5, col=c(red, 'black'), 
           inset=c(0.02, 0.02), legend=c(paste0("\u2113=", 0:1)))
    
    # mode inertia
    #l0 <- gyre[gyre$l == 0,]
    #l1 <- gyre[gyre$l == 1,]
    #l2 <- gyre[gyre$l == 2,]
    l0 <- freqs[freqs$l == 0,]
    l1 <- freqs[freqs$l == 1,]
    l2 <- freqs[freqs$l == 2,]
    plot(NA, axes=FALSE, tck=0, log='y',
         xlim=range(l0$nu),#freq_range,
         #ylim=range(l0$E_norm, l1$E_norm), log='y',
         ylim=range(l0$E),#, l1$E), 
         xlab=expression("Frequency"~nu/mu*Hz),
         ylab=expression("Mode inertia"))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    lines(l0$nu, l0$E)
    lines(l1$nu, l1$E, lty=2, col=red)
    lines(l2$nu, l2$E, lty=3, col=blue)
    abline(v=freq_range[1], lty=3, col='darkgray')
    abline(v=freq_range[2], lty=3, col='darkgray')
    #lines(l0$Re.freq., l0$E_norm)
    #lines(l1$Re.freq., l1$E_norm, lty=2, col=red)
}

png('Rplots.png', width=800, height=600)
hr_echelle(model_number)
dev.off()


for (model_number in sort(mdl_nums[low.mdl])) {
    png(file.path('plots', 
            paste0(formatC(model_number, width=6, format='d', flag=0), '.png')),
        width=800, height=600)
    hr_echelle(model_number)
    dev.off()
}

