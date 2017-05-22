#### Animate asteroseismological quantities across an evolutionary track
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Department of Astronomy, Yale University;
#### Stellar Ages & Galactic Evolution Group,
#### Max-Planck-Institut fur Sonnensystemforschung

library(parallel)
library(parallelMap)
library(data.table)
source(file.path('/', 'scratch', 'seismo', 'bellinger',
    'asteroseismology', 'scripts', 'seismology.R'))

ells <- c(10, 11, 9, 3)

log_dir <- 'LOGS'#'LOGS_SG'
ev.DF <- rbind(
    #read.table(file.path('LOGS_MS', 'history.data'), header=TRUE, skip=5),
    read.table(file.path(log_dir, 'history.data'), header=TRUE, skip=5)
)
logs <- list.files(log_dir)

freq_files <- logs[grep('profile.+-freqs.dat$', logs)]
prof_files <- sub('-freqs.dat', '.data', freq_files)
#mdl_nums <- T#as.numeric(sub('.data', '', sub('^[^0-9]+', '', prof_files)))

parallelStartMulticore(8)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
mdl_nums <- as.numeric(parallelMap(function(pro_file)
        read.table(pro_file, header=TRUE, nrows=1, skip=1)$model_number,
    pro_file=file.path(log_dir, prof_files)))

low.mdl <- T#mdl_nums < 2000+min(ev.DF$model_number)

get_freqs <- function(model_number) {
    idx <- which(mdl_nums == model_number)
    freq_file <- freq_files[idx]
    #prof_file <- prof_files[idx]
    #parse_freqs(file.path(freqs_dir, freq_file), gyre=1)
    #read.table(file.path(log_dir, freq_file), 
    #    col.names=c('l', 'n', 'nu', 'E'))
    parse_freqs(file.path(log_dir, freq_file))#, gyre=T)
}

get_closest <- function(freqs, nu_max, ells) {
    new.freqs <- data.frame()
    for (ell in 1:length(ells)) {
        freqs.ell <- freqs[freqs$l==ell-1,]
        for (ii in 1:ells[ell]) {
            closest <- find_closest2(freqs.ell$nu, nu_max)
            new.freqs <- plyr:::rbind.fill(new.freqs, freqs.ell[closest,])
            freqs.ell <- freqs.ell[-closest,]
        }
    }
    new.freqs
}

collapse_freqs <- function(model_number) {
    #print(model_number)
    DF <- get_freqs(model_number)
    if (nrow(DF) <= 0) return(data.frame(NA))
    #DF <- DF[DF$l == 0,]
    nu_max <- ev.DF[ev.DF$model_number==model_number,]$nu_max
    DF <- get_closest(freqs=DF, nu_max=nu_max, ells=ells)
    freqs <- rbind(DF$nu)
    colnames(freqs) <- paste0('l', DF$l, '.', 'n', DF$n)
    as.data.frame(freqs)
}

track <- ev.DF#[ev.DF$model_number %in% mdl_nums[low.mdl],]

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
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freqs <- get_closest(freqs=freqs, nu_max=nu_max, ells=ells)
    large_sep <- avg(Dnu, freqs)[[1]]
    
    freq_range <- range(freqs$nu) #c(nu_max-5*large_sep, nu_max+5*large_sep)
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
        plot(log_Teff, log_L, axes=F, tcl=0, type='l', lwd=2, 
            xlim=rev(range(log_Teff)), 
            xlab=expression("Temperature" ~ log(T["eff"]/K)),
            ylab=expression("Luminosity" ~ log(L/L["sun"]))))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
        labels=c("O→", "B→", "A→", "F→", "G→", "K→", "M→"))
    with(model, points(log_Teff, log_L))
    legend("bottomright", bty='n', #inset=c(-0.1, -0.05),
        legend=c(paste("Age:", signif(age, 4), "/", 
                       signif(max(ages), 4), "Gyr"),
                 paste("Model:", 
                       model$model_number-min(track$model_number)+1, "/",
                       max(track$model_number)-min(track$model_number)+1)))
    
    
    # u kernel
    lowest.l1 <- freqs[freqs$l==1,][1,]
    u.Y <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]),
      'E_K_u-Y.dat'), header=1)
    K <- u.Y[[paste0('l.', lowest.l1$l, '_n.', lowest.l1$n)]]
    x <- u.Y$x
    if (is.null(K)) K <- x
    plot(NA, axes=F, tcl=0, #log='x',
        xlab=expression("Radius" ~ r/R["star"]),
        ylab=expression("Kernel" ~ K^(u*","~Y)),
        xlim=c(10**-3, 1),
        ylim=range(K[x<0.99]))
    abline(h=0, col='darkgray', lwd=2, lty=2)
    lines(x[-1], K[-1], col='darkred', lwd=2)
    legend("bottom", bty='n', 
        legend=paste0("l=", lowest.l1$l, " n=", lowest.l1$n))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    
    # mode inertia
    ##l0 <- gyre[gyre$l == 0,]
    ##l1 <- gyre[gyre$l == 1,]
    ##l2 <- gyre[gyre$l == 2,]
    #l0 <- freqs[freqs$l == 0,]
    #l1 <- freqs[freqs$l == 1,]
    #l2 <- freqs[freqs$l == 2,]
    #plot(NA, axes=FALSE, tck=0, log='y',
    #     xlim=range(l0$nu),#freq_range,
    #     #ylim=range(l0$E_norm, l1$E_norm), log='y',
    #     ylim=range(l0$E),#, l1$E), 
    #     xlab=expression("Frequency"~nu/mu*Hz),
    #     ylab=expression("Mode inertia"))
    #magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    #lines(l0$nu, l0$E)
    #lines(l1$nu, l1$E, lty=2, col=red)
    #lines(l2$nu, l2$E, lty=3, col=blue)
    #abline(v=freq_range[1], lty=3, col='darkgray')
    #abline(v=freq_range[2], lty=3, col='darkgray')
    ##lines(l0$Re.freq., l0$E_norm)
    ##lines(l1$Re.freq., l1$E_norm, lty=2, col=red)
    
    
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
    for (col_i in 1:ncol(dipoles)) lines(ages, dipoles[,col_i], lty=1, lwd=2)
    for (col_i in 1:ncol(radials)) lines(ages, radials[,col_i], lty=2, lwd=2,
        col=red)
    for (row_i in 1:ncol(dipoles)) 
        points(age, dipoles[selector, row_i], cex=0.5)
    for (row_i in 1:ncol(dipoles)) 
        points(age, radials[selector, row_i], cex=0.5, col=red)
    legend("bottomleft", lty=c(2,1), pch=20, cex=0.5, col=c(red, 'black'), 
           inset=c(0.02, 0.02), legend=c(paste0("\u2113=", 0:1)))
    
    
    
    # propagation
    ell <- 1
    S <- DF$csound**2 / (DF$radius * solar_radius)**2
    lambs <- 10**6/(2*pi) * sqrt(ell*(ell+1) * S)
    
    brunt_N2 <- DF$brunt_N2
    stable <- which(brunt_N2 >= 0)
    indices <- if (any(diff(stable)!=1)) {
        stable[which(diff(stable)!=1)[1]+1] : max(stable)
    } else {
        stable
    }
    brunt <- 10**6/(2*pi) * sqrt(brunt_N2[indices])
    
    plot(NA, axes=F, tcl=0, type='l', 
        log='y', #xaxs='i', #cex=0.5,
        #ylim=freq_range,
        ylim=range(10, 10**5),
             #c(max(1, nu_max-20*large_sep), 
             #  min(10**5, nu_max+20*large_sep)),
        #range(1, 10**5),#range(lambs),#, 0.01, 10**6), 
        #ylim=c(1, 10**5),
        #xlim=range(DF$radius/max(DF$radius)),#range(0.01, max(DF$radius)), 
        xlim=c(10**-3, 1),
        xlab=expression("Radius" ~ r/R["star"]),
        ylab=expression("Frequency" ~ nu/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    abline(h=freq_range[1], lty=3, lwd=1)
    abline(h=freq_range[2], lty=3, lwd=1)
    #for (ell in 1:3) lines(DF$radius/max(DF$radius), lambs[[ell]], lty=ell+1)
    boundaries <- which(diff(DF$gradr - DF$grada < 0) != 0)
    abline(v=DF$radius[boundaries]/max(DF$radius), lty=3, lwd=2, col='darkgray')
    #lines(DF$radius[indices]/max(DF$radius), brunt)
    lines(DF$radius/max(DF$radius), lambs, lty=2, lwd=2, col=blue)
    #lines((DF$radius/max(DF$radius))[indices], brunt, lwd=2)
    lines(DF$radius/max(DF$radius), 10**6/(2*pi) * sqrt(brunt_N2), lwd=2)
    legend("top", lty=c(1,2), lwd=2, col=c('black', blue), #pch=c(1, NA),
        legend=c("Brunt", "Lamb"))
    
    
}

#model_number <- mdl_nums[1]
#options(bitmapType='cairo')
#png('Rplots.png', width=800, height=600)
#hr_echelle(model_number)
#dev.off()

dir.create(file.path('animate'), showWarnings = FALSE)
options(bitmapType='cairo')
#for (model_number in sort(mdl_nums[low.mdl])) {
#    png(file.path('animate', 
#            paste0(formatC(model_number, width=6, format='d', flag=0), '.png')),
#        width=800, height=600)
#    hr_echelle(model_number)
#    dev.off()
#}

parallelMap(function(model_number) {
    png(file.path('animate', 
            paste0(formatC(model_number, width=6, format='d', flag=0), '.png')),
        width=800, height=600)
    hr_echelle(model_number)
    dev.off()
}, model_number=sort(mdl_nums[low.mdl]))


