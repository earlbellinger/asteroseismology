#### HR diagram showing KAGES+LEGACY stars separated by [Fe/H] 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 
#library(emojifont)
library(parallelMap)
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

directory <- file.path('learn-final') #file.path('feh', 'learn-feh')
#directory <- file.path('learn-gaia')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    print(filename)
    
    covs <- read.table(filename, header=T)
    covs <- covs[complete.cases(covs),]
    cov.DF <- with(covs, 
        data.frame(KIC=KIC,
                   M=median(M),         e_M=mad(M),
                   alpha=median(alpha), e_alpha=mad(alpha),
                   age=median(age),     e_age=mad(age)))
    
    if ('radius' %in% names(covs))
        cov.DF <- cbind(cov.DF,
            with(covs, data.frame(R=median(radius), e_R=mad(radius))))
    
    if ('L' %in% names(covs))
        cov.DF <- cbind(cov.DF,
            with(covs, data.frame(L=median(L), e_L=mad(L))))
    
    perturb <- read.table(file.path('perturb', 'feh', #'gaia', 
        sub('\\.', '_perturb.', basename(filename))), header=1)
    perturb.DF <- with(perturb,
            data.frame(Teff=median(Teff), e_Teff=mad(Teff),
                       FeH=median(Fe.H),  e_FeH=mad(Fe.H)))
    
    if ('Dnu0' %in% names(perturb)) 
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(Dnu=median(Dnu0),  e_Dnu=mad(Dnu0))))
    
    if ('dnu02' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(dnu=median(dnu02), e_dnu=mad(dnu02))))
    
    if ('L' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(L=median(L), e_L=mad(L))))
    
    if ('radius' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(R=median(radius), e_R=mad(radius))))
    
    cbind(cov.DF, perturb.DF)
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(file.path(directory))


track.col <- "#222222"#'black'
point.col <- 'darkblue'
point.border <- 'white'#'black'

plot_HR <- function(FeH.bin, ..., 
        make.x=T, make.y=T, make.top=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.5, 1.5, 0.1), lwd=1.5, las=1, cex.axis=text.cex)
    
    xlim <- c(7000, 4800)
    ylim <- c(0.2, 20)
    
    plot(NA, axes=F, log='y', 
        xaxs='i', yaxs='i', 
        xlim=xlim, ylim=ylim, 
        xlab="", ylab="")
    
    
    #spectral.divs <- log10(c(30000, 10000, 7500, 6000, 5200, 3700, 2400))
    spectral.divs <- c(30000, 10000, 7500, 6000, 5200, 3700, 2400)
    rs <- c(175/255, 199/255, 1, 1, 1, 1, 1, 1)
    gs <- c(201, 216, 244, 229, 217, 199, 166)/255
    bs <- c(1, 1, 243/255, 207/255, 178/255, 142/255, 81/255)
    cols <- c(
        rgb(175/255, 201/255, 1),       # O
        rgb(199/255, 216/255, 1),       # B
        rgb(1,       244/255, 243/255), # A 
        rgb(1,       229/255, 207/255), # F 
        rgb(1,       217/255, 178/255), # G 
        rgb(1,       199/255, 142/255), # K 
        rgb(1,       166/255, 81/255))  # M
    #if (F) {
    for (ii in 1:length(spectral.divs)) {
        div <- spectral.divs[ii]
        if (div > xlim[1]) next 
        if (div < xlim[2]) div <- xlim[2]
        if (ii == 1) {
            rect(xlim[1], ylim[1], div, ylim[2], 
                col=cols[ii], border=NA)
        } else {
            prev <- spectral.divs[ii-1]
            if (prev > xlim[1]) prev <- xlim[1]
            rect(prev, ylim[1], div, ylim[2], 
                col=cols[ii], border=NA)
        }
    }
    for (ii in 2:(length(spectral.divs)-1)) {
        div <- spectral.divs[ii]
        gradient.rect(div+0.0025, ylim[1], div-0.0025, ylim[2],
            nslices=10, border=NA, 
            reds=c(rs[ii], rs[ii+1]), 
            greens=c(gs[ii], gs[ii+1]),
            blues=c(bs[ii], bs[ii+1]))
    }
    #}
    
    ## plot evolutionary tracks 
    
    #with(read.table('../grid/sun/sun.dat', header=1), 
    #    lines(Teff, L, lwd=1.5))
    
    track.dir <- if (FeH.bin == 1) 'feh_grid-3'
            else if (FeH.bin == 2) 'feh_grid-15'
            else if (FeH.bin == 3) 'feh_grid0'
            else if (FeH.bin == 4) 'feh_grid02'
    track.subdirs <- list.files(file.path('..', 'grid', 'feh', track.dir),
        full.names=T)
    ZAMS_Teffs <- c()
    ZAMS_Ls <- c()
    for (track.subdir in track.subdirs) {
        track.files <- list.files(track.subdir, full.names=T, recursive=T)
        histories <- track.files[grep('history.data', track.files)]
        hist.DF <- do.call(rbind, Map(function(history.file)
            read.table(history.file, header=1, skip=5), history.file=histories))
        hist.DF <- hist.DF[order(hist.DF$star_age),]
        
        decreasing_L <- with(hist.DF, 
            which((diff(log_L) < 0)# | diff(log_Teff) < 0) 
                & center_h1[-1] > 0.65))
        if (any(decreasing_L)) {
            goes_back_up <- diff(decreasing_L) > 1
            pms <- max(decreasing_L)
            hist.DF <- hist.DF[-1:-pms,]
        }
        
        #above <- 10**hist.DF$log_L > ylim[2]
        #if (any(above)) hist.DF <- hist.DF[1:(which(above)[1]),]
        
        with(hist.DF, lines(10**log_Teff, 10**log_L, lwd=1.33,#1.5, 
            col=track.col))
        
        ZAMS_Teffs <- c(ZAMS_Teffs, 10**hist.DF$log_Teff[1])
        ZAMS_Ls <- c(ZAMS_Ls, 10**hist.DF$log_L[1])
        
        #if (!hist.DF$star_mass[1] * 100 %% 2) {
        if (FeH.bin == 2 && hist.DF$star_mass[1] > 1.3 ||
            FeH.bin == 2 && hist.DF$star_mass[1] < 0.9 || 
            FeH.bin == 3 && hist.DF$star_mass[1] > 1.4 ||
            FeH.bin == 3 && hist.DF$star_mass[1] < 0.9 ||
            FeH.bin == 4 && hist.DF$star_mass[1] > 1.4) next 
        par(family="Helvetica")
        with(hist.DF, 
            text((10**log_Teff[1])*0.98, (10**log_L[1]) * 0.78, 
                labels=star_mass[1],
                pos=2, cex=0.7*text.cex, col='black'))
        par(family=font)
        #}
        
    }
    #print(ZAMS_Teffs)
    #print(ZAMS_Ls)
    
    Teff.seq <- seq(xlim[2], xlim[1], 10)
    lines(Teff.seq, splinefun(ZAMS_Teffs, ZAMS_Ls)(Teff.seq),
        lwd=1.33, lty=2, col=track.col)
    
    ## plot data 
    
    FeH <- if (FeH.bin == 1) DF$FeH <= -0.2
      else if (FeH.bin == 2) DF$FeH >  -0.2 & DF$FeH <= -0.1
      else if (FeH.bin == 3) DF$FeH >  -0.1 & DF$FeH <= 0.1
      else if (FeH.bin == 4) DF$FeH >   0.1
    
    with(DF[FeH,],
        points(Teff, L, 
            pch=21, 
            cex=1.1, 
            #bg=ifelse(Teff > 6000, blue, 
            #   ifelse(Teff > 5200, orange, 'darkred')),
            bg=adjustcolor(point.col, alpha.f=0.75),
            #'darkblue', alpha.f=0.75), 
            lwd=0.75, 
            col=point.border))
    
    if (F) {
    x.mid <- 6850
    y.mid <- 0.78
    
    segments(x.mid+median(DF[FeH,]$e_Teff), y.mid, 
             x.mid-median(DF[FeH,]$e_Teff), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(DF[FeH,]$e_L, na.rm=T), 
             x.mid, y.mid-median(DF[FeH,]$e_L, na.rm=T), 
        lwd=1, col='black')
    }
    
    ## solar symbol
    if (FeH.bin == 3) {
        points(5777, 1, pch=20, cex=1.1/2,   lwd=1.5)
        points(5777, 1, pch=1,  cex=1.1, lwd=1.5)
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]*0.1, xlim[2]-100, ylim[2]*10, col='white', border=NA)
    rect(xlim[1]*1.1, ylim[1], xlim[2]*0.9, ylim[1]*0.9, col='white', border=NA)
    par(xpd=F)
    
    par(family="Helvetica")
    text(7020, 0.265, 
        labels=if (FeH.bin == 1) expression('[Fe/H]' <= -0.2)
          else if (FeH.bin == 2) expression(-0.2~'< [Fe/H]' <= -0.1) 
          else if (FeH.bin == 3) expression(-0.1~'< [Fe/H]' <= 0.1)
          else if (FeH.bin == 4) expression('[Fe/H] >'~0.1), 
        pos=4, cex=0.85*text.cex)
    text(7020, 0.42 + ifelse(FeH.bin == 2 || FeH.bin == 4, 0.02, 0), 
        labels=if (FeH.bin == 1) expression('metal poor')
          else if (FeH.bin == 2) expression('metal deficient') 
          else if (FeH.bin == 3) expression('solar metallicity')
          else if (FeH.bin == 4) expression('metal rich'), 
        pos=4, cex=0.85*text.cex)
    par(family=font)
    
    nxticks <- 4
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    #yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    #yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    #yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp=mgp+c(0, 0.25, 0))
    #xpos <- seq(10**xlim[2], 10**xlim[1], 2000)
    #xpos2 <- seq(10**xlim[2], 10**xlim[1], 500)
    xpos <- xticks#seq(xlim[2], xlim[1], 1000)
    xpos2 <- xticks.minor#seq(xlim[2], xlim[1], 200)
    #axis(side=1, tcl=tcl/2, at=log10(xpos2), labels=F, lwd.ticks=par()$lwd)
    axis(side=1, tcl=tcl/2, at=xpos2, labels=F, lwd.ticks=par()$lwd)
    #axis(side=1, tcl=tcl, at=log10(xpos), labels=xpos, cex.axis=text.cex,
    axis(side=1, tcl=tcl, at=xpos, labels=if (make.x) xpos else F, 
        cex.axis=text.cex, lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, 0.43, 0))
    #axis(2, tcl=tcl, lwd=0, lwd.ticks=par()$lwd, tick=T, at=yticks,
    #    labels=as.logical(make.y))
    #axis(2, tcl=tcl/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
    #    at=yticks.minor, labels=F)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=make.y, lwd.ticks=par()$lwd)
    
    box(lwd=par()$lwd)
    
    if (make.y) {
        mtext(expression(L/L), 2, 2, outer=F, las=0, cex=text.cex)
        par(xpd=NA)
        points(7380, 3.15, pch=20, cex=0.75/3, lwd=1)
        points(7380, 3.15, pch=1,  cex=0.75,   lwd=1)
        par(xpd=F)
    }
    if (make.x) mtext(expression(T["eff"]/K), 1, 2, 
        outer=F, cex=text.cex)
    
    #par(mgp=mgp+c(0.55, 0, 0))
    #if (make.x) title(xlab=expression(T["eff"]/K))
    #par(mgp=mgp+c(0.55, 0, 0))
    ##if (make.y) title(ylab=expression(log[10](L/L["sun"])))
    #if (make.y) title(ylab=expression(L/L["sun"]))
    
    if (make.top) 
        mtext(expression("Spectral Type"), side=3, line=1.3, cex=text.cex)
    
    spectral.labs <- c("O", "B", "A", "F", "G", "K", "M")
    selector <- 1:(length(spectral.divs))
    spectral.Teffs <- sapply(selector, 
    function(ii) {
        div <- spectral.divs[ii]
        if (div > xlim[1]) return(Inf)
        if (div < xlim[2]) div <- xlim[2]
        if (ii == 1) return((xlim[1]+div)/2)
        prev <- spectral.divs[ii-1]
        if (prev > xlim[1]) prev <- xlim[1]
        (div+prev)/2
    })
    axis(3, at=spectral.divs, tcl=tcl, labels=F, cex.axis=text.cex,
        lwd.ticks=par()$lwd)
    par(mgp=mgp+c(0, -0.05, 0))
    if (make.top)
        axis(3, at=spectral.Teffs, labels=spectral.labs[selector], 
            cex.axis=text.cex, tcl=0)
}

for (FeH.bin in 1:4) {
print(FeH.bin)
make_plots(plot_HR, paste0('hr', FeH.bin), 
    FeH.bin=FeH.bin, 
    make.x=FeH.bin >= 3, 
    make.y=FeH.bin == 1 || FeH.bin == 3, 
    make.top=FeH.bin <= 2,
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, 
    wide=F, tall=F, 
    paper_pdf_height=4.17309*1.385,
    cex.paper=0.95,
    use.cairo=T, font='Palatino Linotype') 
}


