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

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    print(filename)
    
    covs <- read.table(filename, header=T)
    covs <- covs[complete.cases(covs),]
    cov.DF <- with(covs, 
        data.frame(KIC=KIC,
                   R=mean(radius),    e_R=sd(radius),
                   M=mean(M),         e_M=sd(M),
                   alpha=mean(alpha), e_alpha=sd(alpha),
                   age=mean(age),     e_age=sd(age)))
    
    if ('L' %in% names(covs))
        cov.DF <- cbind(cov.DF,
            with(covs, data.frame(L=mean(L), e_L=sd(L))))
    
    perturb <- read.table(file.path('perturb', 'feh', 
        sub('\\.', '_perturb.', basename(filename))), header=1)
    perturb.DF <- with(perturb,
            data.frame(Teff=mean(Teff), e_Teff=sd(Teff),
                       FeH=mean(Fe.H),  e_FeH=sd(Fe.H)))
    
    if ('Dnu0' %in% names(perturb)) 
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(Dnu=mean(Dnu0),  e_Dnu=sd(Dnu0))))
    
    if ('dnu02' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(dnu=mean(dnu02), e_dnu=sd(dnu02))))
    
    if ('L' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(L=mean(L), e_L=sd(L))))
    
    if ('Teff' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(Teff=mean(Teff), e_Teff=sd(Teff))))
    
    cbind(cov.DF, perturb.DF)
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(file.path(directory))


track.col <- 'black'
point.col <- blue #'darkblue'
point.border <- 'white'#'black'

plot_CD <- function(FeH.bin, ..., 
        make.x=T, make.y=T, make.top=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.5, 1.5, 0.1), lwd=1.5, las=1, cex.axis=text.cex)
    
    xlim <- c(0, 280)
    ylim <- c(-4, 22)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlim=xlim, ylim=ylim, 
        xlab="", ylab="")
    
    abline(h=0, lwd=1.5, lty=3, col='lightgray')
    
    ## plot evolutionary tracks 
    track.dir <- if (FeH.bin == 1) 'feh_grid2_-3'
            else if (FeH.bin == 2) 'feh_grid2_-15'
            else if (FeH.bin == 3) 'feh_grid2_0'
            else if (FeH.bin == 4) 'feh_grid2_02'
    
    track.files <- list.files(file.path('..', 'grid', #'feh2', 
            track.dir),
        full.names=T)
    track.files <- track.files[grepl('.dat', track.files)]
    
    track.DFs <- Map(function(filename) {
        track.DF <- read.table(filename, header=1)
        decreasing_L <- with(track.DF, 
            which((diff(log_L) < 0 | 
                   diff(log_Teff) < 0) & 
                   X_c[-1] > 0.65))
        if (any(decreasing_L)) {
            goes_back_up <- diff(decreasing_L) > 1
            pms <- max(decreasing_L)
            track.DF <- track.DF[-1:-pms,]
        }
        track.DF
    }, track.files)
    
    Dnu.list <- list() 
    dnu.list <- list() 
    X_cs <- c(0.01, seq(0.1, 0.7, 0.1))
    for (track_i in 1:length(track.DFs)) {
        track.DF <- track.DFs[[track_i]]
        track.DF <- track.DF[order(track.DF$age),]
        trackM <- track.DF$M[1]
        if (trackM < 0.7 || trackM > 1.7) next 
        if (!'Dnu0' %in% names(track.DF) || !'dnu02' %in% names(track.DF)) {
            print(paste0('Skipping', trackM))
            next
        }
        track.DF <- track.DF[complete.cases(track.DF$Dnu0),]
        track.DF <- track.DF[complete.cases(track.DF$dnu02),]
        if (nrow(track.DF) <= 5) {
            print(paste0('Skipping', trackM))
            next
        }
        
        track.DF <- track.DF[track.DF$X_c >= 0.005,] #& track.DF$X_c <= 0.8,]
        Dnus <- seq(min(track.DF$Dnu0), max(track.DF$Dnu0), 0.1)
        Dnu.spl <- smooth.spline(track.DF$Dnu0, track.DF$dnu02, spar=0.5)
        dnus <- predict(Dnu.spl, Dnus)$y
        lines(Dnus, dnus, lwd=1, col="#222222")
            #lines(Dnu0, dnu02, lwd=1.5, col='darkgray'))
        
        #if (FeH.bin == 2 && hist.DF$star_mass[1] > 1.2 || 
        #    FeH.bin == 3 && hist.DF$star_mass[1] > 1.3) next 
        par(family="Helvetica")
        exclude <- c(1.2, 1.4, 1.5, 1.6, 0.7, 1.7)
        text(max(Dnus)-1, max(dnus)-0.95 + if(trackM > 1.1) 0.2 else 0, 
            labels=if (!trackM %in% exclude) trackM else '',
            pos=3, cex=0.6*text.cex, col='black')
        if (trackM == 1.7)
            text(max(Dnus)-12, max(dnus)-0.75, 
                labels=expression(M == 1.7),
                pos=3, cex=0.6*text.cex, col='black')
        par(family=font)
        
        Dnu_Xc_spl <- smooth.spline(track.DF$X_c, track.DF$Dnu0,  spar=0.5)
        dnu_Xc_spl <- smooth.spline(track.DF$X_c, track.DF$dnu02, spar=0.5)
        for (X_c_ii in 1:length(X_cs)) {
            X_c <- X_cs[X_c_ii]
            #if (X_c > max(track.DF$X_c)) X_c <- max(track.DF$X_c)
            if (X_c == 0.7) X_c <- max(track.DF$X_c)
            if (X_c == 0.005) X_c <- min(track.DF$X_c)
            new_Dnu <- predict(Dnu_Xc_spl, X_c)$y
            new_dnu <- predict(dnu_Xc_spl, X_c)$y
            Dnu.list[[X_c_ii]] <- if (length(Dnu.list)>=X_c_ii) 
                c(Dnu.list[[X_c_ii]], new_Dnu) else new_Dnu
            dnu.list[[X_c_ii]] <- if (length(dnu.list)>=X_c_ii) 
                c(dnu.list[[X_c_ii]], new_dnu) else new_dnu
        }
    }
    
    for (X_c_ii in 1:length(X_cs)) {
        lines(Dnu.list[[X_c_ii]], dnu.list[[X_c_ii]], 
            lwd=1, lty=2, col="#222222")
    }
    for (X_c_ii in 1:length(X_cs)) { 
        X_c <- X_cs[X_c_ii]
        par(family="Helvetica")
        if (X_c < 0.6 && X_c > 0.01) text(max(Dnu.list[[X_c_ii]]) - 7, 
            dnu.list[[X_c_ii]][which.max(Dnu.list[[X_c_ii]])] - 0.2, 
            labels=X_c,
            pos=4, cex=0.6*text.cex, col='black')
        if (X_c == 0.6) text(max(Dnu.list[[X_c_ii]]) - 5, 
            dnu.list[[X_c_ii]][which.max(Dnu.list[[X_c_ii]])] - 0.2, 
            labels=bquote(X['c'] == .(X_c)),
            pos=4, cex=0.6*text.cex, col='black')
        if (X_c == 0.7) text(max(Dnu.list[[X_c_ii]]) - 6.5, 
            dnu.list[[X_c_ii]][which.max(Dnu.list[[X_c_ii]])] + 0.6, 
            labels='ZAMS',
            pos=4, cex=0.6*text.cex, col='black')
        if (X_c == 0.01) 
            text(Dnu.list[[X_c_ii]][which.min(dnu.list[[X_c_ii]])] - 7, 
                min(dnu.list[[X_c_ii]]) + 0.4 + ifelse(FeH.bin == 1, 0.2,
                    ifelse(FeH.bin == 2, -0.4, 0)), 
                labels=X_c,
                pos=1, cex=0.6*text.cex, col='black')
        par(family=font)
    }
    
    ## plot data 
    FeH <- if (FeH.bin == 1) DF$FeH <= -0.2
      else if (FeH.bin == 2) DF$FeH >  -0.2 & DF$FeH <= -0.1
      else if (FeH.bin == 3) DF$FeH >  -0.1 & DF$FeH <= 0.1
      else if (FeH.bin == 4) DF$FeH >   0.1
    with(DF[FeH,], 
        points(Dnu, dnu, 
            pch=21, 
            cex=1.1, 
            #bg=adjustcolor(point.col, alpha.f=0.75), 
            bg=adjustcolor(
                ifelse(Teff > 6000, blue, 
                ifelse(Teff > 5200, orange, 'darkred')), alpha.f=0.75),
            lwd=0.75, 
            col=point.border))
    
    x.mid <- 270
    y.mid <- -2.5
    
    segments(x.mid+median(DF[FeH,]$e_Dnu), y.mid, 
             x.mid-median(DF[FeH,]$e_Dnu), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(DF[FeH,]$e_dnu, na.rm=T), 
             x.mid, y.mid-median(DF[FeH,]$e_dnu, na.rm=T), 
        lwd=1, col='black')
    
    
    ## solar symbol
    if (FeH.bin == 3) {
        #points(134.8693, 8.957, pch=20, cex=0.75,   lwd=1.5)
        #points(134.8693, 8.957, pch=1,  cex=1.5, lwd=1.5)
        points(136.5, 8.957, pch=20, cex=1.1/2,   lwd=1.5)
        points(136.5, 8.957, pch=1,  cex=1.1, lwd=1.5)
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]*0.1, xlim[2]*1.1, ylim[2]*10, col='white', border=NA)
    rect(xlim[1]*0.9, ylim[1], xlim[2]*1.1, ylim[1]*0.9, col='white', border=NA)
    par(xpd=F)
    
    par(family="Helvetica")
    text(0, 19.5, 
        labels=if (FeH.bin == 1) expression('[Fe/H]' <= -0.2)
          else if (FeH.bin == 2) expression(-0.2~'< [Fe/H]' <= -0.1) 
          else if (FeH.bin == 3) expression(-0.1~'< [Fe/H]' <= 0.1)
          else if (FeH.bin == 4) expression('[Fe/H] >'~0.1), 
        pos=4, cex=0.85*text.cex)
    text(0, 17 + ifelse(FeH.bin == 2 || FeH.bin == 4, 0.25, 0), 
        labels=if (FeH.bin == 1) expression('metal poor')
          else if (FeH.bin == 2) expression('metal deficient') 
          else if (FeH.bin == 3) expression('solar metallicity')
          else if (FeH.bin == 4) expression('metal rich'), 
        pos=4, cex=0.85*text.cex)
    par(family=font)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=2, labels=make.x, lwd.ticks=par()$lwd,)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=4, labels=make.y, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    if (make.y) 
        mtext(expression(delta*nu/mu*Hz), 2, 2, outer=F, las=0, cex=text.cex)
    if (make.x) 
        mtext(expression(Delta*nu/mu*Hz), 1, 2, outer=F, cex=text.cex)
}

for (FeH.bin in 1:4) {
print(FeH.bin)
make_plots(plot_CD, paste0('cd', FeH.bin), 
    FeH.bin=FeH.bin, 
    make.x=FeH.bin >= 3, 
    make.y=FeH.bin == 1 || FeH.bin == 3, 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, 
    wide=F, tall=F, 
    paper_pdf_height=4.17309*1.385,
    cex.paper=0.95,
    use.cairo=T, font='Palatino Linotype') 
}
