

### LIBRARIES 
source('/scratch/seismo/bellinger/asteroseismology/scripts/utils.R') 
options(scipen=10000)

logs_dir <- 'LOGS_3MS'
prof.idx <- read.table(file.path(logs_dir, 'profiles.index'), skip=1, 
    col.names=c('mdl_num', 'priority', 'prof_num'))
hstry <- read.table(file.path(logs_dir, 'history.data'), header=1, skip=5)



taus <- c()
ages <- c()
logqs <- c()
for (prof_num in prof.idx$prof_num) { 

    prof.path <- file.path(logs_dir, paste0('profile', prof_num, '.data'))
    if (!file.exists(prof.path)) next 
    print(prof.path)
    DF <- read.table(prof.path, header=1, skip=5)
    hstry. <- hstry[hstry$model_number == prof.idx[prof.idx$prof_num == prof_num,]$mdl_num,]
    
    hif.idx <- min(which(DF$neutral_fraction_H <= 0.5))
    
    logq <- log10(1-DF$mass[hif.idx] / hstry.$star_mass)
    tau <- DF$opacity[hif.idx]
    
    ages <- c(ages, hstry.$star_age / 10**9)
    logqs <- c(logqs, logq)
    taus <- c(taus, tau)
    
}



write.table(data.frame(ages, taus), file="taus.dat")

plot_kippenhahn <- function(..., make.x=T, make.y=T, 
        text.cex=1, mgp=utils.mgp, font=utils.font, mar=utils.mar, short=F) {
    
    par(mar=mar+c(0.3, -0.5, 2, -0.1), lwd=1.66, las=1, cex.axis=text.cex)
    
    xlim <- c(0, max(hstry$star_age)/10**9)
    ylim <- c(0.1, 10000)
    
    plot(NA, 
        axes=F, xaxs='i', yaxs='i', 
        type='l', lwd=3, col=1, 
        xlim=xlim, 
        ylim=ylim, 
        xlab="", 
        ylab="",
        log='y')
    
    abline(h=2/3, lwd=par()$lwd, lty=2)
    lines(ages, taus, lwd=3)
    text(diff(xlim)/2, 2/3,
        expression(Photosphere~~(tau==2/3)),
        cex=0.8*text.cex, pos=1, 
        family='Helvetica')
    
    nxticks <- 6
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n=nxticks)
    yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n=nxticks*nxminor)
    yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    
    par(mgp=mgp+c(0, 0.3, 0))
    axis(1, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, tick=T, at=xticks,
        labels=as.logical(make.x))
    axis(1, tcl=-0.346/2, lwd=0, lwd.ticks=par()$lwd, tick=T, 
        at=xticks.minor, labels=F)
    
    magaxis(2, tcl=-0.346, lwd=0, lwd.ticks=par()$lwd, ticks=T, labels=make.y)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make.y) title(ylab=expression(Optical~depth~tau))
    if (make.x) title(xlab=expression(Star~age~Tau/Gyr)) 
    
}

make_plots(plot_kippenhahn, 'kip',
        filepath=file.path('plots', 'kip'),
        cex.paper=0.93, 
        wide=T, slides=F, #make_png=T, make_pdf=F,
        make.x=T,
        make.y=T,
        font="Palatino Linotype", 
        use.cairo=T)

#plot_kippenhahn() 


