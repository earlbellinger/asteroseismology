#### HR diagram showing a solar evolution track 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

#options(scipen=100000)
source('../scripts/utils.R') 
#library(emojifont)
library(parallelMap)

DF <- read.table('SG_US_step_sun/20000.dat', header=1)
DF <- DF[-1,]
#DF.pts <- DF[!apply(is.na(DF[,grep('r01_', names(DF))]), 1, all),]
#DF.pts <- DF[!apply(is.na(DF[,grep('dnu02', names(DF))]), 1, all),]
DF.pts <- DF[!is.na(DF[,grep('dnu02', names(DF))]),]

track.col <- "#222222"#'black'
point.col <- 'darkblue'
point.border <- 'white'#'black'

plot_HR <- function(..., 
        make.x=T, make.y=T, make.top=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.5, -0.2, -0.2), lwd=1.5, las=1, cex.axis=text.cex)
    
    xlim <- c(5840, 5560)#rev(range(DF$Teff)) #c(7000, 4800)
    ylim <- c(0.6, 2.2)#range(DF$L) #c(0.2, 20)
    
    plot(NA, axes=F, #log='y', 
        xaxs='i', yaxs='i', 
        xlim=xlim, ylim=ylim, 
        xlab="", ylab="")
    
    ## plot evolutionary track 
    arrows(DF[nrow(DF.pts),]$Teff, DF[nrow(DF.pts),]$L,
           DF[nrow(DF.pts)+3,]$Teff, DF[nrow(DF.pts)+3,]$L, 
           lwd=1.2, length=0.05, col=track.col)
    with(DF.pts, lines(Teff, L, lwd=1.33, col=track.col))
    points(DF[1,]$Teff, DF[1,]$L, pch=20, cex=0.2, col=track.col)
    
    with(DF.pts,
        points(Teff, L, 
            pch=21, 
            cex=1.1, 
            bg=adjustcolor(ifelse(ev_stage == 1, orange, 
                    ifelse(ev_stage == 2, blue, red)),
                alpha.f=0.4),
            lwd=0.75, 
            col=point.border))
    
    ## solar symbol
    points(5772, 1, pch=20, cex=1.1/2, lwd=1.5)
    points(5772, 1, pch=1,  cex=1.1,   lwd=1.5)
    
    
    par(family="Helvetica")
    with(DF, 
        text(5750, 1.08, labels='Main Sequence',
            pos=4, cex=0.8*text.cex, col='black'))
    with(DF, 
        text(5750, 0.96, labels='~ 7 Gyr',
            pos=4, cex=0.8*text.cex, col='black'))
    with(DF, 
        text(5825, 1.76, labels='Turn-off',
            pos=4, cex=0.8*text.cex, col='black'))
    with(DF, 
        text(5825, 1.64, labels='~ 2 Gyr',
            pos=4, cex=0.8*text.cex, col='black'))
    with(DF, 
        text(5670, 1.79, labels='Sub-giant',
            pos=4, cex=0.8*text.cex, col='black'))
    with(DF, 
        text(5670, 1.67, labels='~ 0.9 Gyr',
            pos=4, cex=0.8*text.cex, col='black'))
    
    segments(5779, 0.845, 5775, 0.935, lwd=1.2, col=track.col)
    with(DF, 
        text(5780, 0.88, labels='4.57 Gyr',
            pos=1, cex=0.8*text.cex, col='black'))
    par(family=font)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    axis(2, pretty(ylim, n=3), tick=T, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    axis(1, pretty(xlim, n=3), tick=T, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    #box(lwd=par()$lwd)
    
    if (make.y) {
        mtext(expression(L/L), 2, 2.1, outer=F, las=0, cex=text.cex)
        par(xpd=NA)
        points(5891, 1.542, pch=20, cex=0.75/3, lwd=1)
        points(5891, 1.542, pch=1,  cex=0.75,   lwd=1)
        par(xpd=F)
    }
    if (make.x) mtext(expression(T["eff"]/K), 1, 2.1, 
        outer=F, cex=text.cex)
    
}

make_plots(plot_HR, 'SG_US_step_sun', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, 
    wide=F, tall=F, 
    paper_pdf_height=4.17309*1.3,
    cex.paper=0.95,
    use.cairo=T, font='Palatino Linotype') 


