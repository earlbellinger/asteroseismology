#### Plot [Fe/H] vs Teff for the LEGACY+KAGES stars 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

require(magicaxis)
source(file.path('..', 'scripts', 'utils.R'))
library(RColorBrewer)

get_star <- function(filename) {
    KIC <- as.numeric(strsplit(basename(filename), '_perturb.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    if (!file.exists(file.path('learn-final', 'covs-SG_US_step', 'final',
        paste0(KIC, '.dat')))) return(NULL)
    print(filename)
    perturb <- read.table(filename, header=1)
    with(perturb,
            data.frame(KIC=KIC,
                       Teff=median(Teff), e_Teff=mad(Teff),
                       FeH=median(Fe.H),  e_FeH=mad(Fe.H)))
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(directory=file.path(file.path('perturb', 'feh')))


plot_unc <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(matrix(c(2,0,1,3), 2, 2, byrow=T), c(3,1), c(1,3))
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    oldpar <- par()
    par(mar=c(oldpar$mar[1], oldpar$mar[2], 0, 0)) 
    
    xlim <- c(-1, 0.5)
    ylim <- c(4900, 6800)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    h1 <- hist(DF$FeH, breaks=seq(xlim[1], xlim[2], 0.1), plot=F)
    h2 <- hist(DF$Teff, breaks=seq(ylim[1], ylim[2], 100), plot=F)
    top <- max(h1$counts, h2$counts)
    #k <- kde2d(df$x, df$y, n=25)
    
    abline(h=5772.8, lty=2, lwd=1.5)
    abline(v=0, lty=2, lwd=1.5)
    points(0, 5772.8, pch=21, cex=1.5, lwd=1.5, col='white', bg='white')
    
    #image(k, col=r) #plot the image
    points(DF$FeH, DF$Teff, 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white') 
    
    ## solar symbol
    points(0, 5772.8, pch=20, cex=0.75, lwd=1.5)
    points(0, 5772.8, pch=1,  cex=1.5,  lwd=1.5)
    
    x.mid <- 0.66 * xlim[2]
    y.mid <- 0.085 * (ylim[2] - ylim[1]) + ylim[1]
    
    segments(x.mid+median(DF$e_FeH), y.mid, 
             x.mid-median(DF$e_FeH), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(DF$e_Teff), 
             x.mid, y.mid-median(DF$e_Teff), 
        lwd=1, col='black')
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    axis(2, pretty(ylim, n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    axis(1, pretty(xlim, n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    mtext(expression(T['eff']/K), 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('[Fe/H]',               1, 2,   outer=F,        cex=text.cex)
    
    par(mar=c(0, oldpar$mar[2], 0.5, 0))
    barplot(h1$counts, axes=F, 
        ylim=c(0, top), space=0, 
        col=adjustcolor(orange, alpha.f=0.75),
        xaxs='i', yaxs='i')
    
    par(mar=c(oldpar$mar[1], 0, 0, 0.5))
    barplot(h2$counts, axes=F, 
        xlim=c(0, top), space=0, 
        col=adjustcolor(orange, alpha.f=0.75), 
        horiz=T,
        xaxs='i', yaxs='i')
    
}

make_plots(plot_unc, 'teffehist', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.4,
    cex.paper=0.95) 




plot_unc <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(matrix(c(2,0,1,3), 2, 2, byrow=T), c(4,1), c(1,4))
    
    par(mar=mar+c(0.3, -0.2, .1, .3), lwd=1.5)
    oldpar <- par()
    par(mar=c(oldpar$mar[1], oldpar$mar[2], 0, 0)) 
    
    ylim <- c(-1.1, 0.6)
    xlim <- rev(c(4800, 7000))
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    h1 <- hist(DF$Teff, breaks=seq(xlim[2], xlim[1], 200), plot=F)
    h2 <- hist(DF$FeH, breaks=seq(ylim[1], ylim[2], 0.1), plot=F)
    top <- max(h1$counts, h2$counts)
    
    
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
    
    rect(xlim[2], ylim[1], xlim[2]-100, ylim[2], col='white', border=NA)
    
    
    #rect(xlim[2], ylim[2], xlim[1], 0.1, 
    #    col=adjustcolor(orange, alpha.f=0.325), border=NA)
    #rect(xlim[2], 0.1, xlim[1], -0.1, 
    #    col=adjustcolor(red, alpha.f=0.2), border=NA)
    #rect(xlim[2], -0.1, xlim[1], -0.2, 
    #    col=adjustcolor(orange, alpha.f=0.2), border=NA)
    #rect(xlim[2], -0.2, xlim[1], ylim[1], 
    #    col=adjustcolor(blue, alpha.f=0.2), border=NA)
    #abline(h=-0.2, lty=2, lwd=1.5, col='darkgray')
    #abline(h=-0.1, lty=2, lwd=1.5, col='darkgray')
    #abline(h=0.1,  lty=2, lwd=1.5, col='darkgray')
    
    #abline(h=0.4437623, lwd=1.5, lty=3, col='darkgray')
    
    abline(v=5772.8, lty=2, lwd=1.5)
    abline(h=0,      lty=2, lwd=1.5)
    points(5772.8, 0, pch=21, cex=1.5, lwd=1.5, 
       #col='white', bg='white')
        col=cols[5], bg=cols[5])
    #"#bda186")
    #cols[5])
    
    # put on the points! 
    points(DF$Teff, DF$FeH, 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white') 
    
    ## solar symbol
    points(5772.8, 0, pch=20, cex=0.75, lwd=1.5)
    points(5772.8, 0, pch=1,  cex=1.5,  lwd=1.5)
    
    
    ## error bar in corner 
    x.mid <- 4975 #0.66 * xlim[1]
    y.mid <- 0.1 * (ylim[2] - ylim[1]) + ylim[1]
    
    segments(x.mid+median(DF$e_Teff), y.mid, 
             x.mid-median(DF$e_Teff), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(DF$e_FeH), 
             x.mid, y.mid-median(DF$e_FeH), 
        lwd=1, col='black')
    
    ## axes 
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    axis(2, pretty(ylim, n=3), tick=F, 
        cex.axis=1.2*text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    axis(1, pretty(xlim, n=3), tick=F, 
        cex.axis=1.2*text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    mtext('[Fe/H]',               2, 2.5, outer=F, las=0, cex=text.cex)
    mtext(expression(T['eff']/K), 1, 2,   outer=F,        cex=text.cex)
    
    par(family="Helvetica", xpd=NA)
    text(6500, ylim[2], labels='F', pos=1, cex=text.cex)
    text(5600, ylim[2], labels='G', pos=1, cex=text.cex)
    text(5000, ylim[2], labels='K', pos=1, cex=text.cex)
    par(family=font, xpd=F)
    
    #abline(h=h2$breaks, lty=2, lwd=1.5, col='gray')
    
    ## bar plots 
    par(mar=c(0, oldpar$mar[2], 0.1, 0))
    barplot(rev(h1$counts), axes=F, 
        ylim=c(0, top), space=0, 
        col=adjustcolor(blue, alpha.f=0.2),
        xaxs='i', yaxs='i')
    
    par(mar=c(oldpar$mar[1], 0, 0, 0))
    barplot(h2$counts, axes=F, 
        xlim=c(0, top), space=0, 
        col=adjustcolor(blue, alpha.f=0.2), 
        horiz=T, 
        ylim=c(0, length(h2$counts)), #+ 
            #if (h2$breaks[length(h2$breaks)] != xlim[2]) 0.5 else 0),
        xaxs='i', yaxs='i')
    
}
make_plots(plot_unc, 'teffehist', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    #paper_pdf_height=4.17309*1.4,
    paper_pdf_width=4.17309*1.385,
    paper_pdf_height=4.17309*1.2,#31,
    cex.paper=0.75)#0.95) 

