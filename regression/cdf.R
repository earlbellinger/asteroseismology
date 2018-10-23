#### Posterior CDFs for Kepler targets 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

require(magicaxis)
source(file.path('..', 'scripts', 'utils.R'))
library(RColorBrewer)

col.pal <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")

#DF <- read.table(file.path('feh', 'learn-feh', 'tables-SG_US_step', 'feh.dat'),
#    header=1, fill=T)
DF <- read.table(file.path('final', 'tables-SG_US_step', 'final.dat'),
    header=1, fill=T)
DF <- DF[-which(DF$Name == '5774694'),] # clip fake Sun 
#DF <- DF[(-nrow(DF)):(-nrow(DF)+5),] # clip Sun etc. 
DF <- DF[grepl('^\\d+$', DF$Name),] # clip Sun etc. 

#DF <- DF[,-which(colnames(DF) == "L" | colnames(DF) == "e_L")]
#col.pal <- adjustcolor(colorRampPalette(col.pal)((ncol(DF)-1)/2))
col.pal <- adjustcolor(colorRampPalette(col.pal)(ncol(DF)/2))

legend.names <- c("mass", "initial helium", "initial Z", 
    expression(alpha["MLT"]), expression(alpha["ov"]), expression(alpha["us"]),
    "diffusion factor", "age", "core hydrogen", "surface gravity", 
    "mean density", "luminosity", "radius", "surface helium")

group <- names(DF[,-1])[-grep('^e_', names(DF[,-1]))] 
unc <- DF[paste0('e_', group)] / DF[group] * 100 
sorted <- order(sapply(unc, median)) 
unc <- unc[,sorted] 
print(sapply(unc, median)) 
print(sapply(unc, fivenum)) 
legend.names <- legend.names[sorted]
legend.pcts <- as.expression(paste0('(', signif(sapply(unc, median), 2), '%)'))
legend.tot <- sapply(1:length(legend.names), function(ii) 
    bquote(.(legend.names[[ii]])~.(legend.pcts[[ii]])))

plot_cdfs <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(family=font, mar=mar+c(0.5, 0, 0, 9), lwd=1.5) 
    xlim <- c(0.05, 200) 
    #xticks <- c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)
    #xticks <- c(0.1, 0.5, 2, 10, 50) 
    xticks <- c(0.1, 1, 10, 100) 
    
    ylim <- c(0, 1)
    oh.one <- pretty(ylim)
    y.nums <- round(seq(0, nrow(unc), length.out=length(oh.one)))
    
    plot(NA, axes=F, log='x',
         xaxs='i', yaxs='i', 
         xlab="", ylab="",
         xlim=xlim, ylim=ylim)
    
    par(xpd=T)
    
    #segments(xlim[1], oh.one, xlim[2], oh.one, 
    #    lty=3, lwd=1.66, col='darkgray')
    #segments(xticks, 0, xticks, 1, 
    #    lty=3, lwd=1.66, col='darkgray')
    #grid(0, NULL, lty=2, col="darkgray", lwd=1.66) 
    #abline(v=xticks, lty=2, col = "darkgray", lwd=1.66) 
    
    for (xtick in 10**pretty(log10(xticks), n=10))
        for (ytick in pretty(ylim, n=20)) 
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, #lwd=0, 
                bg='darkgray', col='darkgray')
    
    for (ytick in oh.one[-1]) 
        for (xtick in 10**pretty(log10(c(0.05, 150)), n=30)[-1]) 
            #10**seq(log10(xlim[1]), log10(xlim[2]), length.out=14))
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, #lwd=0, 
                bg='darkgray', col='darkgray')
    
    #for (ii in 1:length(group)) lines(ecdf(unc[,ii]), col=1, cex=0.5)
    
    for (ii in 1:length(group)) {
        res <- environment(ecdf(unc[,ii]))
        lines(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            lwd=3, col=1)
        lines(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            lwd=1.5, col=col.pal[ii])
        points(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            pch=20, cex=1.4, lwd=1, col=1)
        points(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            pch=20, cex=1, lwd=1, col=col.pal[ii])
        #lines(ecdf(unc[,ii]), col=1,           cex=0.8, col.01line=NULL)
        #lines(ecdf(unc[,ii]), col=col.pal[ii], cex=0.5, col.01line=NULL)
    }
    
    #rect(0.0001, -0.1, 0.06, 0.1, col='white', border='white', density=1)
    segments(0.0001, 0, 0.06, 0, col='white', lwd=3)
    segments(0.0001, -0.01, 1000, -0.01, col='white', lwd=3)
    
    par(family="Helvetica LT Std Light")
    legend("topright", legend=as.expression(legend.tot), #text.font="Helvetica", 
           inset=c(-0.43, -0.08),#c(-0.53, 0),
           bty='n', pch=21, lty=F, 
           col=1, 
           pt.bg=col.pal,
           x.intersp=0.8,
           cex=0.85*text.cex)
    par(xpd=F, family=font)
    
    #magaxis(1:2, labels=F, tcl=-0.25, las=1, lwd.ticks=1.66,
    #        mgp=mgp+c(0, 0.4, 0), family=font, cex.axis=text.cex,
    #        use.par=T)
    magaxis(1, las=1, lwd.ticks=1.66, labels=F, tcl=tcl, 
            mgp=mgp+c(0, 0.4, 0), family=font, cex.axis=text.cex)
    #magaxis(2, las=1, lwd.ticks=1.66, labels=F, tcl=tcl, minorn=0, 
    #        mgp=mgp+c(0, 0.4, 0), family=font, cex.axis=text.cex)
    
    #minor.x <- signif(10**pretty(log10(xlim), n=length(xticks)*5), 2)
    #minor.x <- signif(10**pretty(log10(c(0.01, xticks, 1000)), 
    #    n=length(xticks)*7), 1)
    #minor.x <- minor.x[!minor.x %in% xticks]
    
    #axis(1, at=minor.x, 
    #    tick=T, labels=F, lwd=1.66, tcl=tcl/2)
    
    axis(1, at=xlim, tick=T, labels=F, tcl=0, lwd=1.66)
    
    axis(1, xticks, paste0(xticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, #tcl, 
        mgp=mgp+c(0, 0.3, 0), lwd=1.66)
    
    axis(2, oh.one, y.nums, 1,
         cex.axis=text.cex, tcl=tcl, 
         las=1, mgp=mgp+c(0, 0.4, 0), lwd=1.66)
    
    par(mgp=mgp+c(0.8, 0, 0))
    title(xlab=expression("Relative Uncertainty"))
    par(mgp=mgp+c(1, 0, 0))
    title(ylab=expression("Star"))
}

make_plots(plot_cdfs, paste0('cdf'), filepath='plots/feh', 
    mar=utils.mar + c(0, -1, 0, 2), 
    paper_pdf_height=4.17309*1.45,
    #paper_pdf_width=6.97522*1.8, 
    cex.paper = 0.95, make_png=F, 
    slides=F, wide=T, tall=F, thin=F, 
    use.cairo=T, font="Palatino Linotype")





#col.pal <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")
#col.pal <- c('#34738f', '#122f3d', '#be3e2b', 
#    '#ed8a45', '#eec70e')

get_star <- function(filename) {
    KIC <- as.numeric(strsplit(basename(filename), '_perturb.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    if (!file.exists(file.path('learn-final', 'covs-SG_US_step', 'final',
        paste0(KIC, '.dat')))) return(NULL)
    print(filename)
    perturb <- read.table(filename, header=1)
    perturb.DF <- with(perturb,
            data.frame(KIC=KIC,
                       Teff=median(Teff), e_Teff=mad(Teff),
                       FeH=median(exp(Fe.H)),  
                       e_FeH=median(exp(Fe.H) * mad(Fe.H)),
                       nu_max=median(nu_max), e_nu_max=mad(nu_max)))
    if ('L' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(L=median(L), e_L=mad(L))))
    if ('Dnu0' %in% names(perturb)) 
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(Dnu=median(Dnu0),  e_Dnu=mad(Dnu0))))
    if ('dnu02' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(dnu02=median(dnu02), e_dnu02=mad(dnu02))))
    if ('dnu13' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(dnu13=median(dnu13), e_dnu13=mad(dnu13))))
    #if ('r02' %in% names(perturb))
    #    perturb.DF <- cbind(perturb.DF, 
    #        with(perturb, data.frame(r02=median(r02), e_r02=mad(r02))))
    #if ('r13' %in% names(perturb))
    #    perturb.DF <- cbind(perturb.DF, 
    #        with(perturb, data.frame(r13=median(r13), e_r13=mad(r13))))
    if ('r01' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(r01=median(r01), e_r01=mad(r01))))
    #if ('r10' %in% names(perturb))
    #    perturb.DF <- cbind(perturb.DF, 
    #        with(perturb, data.frame(r10=median(r10), e_r10=mad(r10))))
    perturb.DF
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(directory=file.path(file.path('perturb', 'feh')))

#col.pal <- c('#1C3144', '#D00000', '#FFBA08', '#A2AEBB', '#3F88C5')
#col.pal <- c('#335C67', '#FFF3B0', '#E09F3E', '#9E2A2B', '#540B0E')
#col.pal <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")
col.pal <- rev(c("#F29559", "#DB4D48", "#323031", blue))
col.pal <- c(brewer.pal(11, 'Spectral')[-c(6,7)], "#512854")
#col.pal <- adjustcolor(colorRampPalette(col.pal)((ncol(DF)-1)/2))
col.pal <- adjustcolor(colorRampPalette(col.pal)((ncol(DF))/2))
#col.pal <- adjustcolor(colorRampPalette(col.pal)((ncol(DF)-1)/2))
##col.pal[length(col.pal)] <- adjustcolor('darkgreen', alpha.f=0.75)#'#358F56'
#col.pal[length(col.pal)] <- 'brown'#'purple'#'#156863'
#col.pal <- adjustcolor(colorRampPalette(col.pal)(ncol(DF)/2))

legend.names <- c(expression(T['eff']), 
    '[Fe/H]', 
    expression(nu['max']),
    #expression(L)
    expression(Delta*nu), 
    expression(delta*nu["0,2"]), 
    #expression(r["0,2"]),
    #expression(r["0,1"]), 
    expression(r["1,0"]),
    expression(delta*nu["1,3"]))#, 
    #expression(r["1,3"]))

group <- names(DF[,-1])[-grep('^e_', names(DF[,-1]))] 
unc <- abs(DF[paste0('e_', group)] / DF[group]) * 100 
sorted <- order(sapply(unc, function (x) median(x, na.rm=T))) 
unc <- unc[,sorted] 
print(sapply(unc, median)) 
print(sapply(unc, fivenum)) 
legend.names <- legend.names[sorted]
legend.pcts <- as.expression(paste0('(', 
    signif(sapply(unc, function (x) median(x, na.rm=T)), 2), '%)'))
legend.tot <- sapply(1:length(legend.names), function(ii) 
    bquote(.(legend.names[[ii]])~.(legend.pcts[[ii]])))

plot_cdfs <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(family=font, mar=mar+c(0.5, 0, 0, 9), lwd=1.5) 
    xlim <- c(0.01, 100) 
    xticks <- c(0.01, 0.1, 1, 10, 100) 
    
    ylim <- c(0, 1)
    oh.one <- pretty(ylim)
    y.nums <- round(seq(0, nrow(unc), length.out=length(oh.one)))
    
    plot(NA, axes=F, log='x',
         xaxs='i', yaxs='i', 
         xlab="", ylab="",
         xlim=xlim, ylim=ylim)
    
    par(xpd=T)
    
    for (xtick in 10**pretty(log10(xticks), n=10))
        for (ytick in pretty(ylim, n=20)) 
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, 
                bg='darkgray', col='darkgray')
    
    for (ytick in oh.one[-1]) 
        for (xtick in 10**pretty(log10(c(0.01, 100)), n=30)[-1]) 
            points(xtick, ytick, pch=19, cex=0.1, lwd=0.5, 
                bg='darkgray', col='darkgray')
    
    groups <- c(1:length(group))[order(apply(apply(unc, 2, is.na), 2, sum))]
    for (ii in groups) {
        res <- environment(ecdf(unc[,ii]))
        lines(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            lwd=3, col=1)
        lines(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            lwd=1.5, col=col.pal[ii])
        points(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            pch=20, cex=1.4, lwd=1, col=1)
        points(res$x, res$y*(sum(!is.na(unc[,ii])) / nrow(unc)), 
            pch=20, cex=1, lwd=1, col=col.pal[ii])
        #lines(ecdf(unc[,ii])(seq(0, 1, 0.01)), 
        #      seq(0, 1, 0.01)*(sum(!is.na(unc[,ii])) / nrow(unc)),
        #      col=1,           cex=0.8, col.01line=NULL)
        #lines(ecdf(unc[,ii])(seq(0, 1, 0.01)), 
        #      seq(0, 1, 0.01)*(sum(!is.na(unc[,ii])) / nrow(unc)),
        #      col=col.pal[ii], cex=0.5, col.01line=NULL)
    }
    
    segments(xlim[1]/10, 0, xlim[1], 0, col='white', lwd=3)
    segments(xlim[1]/10, -0.01, xlim[2], -0.01, col='white', lwd=3)
    
    par(family="Helvetica LT Std Light")
    legend("topright", legend=as.expression(legend.tot), 
           inset=c(-0.285, 0.15), 
           bty='n', pch=21, lty=F, 
           col=1, 
           pt.bg=col.pal,
           x.intersp=0.8,
           cex=0.85*text.cex)
    par(xpd=F, family=font)
    
    magaxis(1, las=1, lwd.ticks=1.66, labels=F, tcl=tcl, 
            mgp=mgp+c(0, 0.4, 0), family=font, cex.axis=text.cex)
    
    axis(1, at=xlim, tick=T, labels=F, tcl=0, lwd=1.66)
    
    axis(1, xticks, paste0(xticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, #tcl, 
        mgp=mgp+c(0, 0.3, 0), lwd=1.66)
    
    axis(2, oh.one, y.nums, 1,
         cex.axis=text.cex, tcl=tcl, 
         las=1, mgp=mgp+c(0, 0.4, 0), lwd=1.66)
    
    par(mgp=mgp+c(0.8, 0, 0))
    title(xlab=expression("Relative Uncertainty"))
    par(mgp=mgp+c(1, 0, 0))
    title(ylab=expression("Star"))
}

make_plots(plot_cdfs, paste0('cdf-inputs'), filepath='plots/feh', 
    mar=utils.mar + c(0, -1, 0, 2), 
    paper_pdf_height=4.17309*1.45,
    #paper_pdf_width=6.97522*1.8, 
    cex.paper = 0.95, make_png=F, 
    slides=F, wide=T, tall=F, thin=F, 
    use.cairo=T, font="Palatino Linotype")

