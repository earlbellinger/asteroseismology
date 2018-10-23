#### Assessing the impact of systematic errors in [Fe/H] and Teff together 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 
library(parallelMap)
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

set.seed(0)

dirs <- list.files('both', full.names=T)
baseline <- file.path('feh', 'learn-feh')
bias_dirs <- dirs[grep('bias', dirs)]
imp_dirs <- dirs[grep('imp', dirs)]

feh.biases <- sapply(strsplit(basename(bias_dirs), '_'), function(x) {
    pos <- grepl('feh', x)
    if (!any(pos)) return(0)
    as.numeric(strsplit(x[pos], 'fehbias')[[1]][2])
})
teff.biases <- sapply(strsplit(basename(bias_dirs), '_'), function(x) {
    pos <- grepl('teff', x)
    if (!any(pos)) return(0)
    as.numeric(strsplit(x[pos], 'teffbias')[[1]][2])
})
feh.imps <- sapply(strsplit(basename(imp_dirs), '_'), function(x) {
    pos <- grepl('feh', x)
    if (!any(pos)) return(0)
    as.numeric(strsplit(x[pos], 'fehimp')[[1]][2])
})
teff.imps <- sapply(strsplit(basename(imp_dirs), '_'), function(x) {
    pos <- grepl('teff', x)
    if (!any(pos)) return(0)
    as.numeric(strsplit(x[pos], 'teffimp')[[1]][2])
})

#col.pal <- c("#DB4D48", "#F29559", blue, "#323031")
col.pal <- c("#DB4D48", orange, blue, "#323031")
params <- c('Radius', 'Mass', 'Density', 'Age')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return (NULL)
    
    #perturb <- read.table(file.path('perturb', 'feh', 
    #    sub('\\.', '_perturb.', basename(filename))), header=1)
    #if (any(perturb$Fe) > 0.44) return(NULL)
    
    covs <- read.table(filename, header=T)
    n <- nrow(covs)
    covs <- covs[complete.cases(covs),]
    #print(nrow(covs))
    if (nrow(covs) < n/2) return(NULL) 
    with(covs, 
        data.frame(KIC=KIC,
                   R=median(radius),    e_R=mad(radius),
                   M=median(M),         e_M=mad(M),
                   rho=median(density), e_rho=mad(density),
                   age=median(age),     e_age=mad(age)
        )
    )
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    parallelMap(get_star, list.files(directory, recursive=T, full.names=T)))

baseline.DF <- get_DF(file.path(baseline))

biases.list <- parallelMap(get_DF, bias_dirs)
imps.list <- parallelMap(get_DF, imp_dirs)

unc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    if (ii == 0) {
        teff.imp <- 0
        feh.imp <- 0
        imp.DF <- baseline.DF
    } else {
        teff.imp <- teff.imps[ii]
        feh.imp <- feh.imps[ii]
        imp.DF <- imps.list[[ii]]
    }
    print(nrow(imp.DF))
    #if (nrow(imp.DF) < 50) return(NULL)
    with(imp.DF, data.frame(
        teff.imp=teff.imp,
        feh.imp=feh.imp,
        R=mean(e_R/R),             e_R=mad(e_R/R),
        M=mean(e_M/M),             e_M=mad(e_M/M),
        rho=mean(e_rho/rho),     e_rho=mad(e_rho/rho),
        age=mean(e_age/age),       e_age=mad(e_age/age),
        N=nrow(imp.DF)))
}, ii=0:length(teff.imps)))
unc <- unc[with(unc, order(teff.imp, feh.imp)),]
#unc <- unc[-which(unc$feh.imp == 0.8),]
unc[unc$N < 97,]

acc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    if (ii == 0) {
        teff.bias <- 0
        feh.bias <- 0
        bias.DF <- baseline.DF
    } else {
        teff.bias <- teff.biases[ii]
        feh.bias <- feh.biases[ii]
        bias.DF <- biases.list[[ii]]
    }
    DF <- merge(baseline.DF, bias.DF, by='KIC')
    print(nrow(bias.DF))
    with(DF, data.frame(
        teff.bias=teff.bias,
        feh.bias=feh.bias,
        R=mean(abs(R.x - R.y)/R.x),
        d_R=mad(abs(R.x - R.y)/R.x),
        M=mean(abs(M.x - M.y)/M.x),
        d_M=mad(abs(M.x - M.y)/M.x),
        rho=mean(abs(rho.x - rho.y)/rho.y),
        d_rho=mad(abs(rho.x - rho.y)/rho.y),
        age=mean(abs(age.x - age.y)/age.y),
        d_age=mad(abs(age.x - age.y)/age.y),
        N=nrow(bias.DF)))
}, ii=0:length(teff.biases)))
acc <- acc[with(acc, order(teff.bias, feh.bias)),]
acc[acc$N < 97,]


plot_unc <- function(DF, ii, lab, xlab, ylab, lvls, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- range(DF[,1])
    ylim <- range(DF[,2])
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    param <- DF[,2*ii+1] * 100
    
    contour(unique(DF[,1]), 
            unique(DF[,2]), 
            matrix(param, 
                nrow=length(unique(DF[,1])), 
                ncol=length(unique(DF[,2]))), 
        lwd=1.5, add=T, col=col.pal[ii], labcex=0.9*text.cex,
        levels=lvls)
    
    par(xpd=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) axis(1, pretty(DF[,1], n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    if (make.y) axis(2, pretty(DF[,2], n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    shadowtext(xlim[1], if (xlim[1]==-100) 0.075 else 0.0875, 
        labels=lab, 
        col=1, bg='white', 
        cex=0.0*text.cex, pos=4)
    par(family=font, xpd=F)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=F, las=0, cex=text.cex)
    if (make.x) mtext(xlab, 1, 2, outer=F, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('both_imp', ii), 
    DF=unc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Added Uncertainty [dex]", 
    xlab="Added Uncertainty [K]",
    lab=params[ii],
    lvls=list(
        seq(1, 2, 0.05),
        seq(3, 6, 0.2),
        seq(1, 3, 0.05),
        seq(10, 18, 0.5)
         )[[ii]],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}


for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('both_bias', ii), 
    DF=acc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Systematic Error [dex]", 
    xlab="Systematic Error [K]",
    lab=params[ii],
    lvls=list(
        seq(0, 2, 0.1),
        seq(0, 3, 0.5),
        seq(0, 1, 0.1),
        seq(0, 10, 1)
         )[[ii]],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}

if(F) {
plot_unc <- function(DF, ii, lab, xlab, ylab, lvls, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- range(DF[,1])
    ylim <- range(DF[,2])
    
    param <- DF[,2*ii+1] * 100
    
    filled.contour(unique(DF[,1]), 
            unique(DF[,2]), 
            matrix(param, 
                nrow=length(unique(DF[,1])), 
                ncol=length(unique(DF[,2]))), 
        levels=lvls,
        xaxs='i', yaxs='i',
        key.axes={
            axis(4, cex.axis=text.cex, tcl=0, line=0, mgp=mgp)
            mtext('Explained Variance', side=4, las=3, line=2, cex=text.cex)
        },
        plot.axes={
            contour(
                unique(DF[,1]), 
                unique(DF[,2]), 
                matrix(param, 
                    nrow=length(unique(DF[,1])), 
                    ncol=length(unique(DF[,2]))), 
                add=TRUE, 
                levels=lvls,
                labcex=0.8*text.cex)
            
            magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
                family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
            magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
                family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
            
            if (make.x) axis(1, pretty(DF[,1], n=3), tick=F, 
                cex.axis=text.cex, tcl=0, las=1, 
                mgp=mgp+c(0, 0.25, 0), lwd=1.5)
            if (make.y) axis(2, pretty(DF[,2], n=3), tick=F, 
                cex.axis=text.cex, tcl=0, las=1, 
                mgp=mgp+c(0, 0.35, 0), lwd=1.5)
            
            
            par(family="Helvetica", xpd=NA)
            shadowtext(xlim[1], if (xlim[1]==-100) 0.075 else 0.0875, 
                labels=lab, 
                col=1, bg='white', 
                cex=0.0*text.cex, pos=4)
            par(family=font, xpd=F)
            
            par(xpd=NA)
            rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
            par(xpd=F)
        },
        plot.title={
            #title(xlab='Minimum CZ Depth', cex.lab=text.cex, line=2)
            #title(ylab='Maximum CZ Depth', cex.lab=text.cex, line=2)
            if (make.y) mtext(ylab, 2, 2.5, outer=F, las=0, cex=text.cex)
            if (make.x) mtext(xlab, 1, 2, outer=F, cex=text.cex)
        })
    
    #box(lwd=par()$lwd)
    
}
}




plot_unc <- function(DF, ii, lab, xlab, ylab, lvls, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- range(DF[,1])
    ylim <- range(DF[,2])
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    param <- DF[,2*ii+1] * 100
    #unc <- DF[,2*ii+2] 
    
    #param <- signif(param, nchar(sapply(signif(unc, 2), toString))-3) * 100
    
    .filled.contour(
            unique(DF[,1]), 
            unique(DF[,2]), 
            matrix(param, 
                nrow=length(unique(DF[,1])), 
                ncol=length(unique(DF[,2]))), 
        #lwd=1.5, #add=T, #col=col.pal[ii], 
        #labcex=0.9*text.cex,
        levels=lvls,
        col=colorRampPalette(c('white', 'gray60'))(length(lvls)+1))
    
    if (xlim[1]==-100 && F) {
        segments(pretty(xlim, n=10)-2,
                 0,
                 pretty(xlim, n=10)+2,
                 0,
             col=adjustcolor(1, alpha.f=0.1), lwd=1.5, lty=1)
        segments(0,
                 pretty(ylim, n=10)-0.005,
                 0,
                 pretty(ylim, n=10)+0.005,
             col=adjustcolor(1, alpha.f=0.1), lwd=1.5, lty=1)
    }
    
    contour(unique(DF[,1]), 
            unique(DF[,2]), 
            matrix(param, 
                nrow=length(unique(DF[,1])), 
                ncol=length(unique(DF[,2]))), 
        lwd=1.5, add=T, col=col.pal[ii], labcex=0.9*text.cex,
        levels=lvls)
    
    par(xpd=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) axis(1, pretty(DF[,1], n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    if (make.y) axis(2, pretty(DF[,2], n=3), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    shadowtext(xlim[1]- if (xlim[1]==-100) 2 else 1, 
            if (xlim[1]==-100) 0.08 else 0.09, 
        labels=lab, 
        col=1, bg='white', 
        r=0.2,
        cex=0.0*text.cex, pos=4)
    par(family=font, xpd=F)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=F, las=0, cex=text.cex)
    if (make.x) mtext(xlab, 1, 2, outer=F, cex=text.cex)
}


for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('both_imp', ii), 
    DF=unc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Added Uncertainty [dex]", 
    xlab="Added Uncertainty [K]",
    lab=params[ii],
    lvls=pretty(unc[,2*ii+1] * 100, n=if (ii == 3) 10 else 10),
    #lvls=list(
    #    pretty(unc[,2*ii+1] * 100, n=8), #seq(1, 4, 0.05),
    #    #seq(3, 10, 0.2),
    #    #seq(1, 4, 0.05),
    #    #seq(10, 30, 0.5)
    #     )[[ii]],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}


for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('both_bias', ii), 
    DF=acc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Systematic Error [dex]", 
    xlab="Systematic Error [K]",
    lab=params[ii],
    lvls=pretty(acc[,2*ii+1] * 100, n=8),
    #lvls=list(
    #    seq(0, 2, 0.1),
    #    seq(0, 6, 0.5),
    #    seq(0, 1.2, 0.1),
    #    seq(0, 20, 1)
    #     )[[ii]],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}

