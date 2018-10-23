#### Assessing the impact of systematic errors in [Fe/H]
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 
library(parallelMap)
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

use.spline <- F

set.seed(0)

dirs <- list.files('feh', full.names=T)
baseline <- file.path('feh', 'learn-feh')
bias_dirs <- dirs[grep('bias', dirs)]
imp_dirs <- dirs[grep('imp', dirs)]

biases <- as.numeric(sub('learn-feh_bias', '', basename(bias_dirs)))
imps <- as.numeric(sub('learn-feh_imp', '', basename(imp_dirs)))

col.pal <- c("#DB4D48", "#F29559", blue, "#323031")
params <- c('Radius', 'Mass', 'Density', 'Age')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    
    #perturb <- read.table(file.path('perturb', 'feh', 
    #    sub('\\.', '_perturb.', basename(filename))), header=1)
    #if (any(perturb$Fe.H > 0.44)) return(NULL)
    
    covs <- read.table(filename, header=T)
    n <- nrow(covs)
    covs <- covs[complete.cases(covs),]
    #print(nrow(covs))
    if (nrow(covs) < n/2) return(NULL) 
    #covs <- covs[1:(min(nrow(covs), n)),]
    with(covs, 
        data.frame(KIC=KIC,
                   R=median(radius),    e_R=mad(radius),
                   M=median(M),         e_M=mad(M),
                   #alpha=median(alpha), e_alpha=mad(alpha),
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

#biases.list <- list(biases.list, baseline.DF)
#biases <- c(biases, 0)

unc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    print(ii)
    if (ii == 0) {
        imp <- 0
        imp.DF <- baseline.DF
    } else {
        imp <- imps[ii]
        imp.DF <- imps.list[[ii]]
    }
    #if (nrow(imp.DF) < 50) return(NULL)
    with(imp.DF, data.frame(imp=imp,
        R=mean(e_R),             e_R=sd(e_R),
        M=mean(e_M),             e_M=sd(e_M),
        #alpha=mean(e_alpha), e_alpha=sd(e_alpha),
        rho=mean(e_rho),       e_rho=sd(e_rho),
        age=mean(e_age),       e_age=sd(e_age),
        N=nrow(imp.DF)))
}, ii=0:length(imps)))
unc <- unc[order(unc$imp),]
baseline.uncR <- signif(unc[1,]$R, 1)
baseline.uncM <- signif(unc[1,]$M, 1)
baseline.uncage <- signif(unc[1,]$age, 1)

diff.1 <- unc[which(unc$imp == 0.1),]-unc[1,]
diff.R <- signif(diff.1$R, 1)
diff.M <- signif(diff.1$M, 1)
diff.age <- signif(diff.1$age, 1)
#unc[7,]-unc[1,]

unc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    if (ii == 0) {
        imp <- 0
        imp.DF <- baseline.DF
    } else {
        imp <- imps[ii]
        imp.DF <- imps.list[[ii]]
    }
    #if (nrow(imp.DF) < 50) return(NULL)
    with(imp.DF, data.frame(imp=imp,
        R=mean(e_R/R),             e_R=sd(e_R/R),
        #R_u=quantile(e_R/R, 0.75),   R_l=quantile(e_R/R, 0.25),
        M=mean(e_M/M),             e_M=sd(e_M/M),
        #M_u=quantile(e_M/M, 0.75),   M_l=quantile(e_M/M, 0.25),
        rho=mean(e_rho/rho),     e_rho=sd(e_rho/rho),
        #rho_u=quantile(e_rho/rho, 0.75),   rho_l=quantile(e_rho/rho, 0.25),
        age=mean(e_age/age),       e_age=sd(e_age/age),
        #age_u=quantile(e_age/age, 0.75),   age_l=quantile(e_age/age, 0.25),
        N=nrow(imp.DF)))
}, ii=0:length(imps)))
unc <- unc[order(unc$imp),]
#abs(unc[7,]-unc[1,])/unc[1,]*100

diff.1 <- abs(unc[which(unc$imp == 0.1),]-unc[1,]) / unc[1,] * 100
rel.R <- signif(diff.1$R, 2)
rel.M <- signif(diff.1$M, 2)
rel.age <- signif(diff.1$age, 2)

acc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    if (ii == 0) {
        bias <- 0
        bias.DF <- baseline.DF
    } else {
        bias <- biases[ii]
        bias.DF <- biases.list[[ii]]
    }
    cat(paste(bias, nrow(bias.DF), '\n'))
    #if (nrow(bias.DF) < 50) return(NULL)
    DF <- merge(baseline.DF, bias.DF, by='KIC')
    with(DF, data.frame(bias=bias,
        R=mean(abs(R.x - R.y)/R.x),
        d_R=sd(abs(R.x - R.y)/R.x),
        #R_u=quantile(abs(R.x - R.y)/R.x, 0.75),
        #R_l=quantile(abs(R.x - R.y)/R.x, 0.25),
        M=mean(abs(M.x - M.y)/M.x),
        d_M=sd(abs(M.x - M.y)/M.x),
        #M_u=quantile(abs(M.x - M.y)/M.x, 0.75),
        #M_l=quantile(abs(M.x - M.y)/M.x, 0.25),
        rho=mean(abs(rho.x - rho.y)/rho.y),
        d_rho=sd(abs(rho.x - rho.y)/rho.y),
        #rho_u=quantile(abs(rho.x - rho.y)/rho.x, 0.75),
        #rho_l=quantile(abs(rho.x - rho.y)/rho.x, 0.25),
        age=mean(abs(age.x - age.y)/age.y),
        d_age=sd(abs(age.x - age.y)/age.y),
        #age_u=quantile(abs(age.x - age.y)/age.x, 0.75),
        #age_l=quantile(abs(age.x - age.y)/age.x, 0.25),
        N=nrow(bias.DF)))
}, ii=0:length(biases)))
acc <- acc[order(acc$bias),]
#acc[15,]*100
bias.1 <- acc[which(acc$bias == 0.1),]*100
bias.R <- signif(bias.1$R, 1)
bias.M <- signif(bias.1$M, 1)
bias.age <- signif(bias.1$age, 1)

acc <- acc[acc$N>40,]
signif(min(abs(with(acc[acc$bias>0,], splinefun(R, bias)(0.014))), 
           abs(with(acc[acc$bias<0,], splinefun(R, bias)(0.014)))), 2)
signif(min(abs(with(acc[acc$bias>0,], splinefun(M, bias)(0.035))),
           abs(with(acc[acc$bias<0,], splinefun(M, bias)(0.035)))), 2)
signif(min(abs(with(acc[acc$bias>0,], splinefun(age, bias)(0.11))),
           abs(with(acc[acc$bias<0,], splinefun(age, bias)(0.11)))), 2)

bias.R2 <- acc[order(acc$R),][min(which(acc[order(acc$R),]$R > baseline.uncR)),]$bias
bias.M2 <- acc[order(acc$M),][min(which(acc[order(acc$M),]$M > baseline.uncM)),]$bias
bias.age2 <- acc[order(acc$age),][min(which(acc[order(acc$age),]$age > baseline.uncage)),]$bias

cat(paste0("\n\nWe find that a systematic error of 0.1~dex in [Fe/H] measurements translates on average to differences of only ", bias.age, "\\%, ", bias.M, "\\%, and ", bias.R, "\\% in the resulting stellar ages, masses, and radii, respectively, which are smaller than the reported relative uncertainties for these quantities ($\\sim$~10\\%, 3\\%, 1\\%). \n We furthermore find that increasing the reported uncertainty of [Fe/H] measurements by 0.1~dex (100\\%) increases the uncertainties of stellar ages, masses, and radii on average by only ", diff.age, "~Gyr (", rel.age, "\\%)", ", ", diff.M, "~$\\text{M}_\\odot$ (", rel.M, "\\%), and ", diff.R, "~$\\text{R}_\\odot$ (", rel.R, "\\%)", ", respectively, which are well below the reported uncertainties for these estimates (", baseline.uncage, "~Gyr", ", ", baseline.uncM, "~$\\text{M}_\\odot$, ", baseline.uncR, "~$\\text{R}_\\odot$", ").\n\n"))









plot_unc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(0.0005, 0.1)#range(DF[,1])
    #ylim <- c(0, 20)
    ylim <- c(0, 21)
    #ylim <- c(-25, 25)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- c(0.001, 0.01, 0.1)
    yticks <- pretty(ylim)
    
    param <- DF[,2*ii]   * 100
    p.unc <- DF[,2*ii+1] * 100
    p.u <- param+p.unc
    p.l <- param-p.unc
    #param <- DF[,3*ii-1] * 100
    #p.u   <- DF[,3*ii]   * 100
    #p.l   <- DF[,3*ii+1] * 100
    
    if (use.spline) {
        
        interp <- 10**(seq(log10(min(DF[,1][DF[,1]>0])/10), 
                           log10(max(DF[,1]*10)), 
                        length.out=1000))
        param[param<ylim[1]] <- ylim[1]/10
            fit <- predict(smooth.spline(DF[,1], param, spar=0), interp)$y
        p.u.fit <- predict(smooth.spline(DF[,1], p.u, spar=0), interp)$y
        p.l.fit <- predict(smooth.spline(DF[,1], p.l, spar=0), interp)$y
        polygon(c(interp,      rev(interp)), 
                c(p.u.fit, rev(p.l.fit)), 
                col=adjustcolor(col.pal[ii], alpha.f=0.5), 
                border=NA)
        lines(interp, fit, lwd=3, col=col.pal[ii], lty=1)
        
    } else {
        
        param[param<ylim[1]] <- ylim[1]/10
        xs <- DF[,1]
        xs[xs<xlim[1]] <- xlim[1] / 10
        polygon(c(xs,      rev(xs)), 
                #c(param-p.unc, rev(param+p.unc)),
                c(p.l, rev(p.u)),
                col=adjustcolor(col.pal[ii], alpha.f=0.5), 
                border=NA)
        lines(xs, param, lwd=3, col=col.pal[ii], lty=1)
        
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
    par(xpd=F)
    
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.y) axis(2, yticks[c(T,F)], paste0(yticks[c(T,F)], '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    if (make.x) axis(1, c(0.001, 0.01, 0.1), c('0.001', '0.01', '0.1'), 
        tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    rect(10**(sum(log10(xlim))/2*.85), ylim[2]*0.72, 
         10**(sum(log10(xlim))/2*1.15), ylim[2]*0.9, 
        col='white', border=1, lwd=1.5)
    text(10**(sum(log10(xlim))/2), ylim[2]*0.9, 
        labels=lab, pos=1, cex=0.9*text.cex)
    par(family=font, xpd=F)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=F, las=0, cex=text.cex)
    if (make.x) mtext(xlab, 1, 2, outer=F, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('imp', ii), 
    DF=unc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Relative Uncertainty", 
    xlab="Added Uncertainty [dex]",
    lab=params[ii],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}









plot_acc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(matrix(1:2, ncol=2), width=1, heights=c(2,2), respect=FALSE)
    
    par(oma=mar+c(0.3, -0.3, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(0.1, 0.0001)
    ylim <- c(0, 4)
    
    xticks <- 10**pretty(log10(xlim))
    yticks <- pretty(ylim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim,#c(0.2, 0.0001), 
        ylim=ylim)
    
    param <- DF[,2*ii]   * 100
    p.unc <- DF[,2*ii+1] * 100
    p.u <- param+p.unc
    p.l <- param-p.unc
    
    #param <- DF[,3*ii-1] * 100
    #p.u   <- DF[,3*ii]   * 100
    #p.l   <- DF[,3*ii+1] * 100
    
    param[param<=ylim[1]] <- 0#ylim[1]#/10
    xs <- DF[,1]
    left  <- xs<=0
    right <- xs>=0
    
    xs[abs(xs) <= min(xlim)] <- min(xlim)#10
    
    if (use.spline) {
        
        interp <- 10**(seq(log10(min(DF[,1][DF[,1]>0])/10), 
                           log10(max(DF[,1])), 
                        length.out=1000))
        interp <- c(-rev(interp), 0, interp)
        
            fit <- predict(smooth.spline(xs, param, spar=0), interp)$y
        #unc.fit <- predict(smooth.spline(xs, p.unc, spar=1), interp)$y
        p.u.fit <- predict(smooth.spline(xs, p.u, spar=1), interp)$y
        p.l.fit <- predict(smooth.spline(xs, p.l, spar=1), interp)$y
        
        interp.left  <- abs(interp[interp < 0])
        fit.left     <- fit[interp < 0]
        p.u.fit.left <- p.u.fit[interp < 0]
        p.l.fit.left <- p.l.fit[interp < 0]
        polygon(c(interp.left,           rev(interp.left)), 
                #c(fit.left+unc.fit.left, rev(fit.left-unc.fit.left)), 
                c(p.u.fit.left, rev(p.l.fit.left)), 
                col=adjustcolor(col.pal[ii], alpha.f=0.5), 
                border=NA)
        lines(interp, fit, lwd=3, col=col.pal[ii], lty=1)
        
    } else {
        
        polygon(c(abs(xs[left]),        rev(abs(xs[left]))), 
                #c(c(param-p.unc)[left], rev(c(param+p.unc)[left])),
                c(p.l[left], rev(p.u[left])),
                col=adjustcolor(col.pal[ii], alpha.f=0.5), 
                border=NA)
        lines(abs(xs[left]), param[left], lwd=3, col=col.pal[ii], lty=1)
        
    }
    
    par(xpd=NA)
    rect(xlim[1], ylim[1], 
         xlim[2], ylim[1]-0.1, 
         col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, #minorn=10,
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, #minorn=5,
        family=font, las=1, majorn=5, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) axis(1, c(0.1, 0.001), c('-0.1', '-0.001'), 
        tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    if (make.y) axis(2, yticks, paste0(yticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=TRUE, las=0, cex=text.cex)
    
    ## RIGHT PANEL
    
    xlim <- rev(xlim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim,#c(0.2, 0.0001), 
        ylim=ylim)
    
    polygon(c(abs(xs[right]),      rev(abs(xs[right]))), 
            #c(c(param-p.unc)[right], rev(c(param+p.unc)[right])),
            c(p.l[right], rev(p.u[right])),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)#col.pal[ii], lwd=1.5)#NA)#col.pal[ii], lty=2, lwd=3)
    lines(xs[right], param[right], lwd=3, col=col.pal[ii], lty=1)
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]-0.1, 
         xlim[2]+1, ylim[2]+0.1, 
         col='white', border=NA)
    rect(xlim[1],  ylim[1], 
         xlim[2]*10,  ylim[1]-0.1, 
         col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, #minorn=10,
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) axis(1, c(0.001, 0.1), 
                c('0.001', '0.1'), 
            tick=F, cex.axis=text.cex, tcl=0, las=1, 
            mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    
    par(family="Helvetica", xpd=NA)
    rect(xlim[1]/6.94, ylim[2]*0.72, xlim[1]*6.94, ylim[2]*0.9, 
        col='white', border=1, lwd=1.5)
    text(xlim[1], ylim[2]*0.9, labels=lab, pos=1, cex=0.9*text.cex)
    par(family=font, xpd=F)
    if (make.x) mtext(xlab, 1, 2, outer=TRUE, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_acc, paste0('bias', ii), 
    DF=acc, ii=ii, 
    ylab="Relative Difference",
    xlab="Systematic Error [dex]",
    lab=params[ii],
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}










































plot_unc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(0.1, 500)
    ylim <- c(0, 36)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- c(0.001, 0.01, 0.1)#10**pretty(log10(xlim))
    yticks <- pretty(ylim)#10**pretty(log10(ylim)) #
    
    param <- DF[,2*ii]   * 100
    p.unc <- DF[,2*ii+1] * 100
    p.u <- param+p.unc
    p.l <- param-p.unc
    
    #param <- DF[,3*ii-1] * 100
    #p.u   <- DF[,3*ii]   * 100
    #p.l   <- DF[,3*ii+1] * 100
    
    param[param<ylim[1]] <- ylim[1]/10
    xs <- DF[,1]
    xs[xs<xlim[1]] <- xlim[1] / 10
    polygon(c(xs,      rev(xs)), 
            #c(param-p.unc, rev(param+p.unc)),
            c(p.l, rev(p.u)),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)
    lines(xs, param, lwd=2, col=col.pal[ii], lty=1)
    
    if (ii < 4) {
        param <- scaling_unc[,2*ii]   * 100
        p.unc <- scaling_unc[,2*ii+1] * 100
        param[param<ylim[1]] <- ylim[1]/10
        xs <- scaling_unc[,1]
        xs[xs<xlim[1]] <- xlim[1] / 10
        #if (use.spline) {
        #    fit <- predict(smooth.spline(xs, param, spar=0), interp)$y
        #    lines(interp, fit, lwd=2, col=1, lty=2)
        #} else {
            lines(xs, param, lwd=2, col=1, lty=2)
        #}
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]-10, xlim[2]*1.1, ylim[2]+10, col='white', border=NA)
    #rect(xlim[1], ylim[2], xlim[2]*1.1, ylim[2]*1.1, col='white', border=NA)
    par(xpd=F)
    
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.y) axis(2, yticks[c(T,F)], paste0(yticks[c(T,F)], '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    if (make.x) axis(1, c(0.1, 1, 10, 100), c('0.1', '1', '10', '100'),
        tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    rect(10**(sum(log10(xlim))/2*.5), ylim[2]*0.72, 
         10**(sum(log10(xlim))/2*1.5), ylim[2]*0.9, 
        col='white', border=1, lwd=1.5)
    text(10**(sum(log10(xlim))/2), ylim[2]*0.9, 
        labels=lab, pos=1, cex=0.9*text.cex)
    par(family=font, xpd=F)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=F, las=0, cex=text.cex)
    if (make.x) mtext(xlab, 1, 2, outer=F, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_unc, paste0('imp-teff', ii), 
    DF=unc, ii=ii,
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    ylab="Relative Uncertainty", 
    xlab="Added Uncertainty [K]",
    lab=params[ii],
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}

