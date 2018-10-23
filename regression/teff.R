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

dirs <- list.files('teff', full.names=T)
baseline <- file.path('feh', 'learn-feh')
bias_dirs <- dirs[grep('bias', dirs)]
imp_dirs <- dirs[grep('imp', dirs)]

biases <- as.numeric(sub('learn-teff_bias', '', basename(bias_dirs)))
imps <- as.numeric(sub('learn-teff_imp', '', basename(imp_dirs)))

col.pal <- c("#DB4D48", "#F29559", blue, "#323031")
params <- c('Radius', 'Mass', 'Density', 'Age')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return (NULL)
    
    #perturb <- read.table(file.path('perturb', 'feh', 
    #    sub('\\.', '_perturb.', basename(filename))), header=1)
    #if (any(perturb$Fe.H > 0.44)) return(NULL)
    
    covs <- read.table(filename, header=T)
    n <- nrow(covs)
    covs <- covs[complete.cases(covs),]
    if (nrow(covs) < n/2) return(NULL) 
    #covs <- covs[1:(min(nrow(covs), n)),]
    cov.DF <- with(covs, 
        data.frame(KIC=KIC,
                   R=median(radius),    e_R=mad(radius),
                   M=median(M),         e_M=mad(M),
                   #alpha=median(alpha), e_alpha=mad(alpha),
                   rho=median(density), e_rho=mad(density),
                   age=median(age),     e_age=mad(age)
        )
    )
    
    perturb <- read.table(file.path('perturb', 'feh', 
        sub('\\.', '_perturb.', basename(filename))), header=1)
    perturb.DF <- with(perturb,
            data.frame(Teff=median(Teff), e_Teff=mad(Teff),
                       nu_max=median(nu_max), e_nu_max=mad(nu_max),
                       FeH=median(Fe.H),  e_FeH=mad(Fe.H)))
    
    if ('Dnu0' %in% names(perturb)) 
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(Dnu=median(Dnu0),  e_Dnu=mad(Dnu0))))
    
    cbind(cov.DF, perturb.DF)
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    parallelMap(get_star, list.files(directory, recursive=T, full.names=T)))

baseline.DF <- get_DF(file.path(baseline))

biases.list <- parallelMap(get_DF, bias_dirs)
imps.list <- parallelMap(get_DF, imp_dirs)

unc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    if (ii == 0) {
        imp <- 0
        imp.DF <- baseline.DF
    } else {
        imp <- imps[ii]
        imp.DF <- imps.list[[ii]]
    }
    if (nrow(imp.DF) < 50) return(NULL)
    with(imp.DF, data.frame(imp=imp,
        R=mean(e_R),             e_R=sd(e_R),
        M=mean(e_M),             e_M=sd(e_M),
        #alpha=mean(e_alpha), e_alpha=sd(e_alpha),
        rho=mean(e_rho),       e_rho=sd(e_rho),
        age=mean(e_age),       e_age=sd(e_age)))
}, ii=0:length(imps)))
unc <- unc[order(unc$imp),]
baseline.uncR <- signif(unc[1,]$R, 1)
baseline.uncM <- signif(unc[1,]$M, 1)
baseline.uncage <- signif(unc[1,]$age, 1)

diff.1 <- unc[which(unc$imp == 100),]-unc[1,]
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
    if (nrow(imp.DF) < 50) return(NULL)
    with(imp.DF, data.frame(imp=imp,
        R=mean(e_R/R),             e_R=sd(e_R/R),
        #R_u=quantile(e_R/R, 0.75),   R_l=quantile(e_R/R, 0.25),
        M=mean(e_M/M),             e_M=sd(e_M/M),
        #M_u=quantile(e_M/M, 0.75),   M_l=quantile(e_M/M, 0.25),
        #alpha=mean(e_alpha/alpha), e_alpha=sd(e_alpha/alpha),
        rho=mean(e_rho/rho),     e_rho=sd(e_rho/rho),
        #rho_u=quantile(e_rho/rho, 0.75),   rho_l=quantile(e_rho/rho, 0.25),
        age=mean(e_age/age),       e_age=sd(e_age/age),
        #age_u=quantile(e_age/age, 0.75),   age_l=quantile(e_age/age, 0.25),
        N=nrow(imp.DF)))
}, ii=0:length(imps)))
unc <- unc[order(unc$imp),]
#abs(unc[7,]-unc[1,])/unc[1,]*100

diff.1 <- abs(unc[which(unc$imp == 100),]-unc[1,]) / unc[1,] * 100
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
bias.1 <- acc[which(acc$bias == 100),]*100
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

cat(paste0("\n\nSimilarly, a systematic error of 100~K in effective temperature results in relative differences of ", bias.age, "\\%, ", bias.M, "\\%, and ", bias.R, "\\%, respectively. \n\n Similarly, increasing the reported uncertainty of $T_\\text{eff}$ by 100~K increases those uncertainties by ", diff.age, "~Gyr (", rel.age, "\\%), ", diff.M, "~$\\text{M}_\\odot$ (", rel.M, "\\%), and ", diff.R, "~$\\text{R}_\\odot$ (", rel.R, "\\%).\n\n"))

#cat(paste0("\n\nWe find that a systematic error of 100~K in effective temperature measurements translates on average to differences of ", bias.age, "\\%, ", bias.M, "\\%, and ", bias.R, "\\% in the resulting stellar ages, masses, and radii, respectively, which are smaller than the reported relative uncertainties for these quantities (10\\%, 3\\%, 1\\%). \n We furthermore find that increasing the reported uncertainty of [Fe/H] measurements by 100~K (100\\%) increases the uncertainties of stellar ages, masses, and radii on average by only ", diff.age, "~Gyr (", rel.age, "\\%)", ", ", diff.M, "~$\\text{M}_\\odot$ (", rel.M, "\\%), and ", diff.R, "~$\\text{R}_\\odot$ (", rel.R, "\\%)", ", respectively, which are well below the reported uncertainties for these estimates (", baseline.uncage, "~Gyr", ", ", baseline.uncM, "~$\\text{M}_\\odot$, ", baseline.uncR, "~$\\text{R}_\\odot$", ").\n\n"))



### Scaling relations 
get_scaling <- function(DF, 
        solar_Dnu = 135.1, solar_e_Dnu = 0.1,
        solar_Teff = 5772,    solar_e_Teff = 0.8,
        solar_nu_max = 3090,  solar_e_nu_max = 30,
        n_trials=1000) {
    
    scaling_R <- with(DF, 
        (nu_max/solar_nu_max) * 
        (Dnu/solar_Dnu)**-2 * 
        (Teff / solar_Teff)**0.5)
    
    scaling_M <- with(DF, 
        (nu_max / solar_nu_max)**3 * 
           (Dnu / solar_Dnu)**-4 * 
          (Teff / solar_Teff)**(3/2))
    
    scaling_rho <- with(DF,
           (Dnu / solar_Dnu)**2 * 
          (Teff / solar_Teff))
    #scaling_rho <- scaling_M / scaling_R**3
    
    set.seed(0)
    scaling_dR <- do.call(rbind, Map(function(ii) with(DF, 
            (rnorm(nrow(DF), nu_max, e_nu_max) / 
                rnorm(1, solar_nu_max, solar_e_nu_max)) * 
            (rnorm(nrow(DF), Dnu, e_Dnu) / 
                rnorm(1, solar_Dnu, solar_e_Dnu))**-2 * 
            (rnorm(nrow(DF), Teff, e_Teff) / 
                rnorm(1, solar_Teff, solar_e_Teff))**0.5), 
        ii=1:1000))
    set.seed(0)
    scaling_dM <- do.call(rbind, Map(function(ii) with(DF, 
            (rnorm(nrow(DF), nu_max, e_nu_max) / 
                rnorm(1, solar_nu_max, solar_e_nu_max))**3 * 
            (rnorm(nrow(DF), Dnu, e_Dnu) / 
                rnorm(1, solar_Dnu, solar_e_Dnu))**-4 * 
            (rnorm(nrow(DF), Teff, e_Teff) /
                rnorm(1, solar_Teff, solar_e_Teff))**(3/2)), 
        ii=1:1000))
    #scaling_drho <- scaling_dM / scaling_dR**3
    set.seed(0)
    scaling_drho <- do.call(rbind, Map(function(ii) with(DF, 
            (rnorm(nrow(DF), nu_max, e_nu_max) / 
                rnorm(1, solar_nu_max, solar_e_nu_max))**0 * 
            (rnorm(nrow(DF), Dnu, e_Dnu) / 
                rnorm(1, solar_Dnu, solar_e_Dnu))**2 * 
            (rnorm(nrow(DF), Teff, e_Teff) /
                rnorm(1, solar_Teff, solar_e_Teff))), 
        ii=1:1000))
    
    scaling_dR   <- apply(scaling_dR,   2, function(x) mad(x, na.rm=T))
    scaling_dM   <- apply(scaling_dM,   2, function(x) mad(x, na.rm=T))
    scaling_drho <- apply(scaling_drho, 2, function(x) mad(x, na.rm=T))
    
    data.frame(KIC=DF$KIC,
               R=scaling_R,     e_R=scaling_dR, 
               M=scaling_M,     e_M=scaling_dM, 
               rho=scaling_rho, e_rho=scaling_drho)
}

scaling_unc <- do.call(plyr:::rbind.fill, Map(function(ii) { 
    imp <- if (ii == 0) 0 else imps[ii] 
    imp.DF <- baseline.DF 
    imp.DF$e_Teff <- imp.DF$e_Teff + imp 
    scaling.DF <- get_scaling(imp.DF)
    if (nrow(imp.DF) < 50) return(NULL)
    with(scaling.DF, data.frame(imp=imp,
        R=median(e_R/R),             e_R=mad(e_R/R),
        M=median(e_M/M),             e_M=mad(e_M/M),
        rho=median(e_rho/rho),     e_rho=mad(e_rho/rho)))
}, ii=0:length(imps)))
scaling_unc <- scaling_unc[order(scaling_unc$imp),]

baseline.scaling <- get_scaling(baseline.DF)
scaling_acc <- do.call(plyr:::rbind.fill, Map(function(ii) {
    bias <- if (ii == 0) 0 else biases[ii]
    bias.DF <- baseline.DF
    bias.DF$Teff <- bias.DF$Teff + bias
    DF <- merge(baseline.scaling, get_scaling(bias.DF), by='KIC')
    with(DF, data.frame(bias=bias,
        R=median(abs(R.x - R.y)/R.x),
        d_R=mad(abs(R.x - R.y)/R.x),
        M=median(abs(M.x - M.y)/M.x),
        d_M=mad(abs(M.x - M.y)/M.x),
        rho=median(abs(rho.x - rho.y)/rho.y),
        d_rho=mad(abs(rho.x - rho.y)/rho.y)))
}, ii=0:length(biases)))
scaling_acc <- scaling_acc[order(scaling_acc$bias),]






plot_unc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.3, -0.3, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(0.1, 100)
    ylim <- c(0, 21)
    
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
    rect(xlim[2], ylim[1], xlim[2]*1.05, ylim[2], col='white', border=NA)
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
    rect(10**(sum(log10(xlim))/2*.16), ylim[2]*0.72, 
         10**(sum(log10(xlim))/2*1.84), ylim[2]*0.9, 
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









plot_acc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(matrix(1:2, ncol=2), width=1, heights=c(2,2), respect=FALSE)
    
    par(oma=mar+c(0.3, -0.3, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(100, 0.01)
    ylim <- c(0, 4)
    
    xticks <- 10**pretty(log10(xlim))
    yticks <- pretty(ylim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, 
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
    
    interp <- 10**(seq(log10(min(DF[,1][DF[,1]>0])/10), 
                       log10(max(DF[,1])), 
                    length.out=100))
    interp <- c(-rev(interp), 0, interp)
    #    fit <- predict(smooth.spline(xs, param, spar=0), interp)$y
    #unc.fit <- predict(smooth.spline(xs, p.unc, spar=1), interp)$y
    
    xs[abs(xs) <= min(xlim)] <- min(xlim)#10
    
    polygon(c(abs(xs[left]),        rev(abs(xs[left]))), 
            #c(c(param-p.unc)[left], rev(c(param+p.unc)[left])),
            c(p.l[left], rev(p.u[left])),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)
    lines(abs(xs[left]), param[left], lwd=2, col=col.pal[ii], lty=1)
    
    if (ii < 4) {
        param <- scaling_acc[,2*ii]   * 100
        p.unc <- scaling_acc[,2*ii+1] * 100
        param[param<=ylim[1]] <- 0
        xs <- scaling_acc[,1]
        left  <- xs<=0
        right <- xs>=0
        xs[abs(xs) <= min(xlim)] <- min(xlim)
        
        #if (use.spline || T) {
        fit <- predict(smooth.spline(abs(xs[left]), param[left], spar=0), 
            interp)$y
        lines(interp, fit, lwd=2, col=1, lty=2)
        #} else {
        #    lines(xs[left], param[left], lwd=2, col=1, lty=2)
        #}
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
    
    if (make.x) 
        axis(1, c(100, 1), c('-100', '-1'), 
        tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    if (make.y) axis(2, yticks, paste0(yticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=TRUE, las=0, cex=text.cex)
    
    xlim <- rev(xlim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, 
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
    
    polygon(c(abs(xs[right]),      rev(abs(xs[right]))), 
            #c(c(param-p.unc)[right], rev(c(param+p.unc)[right])),
            c(p.l[right], rev(p.u[right])),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)
    lines(xs[right], param[right], lwd=2, col=col.pal[ii], lty=1)
    
    if (ii < 4) {
        param <- scaling_acc[,2*ii]   * 100
        p.unc <- scaling_acc[,2*ii+1] * 100
        param[param<=ylim[1]] <- 0
        xs <- scaling_acc[,1]
        xs[abs(xs) <= min(xlim)] <- min(xlim)
        #if (use.spline || T) {
        fit <- predict(smooth.spline(abs(xs[right]), param[right], spar=0), 
            interp)$y
        lines(interp, fit, lwd=2, col=1, lty=2)
        #} else {
        #    lines(xs[right], param[right], lwd=2, col=1, lty=2)
        #}
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]-0.1, 
         xlim[2]*10, ylim[2]+0.1, 
         col='white', border=NA)
    rect(xlim[1],  ylim[1], 
         xlim[2]*10,  ylim[1]-0.1, 
         col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) 
        axis(1, c(100, 1), 
                c('100', '1'), 
            tick=F, cex.axis=text.cex, tcl=0, las=1, 
            mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    rect(xlim[1]/13.3, ylim[2]*0.72, xlim[1]*13.3, ylim[2]*0.9, 
        col='white', border=1, lwd=1.5)
    text(xlim[1], ylim[2]*0.9, labels=lab, pos=1, cex=0.9*text.cex)
    par(family=font, xpd=F)
    
    if (make.x) mtext(xlab, 1, 2, outer=TRUE, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_acc, paste0('bias-teff', ii), 
    DF=acc, ii=ii, 
    ylab="Relative Difference",
    xlab="Systematic Error [K]",
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
    rect(10**(sum(log10(xlim))/2*.16), ylim[2]*0.72, 
         10**(sum(log10(xlim))/2*1.84), ylim[2]*0.9, 
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









plot_acc <- function(DF, ii, lab, xlab, ylab, ..., make.x=T, make.y=T,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(matrix(1:2, ncol=2), width=1, heights=c(2,2), respect=FALSE)
    
    par(oma=mar+c(0.3, -0.3, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    xlim <- c(500, 0.1)
    ylim <- c(0, 21)
    
    xticks <- 10**pretty(log10(xlim))
    yticks <- pretty(ylim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, 
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
    
    interp <- 10**(seq(log10(min(DF[,1][DF[,1]>0])/10), 
                       log10(max(DF[,1])), 
                    length.out=100))
    interp <- c(-rev(interp), 0, interp)
    #    fit <- predict(smooth.spline(xs, param, spar=0), interp)$y
    #unc.fit <- predict(smooth.spline(xs, p.unc, spar=1), interp)$y
    
    xs[abs(xs) <= min(xlim)] <- min(xlim)#10
    
    polygon(c(abs(xs[left]),        rev(abs(xs[left]))), 
            #c(c(param-p.unc)[left], rev(c(param+p.unc)[left])),
            c(p.l[left], rev(p.u[left])),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)
    lines(abs(xs[left]), param[left], lwd=2, col=col.pal[ii], lty=1)
    
    if (ii < 4) {
        param <- scaling_acc[,2*ii]   * 100
        p.unc <- scaling_acc[,2*ii+1] * 100
        param[param<=ylim[1]] <- 0
        xs <- scaling_acc[,1]
        left  <- xs<=0
        right <- xs>=0
        xs[abs(xs) <= min(xlim)] <- min(xlim)
        
        #if (use.spline || T) {
        fit <- predict(smooth.spline(abs(xs[left]), param[left], spar=0), 
            interp)$y
        lines(interp, fit, lwd=2, col=1, lty=2)
        #} else {
        #    lines(xs[left], param[left], lwd=2, col=1, lty=2)
        #}
    }
    
    par(xpd=NA)
    rect(xlim[1], ylim[1], 
         xlim[2], ylim[1]-10, 
         col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, #minorn=10,
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=text.cex, #minorn=5,
        family=font, las=1, majorn=5, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) {
        axis(1, c(100, 10, 1), c('-100', '-10', '-1'), 
            tick=F, 
            cex.axis=text.cex, tcl=0, las=1, 
            mgp=mgp+c(0, 0.25, 0), lwd=1.5)
        axis(1, c(10), c('-10'), 
            tick=F, 
            cex.axis=text.cex, tcl=0, las=1, 
            mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    }
    
    if (make.y) axis(2, yticks, paste0(yticks, '%'), tick=F, 
        cex.axis=text.cex, tcl=0, las=1, #tcl, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    if (make.y) mtext(ylab, 2, 2.5, outer=TRUE, las=0, cex=text.cex)
    
    xlim <- rev(xlim)
    
    plot(NA, axes=F, log='x',
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, 
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
    
    polygon(c(abs(xs[right]),      rev(abs(xs[right]))), 
            #c(c(param-p.unc)[right], rev(c(param+p.unc)[right])),
            c(p.l[right], rev(p.u[right])),
            col=adjustcolor(col.pal[ii], alpha.f=0.5), 
            border=NA)
    lines(xs[right], param[right], lwd=2, col=col.pal[ii], lty=1)
    
    if (ii < 4) {
        param <- scaling_acc[,2*ii]   * 100
        p.unc <- scaling_acc[,2*ii+1] * 100
        param[param<=ylim[1]] <- 0
        xs <- scaling_acc[,1]
        xs[abs(xs) <= min(xlim)] <- min(xlim)
        #if (use.spline || T) {
        fit <- predict(smooth.spline(abs(xs[right]), param[right], spar=0), 
            interp)$y
        lines(interp, fit, lwd=2, col=1, lty=2)
        #} else {
        #    lines(xs[right], param[right], lwd=2, col=1, lty=2)
        #}
    }
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]-0.1, 
         xlim[2]*10, ylim[2]+10, 
         col='white', border=NA)
    rect(xlim[1],  ylim[1], 
         xlim[2]*10,  ylim[1]-10, 
         col='white', border=NA)
    par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F, lwd.ticks=par()$lwd)
    
    if (make.x) 
        axis(1, c(100, 10, 1), 
                c('100', '10', '1'), 
            tick=F, cex.axis=text.cex, tcl=0, las=1, 
            mgp=mgp+c(0, 0.25, 0), lwd=1.5)
    
    box(lwd=par()$lwd)
    
    par(family="Helvetica", xpd=NA)
    rect(xlim[1]/6.85, ylim[2]*0.72, 
         xlim[1]*6.85, ylim[2]*0.9, 
        col='white', border=1, lwd=1.5)
    text(xlim[1], ylim[2]*0.9, labels=lab, pos=1, cex=0.9*text.cex)
    par(family=font, xpd=F)
    
    if (make.x) mtext(xlab, 1, 2, outer=TRUE, cex=text.cex)
}

for (ii in 1:length(col.pal)) {
make_plots(plot_acc, paste0('bias-teff', ii), 
    DF=acc, ii=ii, 
    ylab="Relative Difference",
    xlab="Systematic Error [K]",
    lab=params[ii],
    make.x=ii >= 3, 
    make.y=ii == 1 || ii == 3, 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    cex.paper=0.95) 
}


