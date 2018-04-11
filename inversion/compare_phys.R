#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

### CONSTANTS 
M_sun   = 1.988475e33 # g
R_sun   = 6.957e10 # cm

k.pair  = u_Y
rs      = seq(0.05, 0.3, 0.05) 
sampler = c(T)
models  = get_model_list() 
kern.interp.xs = seq(0, 1, 0.001)

num_knots = 60
degree = 3

targ.kern.type <- 'mod_Gauss'

phys <- c('inv', 'inv-ov', 'inv-D', 'inv-Dov')
star.names <- c('3656476', '5184732', '6225718', '8760414', 
    '12069424', '12069449', '4914923', '6116048', '3427720', '5774694', 
    '6106415', '7680114', '8006161', '8179536', '9098294', 
    '10516096', '10963065', '12009504',
    '1435467', '9139163', '8938364', '10454113',
    '11253226')
save.files <- list.files('save')
col.pal <- c('darkgray', orange, blue, red)
col.pal2 <- brewer.pal(11, "Spectral")[c(1:4,8:11)][c(1:4, 7:8)]

#rnorms <- sapply(1:1024, function(x) rnorm(num_knots, 0, 0.1))
#rnorms <- sapply(1:1024, function(x) runif(num_knots, -0.25, 0.25))
rnorms <- sapply(1:1024, function(x) 
    runif(num_knots, -0.3/2, 0.3/2) + runif(1, -0.25/2, 0.25/2))

sunstar <- read.table(file.path('..', 'grid', 'ref_mods', 
        'inv', '5774694', 'LOGS_MS', 
        'profile1-freqs', 'profile1.data.FGONG.dat'),
    header=1)
sunstar.u <- with(sunstar, splinefun(x, u)(kern.interp.xs))
calibrated <- read.table(file.path('..', 'grid', 'calibrate', 
        'LOGS_MS', 'profile1-freqs', 'profile1.data.FGONG.dat'),
    header=1)
calibrated.u <- with(calibrated, splinefun(x, u)(kern.interp.xs))

plot_du.u <- function(du.u, inversion, sampler, star.name, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    par(mar=mar+c(0.2,-0.2,-0.4,-0.5), mgp=mgp+c(0, 0.4, 0))
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0, 0.35), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lty=2, col='black')
    
    intgrl <- function(ii) sintegral(kern.interp.xs, 
        du.u*inversion$avg_kerns[,ii])$value
    ress <- sapply(1:nrow(inversion$result), intgrl)
    resids <- ress - sapply(1:nrow(inversion$result), 
        function(ii) du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
    lines(kern.interp.xs, du.u, type='l', lwd=1.5)
    points(inversion$result$fwhm.mid[sampler], ress[sampler],
        col=red, pch=20, cex=1)
    
    legend("bottomright", bty="n", cex=0.8*text.cex,
        legend=star.name, inset=c(-0.01, -0.05))
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.6, 0, 0))
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_du.u.prof <- function(du.u.prof, inversion, sampler, star.name, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    par(mar=mar+c(0.2,-0.2,-0.4,-0.5), mgp=mgp+c(0, 0.4, 0))
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         #xlim=c(0, 0.35), 
         xlim=c(0, 1),
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lty=2, col='black')
    
    for (ii in 1:nrow(du.u.prof)) {
        lines(kern.interp.xs, du.u.prof[ii,], type='l', lwd=0.1, 
            col=adjustcolor(red, 0.2))
    }
    #intgrl <- function(ii) sintegral(kern.interp.xs, 
    #    du.u*inversion$avg_kerns[,ii])$value
    #ress <- sapply(1:nrow(inversion$result), intgrl)
    #resids <- ress - sapply(1:nrow(inversion$result), 
    #    function(ii) du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
    #lines(kern.interp.xs, du.u, type='l', lwd=1.5)
    #points(inversion$result$fwhm.mid[sampler], ress[sampler],
    #    col=red, pch=20, cex=1)
    
    #legend("bottomright", bty="n", cex=0.8*text.cex,
    #    legend=star.name, inset=c(-0.01, -0.05))
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.6, 0, 0))
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_residuals <- function(residuals, means, errs, inversion, sampler, 
        star.name, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    par(mar=mar+c(0.2,-0.2,-0.4,-0.5), mgp=mgp+c(0, 0.4, 0))
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0, 0.35), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lty=2, col='black')
    for (ii in 1:nrow(residuals)) {
        with(inversion$result[sampler,],
            points(rnorm(length(rs), fwhm.mid, (fwhm.right-fwhm.left)/14),
            residuals[ii,sampler], 
            pch=2, cex=0.01, col='gray'))
    }
    
    with(inversion$result, 
        arrows(fwhm.mid-0.003, err, fwhm.mid-0.003, -err, 
            code=3, angle=90, length=0.01, lwd=1.5, lty=1, col=1))
    with(inversion$result, 
        arrows(fwhm.mid+0.003, means+errs, fwhm.mid+0.003, means-errs, 
            code=3, angle=90, length=0.01, lwd=1.5, lty=1, col=red))
    
    legend('topleft', lty=NA, lwd=1, pch=20, col=c(1, red), bty='n',
        cex=0.8*text.cex, x.intersp=0.01, inset=c(-0.02, -0.04),
        legend=c('Formal error bars', 'Monte Carlo'))
    legend("bottomright", bty="n", cex=0.8*text.cex,
        legend=star.name, inset=c(-0.01, -0.05))
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.6, 0, 0))
    title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    title(ylab=expression("Residuals"))
}

plot_residuals2 <- function(du.u.avg, du.u.true, inversion, sampler, 
        star.name, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times") {
    par(mar=mar+c(1.2,-0.2,-0.4,-0.5), mgp=mgp+c(0, 0.4, 0), pty="s")
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(-0.25, 0.25), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    #abline(h=0, lty=2, col='black')
    segments(-1,-1, 1,1, lty=2)
    for (ii in 1:nrow(inversion$result)) {
        if (!sampler[ii]) next 
        points(du.u.true[,ii], du.u.avg[,ii],
            pch=2, cex=0.2, col=adjustcolor(col.pal2[ii], alpha.f=0.5))
    }
    
    legend.txt <- c()
    for (ii in nrow(inversion$result):1) {
        if (!sampler[ii]) next 
        x <- du.u.true[,ii]
        y <- du.u.avg[,ii]
        fit.lm <- lm(y~x)
        new.x <- c(-0.18, 0.18)
        new.y <- predict(fit.lm, newdata=data.frame(x=new.x))
        #arrows(new.x[1], new.y[1], new.x[2], new.y[2],
        #    lwd=5, length=0.05, code=3, col=1)
        arrows(new.x[1], new.y[1], new.x[2], new.y[2],
            lwd=4, length=0.05, code=3, col='white')
        arrows(new.x[1], new.y[1], new.x[2], new.y[2],
            lwd=3, length=0.05, code=3, col=col.pal2[ii])
        #lines(new.x, ,
        #    lwd=2, col=col.pal2[ii])
        
        legend.txt <- c(legend.txt, bquote(r[0]==.(inversion$result$rs[ii])))
        #value <- signif(summary(fit.lm)$coef[2,][1], 3)
        #while (nchar(value) < 5) value <- paste0(value, 0) 
        #legend.txt <- c(legend.txt, 
        #    paste0(value, " ± ", 
        #           signif(summary(fit.lm)$coef[2,][2], 2)))
    }
    
    legend("topleft", bty="n", pch=20, cex=0.8*text.cex, inset=c(-0.01, -0.025),
        legend=as.expression(legend.txt), col=col.pal2[sampler], x.intersp=0.66)
    legend("bottomright", bty="n", cex=0.8*text.cex,
        legend=star.name, inset=c(-0.01, -0.03))
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=F)
    magaxis(2, tcl=-0.25, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    axis(1, at=c(-0.2, 0.2), labels=T, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(1.7, 0, 0))
    #title(xlab=expression(integral(K^"avg"~frac(du,u)~dr)))
    title(xlab=expression(frac(du,u)(r[0])))
    par(mgp=mgp+c(0.8, 0, 0))
    title(ylab=expression("<"*du/u*">"[r[0]]))
}

S <- function(x, alpha, ks, degree, num_knots) {
    sapply(x, function(x) sum(sapply( 1:(num_knots-(degree+1)), 
        function(ii) alpha[ii] * B(x, ii, degree, ks) )))
}


check_kerns <- function(star.name, inversion, ref.mod, sampler) {
    ## calculate acoustic depth, Aerts et al eq 3.228 
    int_0.r <- function(x, y) as.numeric(cumtrapz(x, y)) # int_0^r 
    int_r.R <- function(x, y) trapz(x, y) - int_0.r(x, y) # int_r^R 
    tau <- int_r.R(kern.interp.xs*R_sun*ref.mod$R, 
        1/ref.mod$cs.spl(kern.interp.xs))
    tau[length(tau)] <- 0 #tau[length(tau)-1]
    ks <- get_knots(tau, num_knots, degree)
    ks <- splinefun(tau, kern.interp.xs)(rev(ks))
    ks[length(ks)] <- 1.001
    
    #du.u <- S(kern.interp.xs, rnorm(num_knots, 0, 0.1), ks, degree, num_knots)
    #du.u <- S(kern.interp.xs, rnorms[,1], ks, degree, num_knots)
    #du.u <- 0.1*kern.interp.xs+0.1
    #du.u <- (sunstar.u - calibrated.u)/sunstar.u
    if (F) {
    plot(kern.interp.xs, du.u, type='l', lwd=1.5,
        xlim=c(0, 0.35), ylim=c(-0.25, 0.25), xaxs='i', yaxs='i',
        xlab="Radius", ylab=expression(delta*u/u))
    intgrl <- function(ii) sintegral(kern.interp.xs, 
        du.u*inversion$avg_kerns[,ii])$value
    ress <- sapply(1:nrow(inversion$result), intgrl)
    resids <- ress - sapply(1:nrow(inversion$result), 
        function(ii) du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
    points(inversion$result$fwhm.mid, ress,
           col=red, pch=20, cex=1)
    dev.off()
    }
    
    du.u <- (sunstar.u - calibrated.u)/sunstar.u
    make_plots(plot_du.u, paste0(star.name, '-true_diff-solar'), #'true_diff', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, 
        cex.paper=1.16, paper_pdf_height=6.75206*2/3, star.name=star.name, 
        du.u=du.u, inversion=inversion, sampler=sampler)
    
    du.u <- S(kern.interp.xs, rnorms[,8], ks, degree, num_knots)
    make_plots(plot_du.u, paste0(star.name, '-true_diff-Bspline'), #'true_diff', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, 
        cex.paper=1.16, paper_pdf_height=6.75206*2/3, star.name=star.name, 
        du.u=du.u, inversion=inversion, sampler=sampler)
    
    
    du.u.prof <- do.call(rbind, parallelMap(function(ii) {
        S(kern.interp.xs, rnorms[,ii], ks, degree, num_knots)
    }, ii=1:ncol(rnorms)))#128))#
    
    du.u.avg <- do.call(rbind, parallelMap(function(ii) {
        #du.u <- S(kern.interp.xs, 
        #    rnorms[,ii], ks, degree, num_knots)
        du.u <- du.u.prof[ii,]
        intgrl <- function(ii) sintegral(kern.interp.xs, 
            du.u*inversion$avg_kerns[,ii])$value
        sapply(1:nrow(inversion$result), intgrl)
    }, ii=1:nrow(du.u.prof)))
    
    du.u.true <- do.call(rbind, Map(function(ii) {
        du.u <- du.u.prof[ii,]
        sapply(1:nrow(inversion$result), 
            function(ii) du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
    }, ii=1:nrow(du.u.prof)))
    
    residuals <- du.u.avg - du.u.true 
    errs <- apply(residuals, 2, std)
    means <- apply(residuals, 2, mean)
    
    make_plots(plot_residuals, paste0(star.name, '-true_diff-residuals'), 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, 
        cex.paper=1.16, paper_pdf_height=6.75206*2/3, sampler=sampler,
        star.name=star.name, 
        residuals=residuals, means=means, errs=errs, inversion=inversion)
    
    make_plots(plot_residuals2, paste0(star.name, '-true_diff-residuals2'), 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, 
        cex.paper=1.16, paper_pdf_height=6.97522, sampler=sampler,
        star.name=star.name, 
        du.u.avg=du.u.avg, du.u.true=du.u.true, inversion=inversion)
    
    make_plots(plot_du.u.prof, 'duprof', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, 
        cex.paper=1.16, paper_pdf_height=6.75206*2/3, du.u.prof=du.u.prof,
        star.name=star.name, sampler=sampler, inversion=inversion)
}






for (star.name_i in 1:length(star.names)) {
    
    star.name <- star.names[star.name_i]
    
    rm.first.kernel <- star.name %in% c('12009504', '8179536', '1435467', 
        '2837475', '9812850', '3427720', '8938364', '9139163', '10454113',
        '10516096', '11253226')
    rm.second.kernel <- star.name %in% c('12009504', '2837475', '8179536',
        '9139163', '9812850', '1435467', '3427720', '10454113', '11253226')
    rm.third.kernel <- star.name %in% c()
    rm.fourth.kernel <- star.name %in% c('9139163', '8938364', '10454113')
    rm.fifth.kernel <- star.name %in% c('1435467', '9139163', '10454113')
    rm.last.kernel <- star.name %in% c('3427720', '5774694', '6106415', 
        '7680114', '8006161', '8179536', '9098294', '10963065', '12009504',
        '1435467', '6603624', '9139163', '8938364', '9812850', '10454113')
    
    for (filename in save.files[grep(star.name, save.files)]) {
        load(file.path('save', filename)) # loads ref.mod etc 
    }
    
    phy_res <- list()
    phy_mass <- list()
    phy_radius <- list()
    
    for (phy_i in 1:length(phys)) {
        phy <- phys[phy_i]
        model.filename <- file.path('..', 'grid', 'ref_mods', phy, star.name, 
            'LOGS_MS', 'profile1-freqs', 'profile1.data.FGONG.dat')
        if (!file.exists(model.filename)) next 
        model.prof <- read.table(model.filename, header=1)
        model <- list(k.pair=k.pair, 
            f1.spl=with(model.prof, splinefun(x, u)),
            f1.exp=k.pair$f1.exp, 
            f2.exp=k.pair$f2.exp)
        
        inversion <- lists_to_inversion(model=model, rs=rs, 
            inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
            cross.lists=cross.lists, inv.params=inv.params, 
            kern.interp.xs=kern.interp.xs) 
        
        mass <- max(model.prof$m)/M_sun
        mass <- signif(mass,3)
        
        radius <- model.prof$r[which.min(abs(model.prof$x-1))]/R_sun
        radius <- signif(radius,3)
        
        penalties <- sapply(1:nrow(inversion$result), function(jj) {
            avg_kern <- inversion$avg_kerns[,jj]
            cross_kern <- inversion$cross_kerns[,jj]
            max(abs(avg_kern[kern.interp.xs > 0.5]), 
                abs(cross_kern))
        })
        
        sampler <- penalties < 3
        if (rm.first.kernel) sampler[1] <- F
        if (rm.second.kernel) sampler[2] <- F
        if (rm.third.kernel) sampler[3] <- F
        if (rm.fourth.kernel) sampler[4] <- F
        if (rm.fifth.kernel) sampler[5] <- F
        if (rm.last.kernel) sampler[length(sampler)] <- F 
        
        if (!any(sampler)) next 
        
        phy_res[[phy]] <- inversion$result[sampler,]
        phy_mass[[phy]] <- mass
        phy_radius[[phy]] <- radius
        
        if (phy_i == 1) {
            check_kerns(star.name, inversion, ref.mod, sampler) 
            
            make_plots_inversion_all(model, inversion, 
                kern.interp.xs=kern.interp.xs,
                k.str=paste0('-', star.name, '-', phy), 
                cross.inset="bottomright",
                inversion_ylim = c(-0.4, 0.4),
                xlim=c(0, 0.4),
                caption=star.name, 
                col.pal=col.pal2[sampler],#col.pal[phy_i], 
                caption2=as.expression(bquote('M' == .(mass))),
                sampler=sampler,
                slides=F, wide=F, tall=F, 
                suppress_cross=phy_i>1,
                suppress_avg=phy_i>1,
                filepath=file.path('plots', 'phys'),
                #core.bound=core.bound, 
                cross_kern_ylim=c(-0.8, 0.3))
        }
    }
    
    #caption2 <- sapply(1:length(phys), function(ii) {
    #    paste0("M=", phy_mass[[ii]], ", R=", phy_radius[[ii]], 
    #        if (ii==2) " (ov)" else if (ii==3) " (D)" else if (ii==4) " (D, ov)" else "")
    #})
    caption2 <- c('No diffusion, no overshoot',
        'No diffusion, overshoot', 
        'Diffusion, no overshoot',
        'Diffusion, overshoot')
    
    make_plots(plot_inversion_multi, star.name, 
        filepath=file.path('plots', 'phys-comp'),
        slides=F, wide=F, tall=F, #make_png=F, 
        inversion.list=phy_res,
        kern.interp.xs=kern.interp.xs,
        #ylim = c(-0.4, 0.4),
        ylim = c(-0.25, 0.25),
        xlim=c(0, 0.35),
        caption=star.name, 
        #caption2=if (star.name == '3427720' || star.name == '8179536') 
        #    caption2 else NULL,
        caption2=NULL, 
        col.pal=col.pal,
        cex.paper=1.16,
        paper_pdf_height=6.75206*2/3)
        #caption2=as.expression(bquote('M' == .(mass))))
}



