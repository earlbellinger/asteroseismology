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

num_knots = 20
degree = 3
n_trials = 1024
dir.create('save-poly', showWarnings = FALSE)


cex.paper <- 1.26#3
cex.sun <- 0.96
paper_pdf_height <- 5.35#6.75206*2/3*6/5
paper_pdf_width <- 4.17309*6/4
star.name.inset <- c(-0.13, -0.05)

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

left.stars <- c("1435467", "5184732", "7680114", "8938364", "10516096",
    "5774694", "12069424")
bot.stars  <- c("10516096", "10963065", "11253226", "12009504",
    "5774694", "12069424", "12069449")

#rnorms <- sapply(1:1024, function(x) rnorm(num_knots, 0, 0.1))
#rnorms <- sapply(1:1024, function(x) runif(num_knots, -0.25, 0.25))
# generate random functions 
rnorms <- sapply(1:n_trials, function(x) 
    runif(num_knots, -0.3/2, 0.3/2) + runif(1, -0.25/2, 0.25/2))
# randomly set them to be equal to the previous value 
for (ii in 2:num_knots) {
    replacements <- rbinom(n_trials, 1, 0.5)
    indices <- replacements * 1:n_trials
    indices <- indices[which(indices > 0)]
    rnorms[ii, indices] <- rnorms[ii-1, indices] 
}

sunstar <- read.table(file.path('..', 'grid', 'ref_mods', 
        'inv', '5774694', 'LOGS_MS', 
        'profile1-freqs', 'profile1.data.FGONG.dat'),
    header=1)
sunstar.u <- with(sunstar, splinefun(x, u)(kern.interp.xs))

modmix <- read.table(file.path('models', 'BPB2000', 'modmix.dat'),
    header=1)
modmix.u <- with(modmix, splinefun(x, u)(kern.interp.xs))

calibrated <- read.table(file.path('..', 'grid', 'calibrate', 'no_diffusion_final',
        'LOGS_3MS', 'profile1-freqs', 'profile1.data.FGONG.dat'),
    header=1)
calibrated.u <- with(calibrated, splinefun(x, u)(kern.interp.xs))

calibrated_diff <- read.table(file.path('..', 'grid', 'calibrate', 'diffusion_final',
        'LOGS_3MS', 'profile1-freqs', 'profile1.data.FGONG.dat'),
    header=1)
calibrated.diff.u <- with(calibrated_diff, splinefun(x, u)(kern.interp.xs))

no_overshoot <- read.table(file.path('..', 'grid', 'overshoot', 'nos',
        'LOGS_3MS', 'profile1.data.FGONG.dat'),
    header=1)
no_overshoot.u <- with(no_overshoot, splinefun(x, u)(kern.interp.xs))

overshoot <- read.table(file.path('..', 'grid', 'overshoot', 'os',
        'LOGS_3MS', 'profile1.data.FGONG.dat'),
    header=1)
overshoot.u <- with(overshoot, splinefun(x, u)(kern.interp.xs))


plot_du.u <- function(du.u, inversion, sampler, star.name, errbars=F, 
        make_ylab=T, make_xlab=T, small=F, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font=utils.font, tcl=utils.tcl) {
    
    col.pal <- brewer.pal(11, "Spectral")[c(1:4,8:11)]
    col.pal <- col.pal[c(1:4, 7:8)]
    col.pal <- adjustcolor(col.pal, alpha.f=0.9)
    
    par(mar=mar+c(0.3 + ifelse(small, 0.3, 0),
                 -0.2 + ifelse(small, 0.2, 0), 
                 -0.5, -0.3), 
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0, 0.34), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lwd=par()$lwd, lty=2, col='black')
    
    intgrl <- function(ii) sintegral(kern.interp.xs, 
        du.u*inversion$avg_kerns[,ii])$value
    ress <- sapply(1:nrow(inversion$result), intgrl)
    resids <- ress - sapply(1:nrow(inversion$result), 
        function(ii) du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
    lines(kern.interp.xs, du.u, type='l', lwd=par()$lwd*2, lty=2, col="darkgray")#"#e76f51")
    
    if (errbars) {
        with(inversion$result[sampler,], 
            arrows(fwhm.left, ress[sampler], fwhm.right, ress[sampler], 
                col=col.pal[sampler],#adjustcolor("#264653", alpha.f=1),#'darkgray', 
                lwd=3, length=0.05, angle=90, code=3))
        with(inversion$result[sampler,], 
            arrows(fwhm.mid, ress[sampler]-err, fwhm.mid, ress[sampler]+err, 
                col=col.pal[sampler],#adjustcolor("#264653", alpha.f=1),#'darkgray', 
                lwd=3, length=0.05, angle=90, code=3))
    }
    #points(inversion$result$fwhm.mid[sampler], ress[sampler],
    #    col=1, pch=20, cex=1)
    points(inversion$result$fwhm.mid[sampler], ress[sampler],
        col=1, pch=20, cex=1)
    
    par(family='Helvetica')# LT Std Light')
    legend("bottomright", bty="n", cex=0.9*text.cex,
        legend=star.name, inset=c(-0.01, -0.05))
    par(family=font)
    
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4 + ifelse(small, 0.1, 0), 0), 
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, majorn=3, labels=make_xlab)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), 
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=make_ylab, tick=F, las=1, 
        cex.axis=text.cex)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55 + ifelse(small, 0.3, 0), 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_du.u.prof <- function(du.u.prof, inversion, sampler, star.name, 
        make_ylab=T, make_xlab=T, small=F, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font=utils.font, tcl=utils.tcl) {
    par(mar=mar+c(0.3 + ifelse(small, 0.3, 0),
                 -0.2 + ifelse(small, 0.2, 0), -0.5, -0.3), 
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         #xlim=c(0, 0.35), 
         xlim=c(0, 1),
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lty=2, col='black')
    
    for (ii in 1:nrow(du.u.prof)) {
        lines(kern.interp.xs, du.u.prof[ii,], lwd=0.1, 
            col=adjustcolor(1, 0.2))
    }
    offset <- 4
    cols <- c("#6699cc", "#fff275", "#ff8c42", "#ff3c38", "#a23e48")
    for (ii in 1:5+offset) {
        lines(kern.interp.xs, du.u.prof[ii,], lwd=1.5,
            col=cols[ii-offset])
            #col=c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")[ii-offset])
            #col=c("#ed6a5a", "#f4f1bb", "#9bc1bc", "#5ca4a9", "#e6ebe0")[ii-offset])
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
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4 + ifelse(small, 0.1, 0), 0), 
        cex.axis=text.cex, 
        family=font, majorn=3, labels=make_xlab)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=make_ylab, tick=F, las=1, cex.axis=text.cex)
    
    par(mgp=mgp+c(0.55 + ifelse(small, 0.3, 0), 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=bquote(delta*.(f1.exp)/.(f1.exp)))
}

plot_residuals_line <- function(residuals, means, errs, inversion, sampler, 
        star.name, make_ylab=T, make_xlab=T, small=F, 
        caption.inset=star.name.inset, make.legend=F, 
        ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font=utils.font, tcl=utils.tcl) {
    par(mar=mar+c(0.3 + ifelse(small, 0.3, 0),
                 -0.2 + ifelse(small, 0.2, 0), -0.5, -0.3), 
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0, 0.34), 
         ylim=c(0, 0.22),#ifelse(small, 0.22, 0.16)), 
         xlab="",
         ylab="")
    
    col.pal <- c("#282f44", "#f49e4c", "#fdd692") #"#0a122a")#c(red, 'darkgray')
    
    lines(inversion$result[sampler,]$fwhm.mid, errs[sampler], 
        lwd=3, lty=1, col=col.pal[2])
    with(inversion$result[sampler,], lines(fwhm.mid, err, 
        lwd=3, lty=2, col=col.pal[1]))
    
    points(inversion$result[sampler,]$fwhm.mid, errs[sampler], 
        pch=21, col=col.pal[2], bg=col.pal[3], lwd=2, cex=1.1)
    with(inversion$result[sampler,], points(fwhm.mid, err, 
        pch=22, col=col.pal[1], bg=col.pal[3], lwd=2, cex=1.2))
    
    
    par(family='Helvetica')# LT Std Light')
    legend("topleft", bty="n", cex=0.9*text.cex,
        legend=star.name, inset=caption.inset)
    if (F){
    if (make.legend) legend('top', ncol=2, bty='n', 
        lwd=3, lty=c(2, 1), pch=NA, col=col.pal, 
        legend=c('Uncertainties', 'Residuals'), text.col='white',
        cex=0.9*text.cex)
    if (make.legend) legend('top', ncol=2, bty='n', 
        lwd=2, lty=NA, pch=c(22, 21), col=col.pal, pt.bg=col.pal[3], 
        pt.cex=c(0.9, 1),
        legend=c('Uncertainties', 'Residuals'),
        cex=0.9*text.cex)
    }
    
    if (make.legend) legend('topleft', bty='n', inset=c(0.03, -0.01), # LINES
        lwd=3, lty=c(2, 1), pch=NA, col=col.pal[1], x.intersp=0.4, seg.len=2.7,
        legend=c('Uncertainties'), text.col='white', cex=0.9*text.cex)
    if (make.legend) legend('topleft', bty='n', inset=c(0.03, -0.01), # POINTS
        lwd=2, lty=NA, pch=c(22, 21), col=col.pal[1], x.intersp=0.4, 
        pt.bg=col.pal[3], pt.cex=c(1.1), legend=c('Uncertainties'),
        cex=0.9*text.cex, seg.len=2.7)
    if (make.legend) legend('topright', bty='n', inset=c(0.015, -0.01), # LINES
        lwd=3, lty=c(1), pch=NA, col=col.pal[2], x.intersp=0.6, seg.len=2,
        legend=c('Residuals'), text.col='white', cex=0.9*text.cex)
    if (make.legend) legend('topright', bty='n', inset=c(0.015, -0.01), # POINTS
        lwd=2, lty=NA, pch=c(21), col=col.pal[2], x.intersp=0.6, seg.len=2,
        pt.bg=col.pal[3], pt.cex=c(1.2), legend=c('Residuals'), 
        cex=0.9*text.cex)
    par(family=font)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4 + ifelse(small, 0.1, 0), 0), 
        cex.axis=text.cex, 
        family=font, majorn=3, labels=make_xlab, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=2, labels=make_ylab, lwd.ticks=par()$lwd)
    #axis(2, at=c(-0.2, 0.2), labels=make_ylab, tick=F, las=1, 
    #    cex.axis=text.cex)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55 + ifelse(small, 0.3, 0), 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=expression(delta*u/u))#"Residuals"))
}

plot_residuals <- function(residuals, means, errs, inversion, sampler, 
        star.name, 
        make_ylab=T, make_xlab=T, small=F, 
        caption.inset=star.name.inset, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font=utils.font, tcl=utils.tcl) {
    par(mar=mar+c(0.3 + ifelse(small, 0.3, 0),
                 -0.2 + ifelse(small, 0.2, 0), -0.5, -0.3), 
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0, 0.34), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    abline(h=0, lty=2, col='black', lwd=par()$lwd)
    for (ii in 1:nrow(residuals)) {
        with(inversion$result[sampler,],
            points(rnorm(length(rs), fwhm.mid, (fwhm.right-fwhm.left)/14),
                residuals[ii,sampler], 
                pch=2, cex=0.01, col='gray'))
    }
    
    with(inversion$result[sampler,], 
        arrows(fwhm.mid-0.003, err, fwhm.mid-0.003, -err, 
            code=3, angle=90, length=0.01, lwd=2, lty=1, col=1))
    with(inversion$result[sampler,], 
        arrows(fwhm.mid+0.003, means[sampler]+errs[sampler], 
               fwhm.mid+0.003, means[sampler]-errs[sampler], 
            code=3, angle=90, length=0.01, lwd=2, lty=1, col=red))
    
    #legend('topleft', lty=NA, lwd=1, pch=20, col=c(1, red), bty='n',
    #    cex=0.8*text.cex, x.intersp=0.01, inset=c(-0.02, -0.04),
    #    legend=c('Formal error bars', 'Monte Carlo'))
    #legend("bottomright", bty="n", cex=0.8*text.cex,
    #    legend=star.name, inset=c(-0.01, -0.05))
    legend("topleft", bty="n", cex=0.9*text.cex,
        legend=star.name, inset=caption.inset)
    par(family=font)
    
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4 + ifelse(small, 0.1, 0), 0), 
        cex.axis=text.cex, 
        family=font, majorn=3, labels=make_xlab, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F, lwd.ticks=par()$lwd)
    axis(2, at=c(-0.2, 0.2), labels=make_ylab, tick=F, las=1, 
        cex.axis=text.cex)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55 + ifelse(small, 0.3, 0), 0, 0))
    if (make_xlab) title(xlab=expression("Radius"~r/R))
    par(mgp=mgp+c(0.8, 0, 0))
    f1.exp <- inversion$k.pair$f1.exp
    if (make_ylab) title(ylab=expression("Residuals"))
}

plot_residuals2 <- function(du.u.avg, du.u.true, inversion, sampler, 
        star.name, make_legend=F, 
        make_ylab=T, make_xlab=T, small=F, 
        caption.inset=star.name.inset, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, 
        font=utils.font, tcl=utils.tcl) {
    #par(mar=mar+c(1.2,-0.2,-0.4,-0.3), mgp=mgp+c(0, 0.4, 0)) #, pty="s")
    par(mar=mar+c(0.3 + ifelse(small, 0.3, 0),
                 -0.2 + ifelse(small, 0.2, 0), -0.5, -0.3), 
        mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(-0.25, 0.25), 
         ylim=c(-0.25, 0.25),
         xlab="",
         ylab="")
    #abline(h=0, lty=2, col='black')
    segments(-1, -1, 1, 1, lty=2, lwd=par()$lwd)
    for (ii in 1:nrow(inversion$result)) {
        if (!sampler[ii]) next 
        points(du.u.true[,ii], du.u.avg[,ii],
            pch=2, cex=0.2, col=adjustcolor(col.pal2[ii], alpha.f=0.5))
    }
    
    par(family='Helvetica')# LT Std Light')
    if (F) {
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
    if (make_legend)
        legend("topleft", bty="n", pch=20, cex=0.8*text.cex, 
            inset=c(-0.01, -0.025), legend=as.expression(legend.txt), 
            col=col.pal2[sampler], x.intersp=0.66)
    }
    
    legend("topleft", bty="n", cex=0.9*text.cex,
        legend=star.name, inset=caption.inset)#c(-0.01, -0.03))
    par(family=font)
    
    magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4 + ifelse(small, 0.2, 0), 0), 
        cex.axis=text.cex, 
        family=font, majorn=3, labels=F)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=F)
    axis(2, at=c(-0.2, 0.2), labels=make_ylab, tick=F, las=1, cex.axis=text.cex)
    axis(1, at=c(-0.2, 0.2), labels=make_xlab, tick=F, las=1, cex.axis=text.cex)
    
    #par(mgp=mgp+c(1.7, 0, 0))
    par(mgp=mgp+c(0.55 + ifelse(small, 0.3, 0), 0, 0))
    #title(xlab=expression(integral(K^"avg"~frac(du,u)~dr)))
    #title(xlab=expression(frac(du,u)(r[0])))
    if (make_xlab) title(xlab=expression("true"~du/u(r[0])))
    par(mgp=mgp+c(0.8, 0, 0))
    if (make_ylab) title(ylab=expression("<"*du/u*">"[r[0]]))
}

S <- function(x, alpha, ks, degree, num_knots) {
    sapply(x, function(x) sum(sapply( degree:(num_knots+degree-1),
        function(ii) alpha[ii-degree+1] * B(x, ii, degree, ks) )))
    #(degree+1):(num_knots+degree-1),
    #1:(num_knots-(degree+1)), 
        #function(ii) alpha[ii-degree] * B(x, ii, degree, ks) )))
}


check_kerns <- function(star.name, inversion, ref.mod, sampler, sun=F) {
    ## calculate acoustic depth, Aerts et al eq 3.228 
    int_0.r <- function(x, y) as.numeric(cumtrapz(x, y)) # int_0^r 
    int_r.R <- function(x, y) trapz(x, y) - int_0.r(x, y) # int_r^R 
    tau <- int_r.R(kern.interp.xs*R_sun*ref.mod$R, 
        1/ref.mod$cs.spl(kern.interp.xs))
    int.ks <- seq(min(tau), max(tau), length.out=20)
    ks <- c(-1, -1, -1, splinefun(tau, kern.interp.xs)(rev(int.ks)), 2,2,2)
    
    du.u <- (calibrated.u - calibrated.diff.u)/calibrated.u
    make_plots(plot_du.u, paste0(star.name, '-true_diff-solar'), #'true_diff', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        errbars=T, small=T, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        star.name=star.name, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        du.u=du.u, inversion=inversion, sampler=sampler)
    if (sun) make_plots(plot_du.u, paste0(star.name, '-true_diff-solar'), 
        filepath=file.path('plots', 'sun'),
        slides=F, wide=F, tall=F, make_png=F, 
        errbars=T, 
        cex.paper=cex.sun, 
        star.name="",
        du.u=du.u, inversion=inversion, sampler=sampler)
    
    du.u <- (calibrated.u - modmix.u)/calibrated.u
    make_plots(plot_du.u, paste0(star.name, '-true_diff-modmix'), #'true_diff', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        star.name=star.name, 
        errbars=T, small=T, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        du.u=du.u, inversion=inversion, sampler=sampler)
    if (sun) 
    make_plots(plot_du.u, paste0(star.name, '-true_diff-modmix'), #'true_diff', 
        filepath=file.path('plots', 'sun'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.sun, 
        star.name="",
        errbars=T, 
        du.u=du.u, inversion=inversion, sampler=sampler)
    
    if (F) {
    du.u <- S(kern.interp.xs, rnorms[,1], ks, degree, num_knots)
    make_plots(plot_du.u, paste0(star.name, '-true_diff-Bspline'), #'true_diff', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        errbars=T, small=T, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        star.name=star.name, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        du.u=du.u, inversion=inversion, sampler=sampler)
    }
    
    du.u.file <- file.path('save-poly', paste0('du.u.prof-', star.name))
    du.u.avg.file <- file.path('save-poly', paste0('du.u.avg-', star.name))
    du.u.true.file <- file.path('save-poly', paste0('du.u.true-', star.name))
    if (file.exists(du.u.file)) {
        load(du.u.file)
        load(du.u.avg.file)
        load(du.u.true.file)
    } else {
        du.u.prof <- do.call(rbind, parallelMap(function(ii) {
            S(kern.interp.xs, rnorms[,ii], ks, degree, num_knots)
        }, ii=1:ncol(rnorms)))#128))#
        save(du.u.prof, file=du.u.file)
        
        du.u.avg <- do.call(rbind, Map(function(ii) {
            #du.u <- S(kern.interp.xs, 
            #    rnorms[,ii], ks, degree, num_knots)
            du.u <- du.u.prof[ii,]
            intgrl <- function(ii) sintegral(kern.interp.xs, 
                du.u*inversion$avg_kerns[,ii])$value
            sapply(1:nrow(inversion$result), intgrl)
        }, ii=1:nrow(du.u.prof)))
        save(du.u.avg, file=du.u.avg.file)
        
        du.u.true <- do.call(rbind, Map(function(ii) {
            du.u <- du.u.prof[ii,]
            sapply(1:nrow(inversion$result), function(ii) 
                du.u[which(kern.interp.xs == inversion$result$fwhm.mid[ii])])
        }, ii=1:nrow(du.u.prof)))
        save(du.u.true, file=du.u.true.file)
    }
    
    if (F) {
    make_plots(plot_du.u.prof, 'duprof', 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_pdf=F,
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        du.u.prof=du.u.prof, small=T, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        star.name=star.name, sampler=sampler, inversion=inversion)
    }
    
    residuals <- du.u.avg - du.u.true 
    errs <- apply(residuals, 2, std)
    means <- apply(residuals, 2, mean)
    
    print(star.name)
    print(sampler)
    print(inversion$result$err*1.5 > errs)
    
    if (F) {
    make_plots(plot_residuals, paste0(star.name, '-true_diff-residuals'), 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        sampler=sampler, small=T, 
        star.name=star.name, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        residuals=residuals, means=means, errs=errs, inversion=inversion)
    if (sun) make_plots(plot_residuals, 
        paste0(star.name, '-true_diff-residuals'), 
        filepath=file.path('plots', 'sun'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.sun, 
        sampler=sampler,
        star.name="", 
        residuals=residuals, means=means, errs=errs, inversion=inversion)
    }
    
    make_plots(plot_residuals_line, 
        paste0(star.name, '-true_diff-residuals-line'), 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        sampler=sampler, small=T, 
        star.name=star.name, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        residuals=residuals, means=means, errs=errs, inversion=inversion)
    if (sun) make_plots(plot_residuals_line, 
        paste0(star.name, '-true_diff-residuals-line'), 
        filepath=file.path('plots', 'sun'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.sun, 
        sampler=sampler,
        paper_pdf_height=5.15, 
        star.name="", make.legend=T, 
        residuals=residuals, means=means, errs=errs, inversion=inversion)
    
    if (F) make_plots(plot_residuals2, paste0(star.name, '-true_diff-residuals2'), 
        filepath=file.path('plots', 'tests'),
        slides=F, wide=F, tall=F, make_png=F, 
        cex.paper=cex.paper, #paper_pdf_height=6.75206*2/3, 
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width,
        sampler=sampler, small=T, 
        star.name=star.name, 
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        du.u.avg=du.u.avg, du.u.true=du.u.true, inversion=inversion)
}



DF <- NULL

for (star.name_i in 1:length(star.names)) {
    
    star.name <- star.names[star.name_i]
    print(star.name)
    
    sun <- star.name == '5774694'
    cyg <- star.name == '12069424' || star.name == '12069449'
    
    rm.first.kernel <- star.name %in% c('12009504', '8179536', '1435467', 
        '2837475', '9812850', '3427720', '8938364', '9139163', '10454113',
        '10516096', '11253226', '8006161')
    rm.second.kernel <- star.name %in% c('12009504', '2837475', '8179536',
        '9139163', '9812850', '1435467', '3427720', '10454113', '11253226',
        '7680114', '8006161')
    rm.third.kernel <- star.name %in% c()
    rm.fourth.kernel <- star.name %in% c('9139163', '8938364', '10454113', '5774694')
    rm.fifth.kernel <- star.name %in% c('1435467', '9139163', '10454113')
    rm.last.kernel <- star.name %in% c('3427720', '5774694', '6106415', 
        '7680114', '8006161', '8179536', '9098294', '10963065', '12009504',
        '1435467', '6603624', '9139163', '8938364', '9812850', '10454113')
    
    for (filename in save.files[grep(star.name, save.files)]) {
        load(file.path('save', filename)) 
        # loads ref.mod, inv.lists, avg.kerns.lists, cross.lists, MRs, inv.params
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
            kern.interp.xs=kern.interp.xs, fctnl='median') 
        
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
        
        sampler <- penalties < 4
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
            check_kerns(star.name, inversion, ref.mod, sampler, sun) 
            
            #if (F) 
            make_plots_inversion_all(model, inversion, 
                filepath=file.path('plots', 'phys'),
                k.str=paste0('-', star.name, '-', phy), 
                cross.inset="bottomright",
                caption=star.name, 
                caption.inset=star.name.inset,
                kern.capt=T, 
                small=T, 
                #col.pal=col.pal2[sampler],#col.pal[phy_i], 
                #caption2=as.expression(bquote('M' == .(mass))),
                sampler=sampler,
                suppress_cross=phy_i>1,
                suppress_avg=phy_i>1,
                #core.bound=core.bound, 
                kern.interp.xs=kern.interp.xs,
                xlim=c(0, 0.34),
                make_xlabs=rep(star.name %in% bot.stars, 3),
                make_ylabs=rep(star.name %in% left.stars, 3),
                inversion_ylim=c(-0.4, 0.4),
                cross_kern_ylim=c(-0.7, 0.3),
                avg_kern_ylim=c(-5, 20),
                slides=F, wide=F, tall=F, make_png=F, 
                cex.paper=cex.paper,
                paper_pdf_height=paper_pdf_height,
                paper_pdf_width=paper_pdf_width)#, make.inset=T)
            
            if (sun) {
                make_plots_inversion_all(model, inversion, 
                    filepath=file.path('plots', 'sun'),
                    k.str=paste0('-', star.name, '-', phy), 
                    cross.inset="bottomright",
                    #caption=star.name, 
                    #caption.inset=star.name.inset,
                    #kern.capt=T, 
                    sampler=sampler,
                    suppress_cross=phy_i>1,
                    suppress_avg=phy_i>1,
                    #core.bound=core.bound, 
                    kern.interp.xs=kern.interp.xs,
                    xlim=c(0, 0.34),
                    make_xlabs=c(T,T,T), #rep(star.name %in% bot.stars, 3),
                    make_ylabs=c(T,T,T), #rep(star.name %in% left.stars, 3),
                    inversion_ylim=c(-0.4, 0.4),
                    cross_kern_ylim=c(-0.7, 0.3),
                    avg_kern_ylim=c(-4, 20),
                    slides=F, wide=T, tall=T, make_png=T, 
                    cex.paper=cex.sun)#,
            }
            
            DF.obs <- read.table(file.path('..', 'regression', 'data',
                'inversions', paste0(star.name, '-obs.dat')), header=1)
            DF.rf  <- read.table(file.path('..', 'regression', 'learn-inv',
                'covs-inv', 'inversions', paste0(star.name, '.dat')), header=1)
            
            interf <- F
            if ('radius' %in% names(DF.rf)) {
                R <- mean(DF.rf$radius)
                e_R <- sd(DF.rf$radius)
            } else {
                R <- DF.obs[DF.obs$name == 'radius',]$value
                e_R <- DF.obs[DF.obs$name == 'radius',]$uncertainty
                interf <- T
            }
            
            params <- data.frame(
                KIC=star.name,
                age=mean(DF.rf$age),
                    e_age=sd(DF.rf$age),
                M=mean(DF.rf$M),
                    e_M=sd(DF.rf$M),
                R=R,
                    e_R=e_R,
                FeH=DF.obs[DF.obs$name == 'Fe/H',]$value,
                    e_FeH=DF.obs[DF.obs$name == 'Fe/H',]$uncertainty,
                alpha=mean(DF.rf$alpha),
                    e_alpha=sd(DF.rf$alpha)
            )
            
            inv.res <- data.frame(row.names=1)
            for (res.i in c(1:length(sampler))[sampler]) {
                inv.res0 <- with(inversion$result[res.i,], data.frame(
                    r0=df_dr,
                    e_r0=err,
                    r02=f/10**15,
                    e_r02=f.err/10**15,
                    rmax=fwhm.mid,
                    fwhm=(fwhm.right-fwhm.left)/2))
                r0 <- paste0('tr', gsub('\\.', '', inversion$result$rs[res.i]))
                r02 <- paste0('fr', gsub('\\.', '', inversion$result$rs[res.i]))
                rmax <- paste0('rmax', gsub('\\.', '', inversion$result$rs[res.i]))
                names(inv.res0) <- c(r0, paste0('e_', r0),
                    r02, paste0('e_', r02),
                    rmax, paste0('e_', rmax))
                inv.res <- cbind(inv.res, inv.res0)
            }
            
            DF <- plyr:::rbind.fill(DF, cbind(params, inv.res))
        }
    }
    
    #caption2 <- sapply(1:length(phys), function(ii) {
    #    paste0("M=", phy_mass[[ii]], ", R=", phy_radius[[ii]], 
    #        if (ii==2) " (ov)" else if (ii==3) " (D)" else if (ii==4) " (D, ov)" else "")
    #})
    #caption2 <- c('No diffusion, no overshoot',
    #    'No diffusion, overshoot', 
    #    'Diffusion, no overshoot',
    #    'Diffusion, overshoot')
    caption2 <- c("Vanilla", "Overshoot", "Diffusion", "Both")
    
    #if (F) 
    make_plots(plot_inversion_multi, star.name, 
        filepath=file.path('plots', 'phys-comp'),
        slides=F, wide=F, tall=F, make_png=F, small=T, 
        inversion.list=phy_res,
        kern.interp.xs=kern.interp.xs,
        #ylim = c(-0.4, 0.4),
        ylim = c(-0.25, 0.25),
        xlim=c(0, 0.34),
        caption=star.name, 
        caption.inset=star.name.inset,
        #caption2=if (star.name == '3427720' || star.name == '8179536') 
        #    caption2 else NULL,
        make_ylab=star.name %in% left.stars,
        make_xlab=star.name %in% bot.stars,
        caption2=NULL, 
        col.pal=col.pal, 
        cex.paper=cex.paper,
        paper_pdf_height=paper_pdf_height,
        paper_pdf_width=paper_pdf_width)#*6/8)
        #caption2=as.expression(bquote('M' == .(mass))))
    if (sun || cyg) make_plots(plot_inversion_multi, star.name, 
        filepath=file.path('plots', 'sun'),
        slides=F, wide=F, tall=F, make_png=F, 
        paper_pdf_height=5.15,
        inversion.list=phy_res,
        kern.interp.xs=kern.interp.xs,
        sampler=sampler, 
        caption2=if(sun) caption2 else NULL,
        ylim = c(-0.25, 0.25),
        xlim=c(0, 0.34),
        col.pal=col.pal, 
        cex.paper=cex.sun)
}

if (F) {

DF <- DF[-which(DF$KIC == '5774694'),]
DF <- DF[order(as.numeric(paste0(DF$KIC))),]
#DF <- DF[order(as.numeric(levels(DF$KIC))),]
#DF <- DF[order(DF$KIC),]
write.table(DF, file='inv-res.dat', quote=F, row.names=F, sep='\t')
write.table(DF[,-grep('tr', names(DF))], 
    file='inv-res2.dat', quote=F, row.names=F, sep='\t')
write.table(DF[,-grep('tr', names(DF))][,-c(2:11)], 
    file='inv-res3.dat', quote=F, row.names=F, sep='\t')
write.table(DF[,grep('KIC|rmax|e_rmax', names(DF))],
    file='inv-res-rmax.dat', quote=F, row.names=F, sep='\t')

# (.+?)\t(.+?)\t        \1\t&\t\2\t$\\\pm$\t
# \$\\pm\$\te_(.+?)\t   \t\t
# (\$\\pm\$\t\d\.(0*)\d\d)(\d*?)\t     \1\t
# NA   -
# \n   \\\\\n
# -\t$\pm$\t-  -
# -	$\pm$	-    -


## REDO INVERSION 
if (F) {
star.name <- '12069449'#'5184732'
for (filename in save.files[grep(star.name, save.files)]) {
    load(file.path('save', filename)) 
    # loads ref.mod, inv.lists, avg.kerns.lists, cross.lists, MRs, inv.params
}
#rs = seq(0.05, 0.3, 0.05) 
#if (star.name == '5184732' || star.name == '12069449') {
rs = seq(0.05, 0.3, 0.01) 

model <- ref.mod
model.names <- names.list[[star.name]]
target.name <- star.name
freqs <- get_freqs(paste0('KIC_', target.name))
model.list <- parallelMap(function(model.name) { print(model.name); 
        get_model(freqs=freqs, model.name=model.name, 
            target.name=target.name, k.pair=k.pair, square.Ks=T) }, 
    model.name=model.names)
names(model.list) <- model.names

set.seed(0)

inv.lists <- list()
avg.kerns.lists <- list()
cross.lists <- list()
for (trial_i in 1:length(inv.params$error.sups)) {
#inv.lists <- Map(function(trial_i) {
    print(trial_i)
    cross.term <- inv.params$cross.terms[trial_i]
    width <- inv.params$widths[trial_i]
    error.sup <- inv.params$error.sups[trial_i]
    star.M <- MRs$Ms[trial_i]
    star.R <- MRs$Rs[trial_i]
    
    nu.star <- rnorm(length(freqs$nu), freqs$nu, freqs$dnu)
    
    #radial_velocity <- -28.02
    #doppler_beta <- radial_velocity/speed_of_light
    #doppler_shift <- sqrt((1+doppler_beta)/(1-doppler_beta))
    #nu.star <- nu.star * doppler_shift
    
    repeat {
        trial.M1 <- rnorm(1)#, initial.M, sigma.M)
        trial.R1 <- rnorm(1)#, initial.R, sigma.R)
        #if (trial.M > initial.M - 10*sigma.M &&
        #    trial.M < initial.M + 10*sigma.M &&
        #    trial.R > initial.R - 10*sigma.R && 
        #    trial.R < initial.R + 10*sigma.R) break
        break 
    }
    for (perturbed.model.name in model.names) {
        modl <- model.list[[perturbed.model.name]]
        modl$nus$nu.y <- nu.star
        modl$F_surf <- get_F_surf(modl$nus, num.knots=0, use.BG=T, 
            nu_ac=modl$nu_ac)
        model.list[[perturbed.model.name]] <- modl
    }
    
    
    inv.list <- parallelMap(function(model_i) {
        model <- model.list[[model_i]]
        invert.OLA(model=model, rs=rs, 
            cross.term=cross.term, error.sup=error.sup, width=width, 
            targ.kern.type=targ.kern.type, 
            get_cross_kerns=T, get_avg_kerns=T, 
            kern.interp.xs=kern.interp.xs,
            subtract.mean=F, 
            dM=model$M-star.M, dR=model$R-star.R, 
            perturb=F, num_realizations=1, F_surf=model$F_surf)
    }, model_i=1:length(model.list))
    names(inv.list) <- model.names
    
    #inv.list
    #Map(function(inv) inv$result, inv=inv.list)
    
            
    inv.lists[[trial_i]] <- Map(function(inv) inv$result, inv=inv.list)
    avg.kerns.lists[[trial_i]] <- Map(function(inv) inv$avg_kern, inv=inv.list)
    cross.lists[[trial_i]] <- Map(function(inv) inv$cross_kern, inv=inv.list)
    
}#, trial_i=1:length(inv.params$error.sups))#9)#
#}
}


}
