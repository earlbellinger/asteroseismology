#### Animate the construction of averaging kernels 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

source('../scripts/utils.R') 
source('../scripts/seismology.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')

col.pal <- c(blue, "darkgray", "#F29559", "#DB4D48")

### CONSTANTS 
k.pair  = u_Y
rs      = seq(0.05, 0.3, 0.05) 
sampler = c(T)
models  = get_model_list() 
kern.interp.xs = seq(0, 1, 0.001)
targ.kern.type <- 'mod_Gauss'

targ.mode <- '-p_CygAwball-m_CygA-e_CygA-r_diffusion-mod_Gauss'
load(file.path('save-paper', paste0('ref.mod',         targ.mode)))
load(file.path('save-paper', paste0('inv.lists',       targ.mode)))
load(file.path('save-paper', paste0('avg.kerns.lists', targ.mode)))
load(file.path('save-paper', paste0('cross.lists',     targ.mode)))
load(file.path('save-paper', paste0('MRs',             targ.mode)))
load(file.path('save-paper', paste0('inv.params',      targ.mode)))

ii <- 1
r0i <- 2

cross.term <- inv.params$cross.terms[ii]
error.sup <- inv.params$error.sups[ii]
width <- inv.params$widths[ii]
star.M <- MRs$Ms[ii]
star.R <- MRs$Rs[ii]

inversion <- invert.OLA(model=ref.mod, 
    rs=rs, 
    cross.term=cross.term, 
    error.sup=error.sup, 
    width=width,
    use.BG=T, 
    BG_pows=c(-2, 2), 
    num.knots=0, 
    targ.kern.type=targ.kern.type, 
    num_realizations=1, 
    r_f=0.2, 
    perturb=F, 
    dM=ref.mod$M-star.M, 
    dR=ref.mod$R-star.R, 
    subtract.mean=F, 
    kern.interp.xs=kern.interp.xs,
    get_avg_kerns=T, 
    get_cross_kerns=T, 
    verbose=T)

inv.coefs <- inversion$inv.coefs
colnames(inv.coefs) <- paste0('c', 1:ncol(inv.coefs))
nus <- cbind(ref.mod$nus, inv.coefs)


radii <- ref.mod$r * ref.mod$R * solar_radius
r_ts <- sapply(ref.mod$modes, function(mode) {
    ell <- as.numeric(strsplit(strsplit(mode, '_')[[1]][1], '\\.')[[1]][2])
    nn <- as.numeric(strsplit(strsplit(mode, '_')[[1]][2], '\\.')[[1]][2])
    nu <- ref.mod$nus[with(ref.mod$nus, l==ell&n==nn),]$nu.x
    #if (ell == 0) return(0)
    ref.mod$r[which.min((
            ref.mod$cs.spl(radii / (ref.mod$R * solar_radius))**2 / (radii**2) 
            - (10**-6*nu*(2*pi))**2/(ell+1/2)#ell*(ell+1))
        )**2)]
})

nus <- cbind(nus, data.frame(r_t = r_ts))

nu_max <- 2201

seis.DF <- seismology(with(nus, data.frame(l=l, n=n, nu=nu.y, dnu=dnu)), 
    nu_max=nu_max)

kerns <- sapply(1:nrow(nus), 
    function(ii) max(abs(
        inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]]
)))
kern.max <- max(kerns)
#mode.order <- order(inv.coefs[, r0i], decreasing=T)
mode.order <- order(sapply(nus$nu.y, function(x) (x - nu_max)**2))
#mode.order <- order(nus$nu.y, decreasing=F)
#mode.order <- 1:nrow(inv.coefs[,r0i])


if (F) {
avg_kern <- colSums(t(inv.coefs[, r0i] * ref.mod$k1[,rownames(inv.coefs)]))
remaining <- 1:nrow(nus)
mode.order <- c()
for (ii in 1:nrow(nus)) {
    best_val <- Inf
    best_jj <- -1
    #best_newk <- avg_kern
    for (jj in remaining) {
        rem.idx <- which(remaining == jj)
        new_k <- colSums(t(inv.coefs[remaining[-rem.idx], r0i] * 
            ref.mod$k1[,rownames(inv.coefs)[remaining[-rem.idx]]]))
        
        val <- sum((avg_kern - new_k)**2)
        if (val < best_val) {
            best_val <- val
            #best_newk <- new_k
            best_jj <- jj
        }
    }
    #avg_kern <- best_newk
    remaining <- remaining[-which(remaining == best_jj)]
    mode.order <- c(mode.order, best_jj)
}
}

plot_echelle <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    ylim <- c(1500, 3500)#range(nus$nu.y)#c(, 3200)
    xlim <- c(0, seis.DF$Dnu0)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim) 
    
    xs <- ref.mod$nus$nu.y%%seis.DF$Dnu0
    xs[xs < 0] <- xs[xs < 0] + seis.DF$Dnu0
    
    with(ref.mod$nus, 
        points(xs, nu.y,
            pch=21, cex=1.5, lwd=1.5, 
            bg=col.pal[l+1], 
            col=ifelse(1:length(nu.y) %in% mode.order[1:modes], 1, 'white')))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    text(x=1+with(ref.mod$nus, max(nu.y[l==2])%%seis.DF$Dnu0), 
        y=3400, cex=text.cex, labels='2', col=col.pal[2+1], pos=1)
    text(x=1+with(ref.mod$nus, max(nu.y[l==0])%%seis.DF$Dnu0), 
        y=3400, cex=text.cex, labels='0', col=col.pal[0+1], pos=1)
    text(x=1+with(ref.mod$nus, max(nu.y[l==3])%%seis.DF$Dnu0), 
        y=3400, cex=text.cex, labels='3', col=col.pal[3+1], pos=1)
    text(x=1+with(ref.mod$nus, max(nu.y[l==1])%%seis.DF$Dnu0), 
        y=3400, cex=text.cex, labels='1', col=col.pal[1+1], pos=1)
    par(family=font)
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression((nu~mod~Delta*nu)/mu*Hz))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(nu/mu*Hz))
}

if (F) make_plots(plot_echelle, 'echelle', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 


plot_inv.coefs_ltp <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', xlab="", ylab="", 
        xlim=c(0, 0.1),
        ylim=c(-1.2, 1.2)) #c(-350, 350))
    
    #abline(v=0.05, lwd=1.5, lty=2)
    
    with(nus, 
        segments(r_t, 
                 0, 
                 r_t, 
                 inv.coefs[, r0i]/max(abs(inv.coefs[, r0i])), 
            col=col.pal[l+1], 
            lwd=3))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2223, y=-0.68, cex=0.9*text.cex, labels='0', col=col.pal[1])
    #text(x=0.04, y=0.5,   cex=0.9*text.cex, labels='1', col=col.pal[2])
    #text(x=0.085, y=0.5, cex=0.9*text.cex, labels='2', col=col.pal[3])
    #text(x=0.142, y=-0.2, cex=0.9*text.cex, labels='3', col=col.pal[4])
    par(family=font)
    
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(Lower~turning~point~r[t]/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Inversion~coefficient))
}

if (F) make_plots(plot_inv.coefs_ltp, 'inv_coefs', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 



plot_uykern <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    ylim <- c(-0.1, 0.1)
    xlim <- c(0, 1)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim) 
    
    #for (ii in 1:nrow(nus)) {
    #    lines(ref.mod$k1$x, 
    #          inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]],
    #          col=adjustcolor(1, alpha.f=0.25))
    #}
    if (F) {
    for (ii in mode.order[1:modes]) {
        lines(ref.mod$k1$x, 
            inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]] / kern.max,
            col=if (ii == mode.order[modes]) 
                col.pal[nus$l[ii]+1] else adjustcolor(1, alpha.f=0.25),
            lwd=if (ii == mode.order[modes]) 
                1.5 else 0.5)
    }
    }
    
    
    for (ii in 1:nrow(nus)) {
        lines(ref.mod.tmp$k1$x, 
            inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]] / kern.max,
            col=if (ii == nrow(ref.mod.tmp$nus)) 
                col.pal[nus$l[ii]+1] else adjustcolor(1, alpha.f=0.25),
            lwd=if (ii == nrow(ref.mod.tmp$nus)) 
                1.5 else 0.5)
    }
    
    
    abline(h=0, lwd=1.5, lty=2)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2003, y=0.725, cex=0.9*text.cex, labels='2', col=col.pal[3])
    par(family=font)
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(r/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression((u*','*Y)~Kernel))
}

if (F) make_plots(plot_uykern, 'uykern', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 



plot_yukern <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    ylim <- c(-0.1, 0.1)
    xlim <- c(0, 1)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim) 
    
    #for (ii in 1:nrow(nus)) {
    #    lines(ref.mod$k1$x, 
    #          inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]],
    #          col=adjustcolor(1, alpha.f=0.25))
    #}
    
    if (F) {
    for (ii in mode.order[1:modes]) {
        #lines(ref.mod$k2$x, 
        #      inv.coefs[ii, r0i] * ref.mod$k2[,rownames(inv.coefs)[ii]],
        #      col=col.pal[nus$l[ii]+1],
        #      lwd=1.5)
        lines(ref.mod$k2$x, 
            inv.coefs[ii, r0i] * ref.mod$k2[,rownames(inv.coefs)[ii]] / kern.max,
            col=if (ii == mode.order[modes]) 
                col.pal[nus$l[ii]+1] else adjustcolor(1, alpha.f=0.25),
            lwd=if (ii == mode.order[modes]) 
                1.5 else 0.5)
    }
    }
    
    
    for (ii in 1:nrow(nus)) {
        #lines(ref.mod$k2$x, 
        #      inv.coefs[ii, r0i] * ref.mod$k2[,rownames(inv.coefs)[ii]],
        #      col=col.pal[nus$l[ii]+1],
        #      lwd=1.5)
        lines(ref.mod.tmp$k2$x, 
            inv.coefs[ii, r0i] * ref.mod$k2[,rownames(inv.coefs)[ii]] / kern.max,
            col=if (ii == nrow(ref.mod.tmp$nus)) 
                col.pal[nus$l[ii]+1] else adjustcolor(1, alpha.f=0.25),
            lwd=if (ii == nrow(ref.mod.tmp$nus)) 
                1.5 else 0.5)
    }
    
    
    abline(h=0, lwd=1.5, lty=2)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2003, y=0.725, cex=0.9*text.cex, labels='2', col=col.pal[3])
    par(family=font)
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(r/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression((Y*','*u)~Kernel))
}

if (F) make_plots(plot_yukern, 'yukern', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 



plot_avg_kern <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    ylim <- c(-1, 1)
    xlim <- c(0, 1)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim) 
    
    #with(nus, 
    #    points(nu.y%%seis.DF$Dnu0 - 20, nu.y,
    #        pch=21, cex=1.5, lwd=1.5, 
    #        bg=c(1, blue, orange, red)[nus$l+1], 
    #        col='white'))
    
    #abline(h=0, lwd=1.5, lty=2)
    
    if (F) {
    kerns <- sapply(mode.order[1:modes], 
        function(ii) inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]])
    avg_kern <- colSums(t(kerns))
    lines(ref.mod$k1$x, avg_kern / (1.05*max(abs(avg_kern))), 
        lwd=1.5, col=blue)
    }
    
    lines(kern.interp.xs, #ref.mod.tmp$k1$x, 
        inversion$avg_kerns[, r0i] / (1.05*max(abs(inversion$avg_kerns[, r0i]))), 
        lwd=1.5, col=blue)
    
    #lines(kern.interp.xs, 
    #    inversion$avg_kerns[, r0i],
    #    lwd=1.5, col=blue)
    
    lines(ref.mod$k1$x,
        inversion$targ_kerns[, r0i] / (1.05*max(abs(inversion$targ_kerns[, r0i]))),
        lwd=1.5, lty=2, col='darkgray')
    
    abline(h=0, lwd=1.5, lty=2)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2003, y=0.725, cex=0.9*text.cex, labels='2', col=col.pal[3])
    par(family=font)
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(r/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Averaging~Kernel))
}

if (F) make_plots(plot_avg_kern, 'avg_kern', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 


plot_cross_kern <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,0.1,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    ylim <- c(-1, 1)
    xlim <- c(0, 1)
    
    plot(NA, axes=F, 
        xaxs='i',  yaxs='i', 
        xlab="",   ylab="", 
        xlim=xlim, ylim=ylim) 
    
    if (F) {
    kerns <- sapply(mode.order[1:modes], 
        function(ii) inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]])
    avg_kern <- colSums(t(kerns))
    
    kerns <- sapply(mode.order[1:modes], 
        function(ii) inv.coefs[ii, r0i] * ref.mod$k2[,rownames(inv.coefs)[ii]])
    lines(ref.mod$k1$x, colSums(t(kerns)) / (1.05*max(abs(avg_kern))), 
        lwd=1.5, col=orange)
    }
    
    lines(kern.interp.xs, #ref.mod.tmp$k1$x, 
        inversion$cross_kerns[, r0i] / (1.05*max(abs(inversion$avg_kerns[, r0i]))), 
        lwd=1.5, col=blue)
    
    abline(h=0, lwd=1.5, lty=2)
    
    #lines(kern.interp.xs, 
    #    inversion$avg_kerns[, r0i], 
    #    lwd=1.5, col=blue)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2003, y=0.725, cex=0.9*text.cex, labels='2', col=col.pal[3])
    par(family=font)
    #legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
    #    lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    #abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(r/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Cross-Term~Kernel))
}

if (F) make_plots(plot_cross_kern, 'cross_kern', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 








plot_all <- function(modes=5, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mfrow=c(3,2), oma=c(1, 1, 1, 1))
    plot_echelle(modes, ...)
    plot_inv.coefs_ltp(modes, ...)
    plot_uykern(modes, ...)
    plot_yukern(modes, ...)
    plot_avg_kern(modes, ...)
    plot_cross_kern(modes, ...)
    
}

if (F) make_plots(plot_all, 'all', 
    filepath=file.path('plots', 'inv_anim'), 
    slides=F, make_pdf=F, thin=F, short=F,
    make_png=T, 
    paper_pdf_width=513.11743 / 72.27,
    paper_pdf_height=630 / 72.27,
    cex.paper=1.16) 


for (modes in 3:nrow(ref.mod$nus)) {

    ref.mod.tmp <- ref.mod 
    
    ref.mod.tmp$modes <- ref.mod.tmp$modes[mode.order[1:modes]]
    ref.mod.tmp$k1 <- ref.mod.tmp$k1[,c('x', ref.mod$modes[mode.order[1:modes]])]
    ref.mod.tmp$k2 <- ref.mod.tmp$k2[,c('x', ref.mod$modes[mode.order[1:modes]])]
    ref.mod.tmp$nus <- ref.mod.tmp$nus[mode.order[1:modes],]
    
    inversion <- invert.OLA(model=ref.mod.tmp, 
        rs=rs, 
        cross.term=cross.term, 
        error.sup=error.sup, 
        width=width,
        use.BG=T, 
        BG_pows=c(-2, 2), 
        num.knots=0, 
        targ.kern.type=targ.kern.type, 
        num_realizations=1, 
        r_f=0.2, 
        perturb=F, 
        dM=ref.mod$M-star.M, 
        dR=ref.mod$R-star.R, 
        subtract.mean=F, 
        kern.interp.xs=kern.interp.xs,
        get_avg_kerns=T, 
        get_cross_kerns=T, 
        verbose=T)
    
    inv.coefs <- inversion$inv.coefs
    colnames(inv.coefs) <- paste0('c', 1:ncol(inv.coefs))
    nus <- cbind(ref.mod.tmp$nus, inv.coefs)
    
    radii <- ref.mod$r * ref.mod$R * solar_radius
    r_ts <- sapply(ref.mod.tmp$modes, function(mode) {
        ell <- as.numeric(strsplit(strsplit(mode, '_')[[1]][1], '\\.')[[1]][2])
        nn <- as.numeric(strsplit(strsplit(mode, '_')[[1]][2], '\\.')[[1]][2])
        nu <- ref.mod$nus[with(ref.mod$nus, l==ell&n==nn),]$nu.x
        #if (ell == 0) return(0)
        ref.mod$r[which.min((
                ref.mod$cs.spl(radii / (ref.mod$R * solar_radius))**2 / (radii**2) 
                - (10**-6*nu*(2*pi))**2/(ell+1/2)#ell*(ell+1))
            )**2)]
    })
    nus <- cbind(nus, data.frame(r_t = r_ts))
    
    kerns <- sapply(1:nrow(nus), 
        function(ii) max(abs(
            inv.coefs[ii, r0i] * ref.mod$k1[,rownames(inv.coefs)[ii]]
    )))
    kern.max <- max(kerns)
    
    print(modes)
    make_plots(plot_all, formatC(modes, width=3, format='d', flag=0), 
        filepath=file.path('plots', 'inv_anim', 'all'), 
        modes=modes,
        slides=F, make_pdf=F, thin=F, short=F,
        make_png=T, 
        paper_pdf_width=513.11743 / 72.27,
        paper_pdf_height=630 / 72.27,
        cex.paper=1.16) 
}


# ffmpeg -start_number 1  -framerate 3 -i tall/%03d.png animate.avi

