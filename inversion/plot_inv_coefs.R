#### Inversion coefficients vs. frequency 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')

col.pal <- c(blue, "#323031", "#F29559", "#DB4D48")

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
    get_avg_kerns=F, 
    get_cross_kerns=F, 
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

plot_inv.coefs <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,-0.4,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', xlab="", ylab="", 
        xlim=c(1500, 3200),
        ylim=c(-1.2, 1.2)) #c(-350, 350))
    
    with(nus, segments(nu.y, 0, nu.y, c1/max(abs(c1)), 
        col=col.pal[l+1], lwd=3))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    par(family="Helvetica")
    #text(x=2003, y=0.725, cex=0.9*text.cex, labels='2', col=col.pal[3])
    #text(x=2223, y=-0.68, cex=0.9*text.cex, labels='0', col=col.pal[1])
    #text(x=2580, y=-0.55, cex=0.9*text.cex, labels='3', col=col.pal[4])
    #text(x=2583, y=0.5, cex=0.9*text.cex, labels='1', col=col.pal[2])
    par(family=font)
    legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
        lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(Frequency~nu/mu*Hz))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Inv.~coef.))
}

make_plots(plot_inv.coefs, 'CygA_r005', 
    filepath=file.path('plots', 'invcoefs'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 

plot_inv.coefs_ltp <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,-0.4,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', xlab="", ylab="", 
        xlim=c(0, 0.2),
        ylim=c(-1.2, 1.2)) #c(-350, 350))
    
    abline(v=0.05, lwd=1.5, lty=2)
    
    with(nus, segments(r_t, 0, r_t, c1/max(abs(c1)), 
        col=col.pal[l+1], lwd=3))
    
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
    
    legend('bottomright', bty='n', legend=bquote(r[0] == 0.05),
        lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(Lower~turning~point~r[t]/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Inv.~coef.))
}

make_plots(plot_inv.coefs_ltp, 'CygA_ltp_r005', 
    filepath=file.path('plots', 'invcoefs'), 
    slides=F, make_png=F, 
    cex.paper=1.16) 









num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)
#set.seed(0)
n_iter <- length(inv.params$cross.terms) #3#
#inv.coefs.list <- list()
#for (ii in 1:n_iter) {
inv.coefs.list <- parallelMap(function(ii) {
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
        get_avg_kerns=F, 
        get_cross_kerns=F, 
        verbose=T)
    
    inversion$inv.coefs
    #inv.coefs.list[[ii]] <- inversion$inv.coefs
    #if (ii == 1) {
    #    inv.coefs <- inversion$inv.coefs
    #} else {
    #    inv.coefs <- inv.coefs + inversion$inv.coefs
    #}
}, ii=1:n_iter)
#inv.coefs <- apply(inv.coefs.list, 1, function(x) x[,1])

for (jj in 1:length(rs)) {
    inv.coefs.r <- Map(function(ii) inv.coefs.list[[ii]][,jj],
        ii=1:n_iter)
    inv.coefs.r <- data.frame(matrix(unlist(inv.coefs.r), 
        nrow=n_iter, byrow=T))
    inv.coefs <- apply(inv.coefs.r, 2, mean)
    inv.coefs.sd <- apply(inv.coefs.r, 2, sd)
    #colnames(inv.coefs) <- paste0('c', 1:ncol(inv.coefs))
    
    #inv.coefs <- inv.coefs / n_iter
    #colnames(inv.coefs) <- paste0('c', 1:ncol(inv.coefs))
    
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
    
    nus <- cbind(ref.mod$nus, inv.coefs, inv.coefs.sd)
    nus <- cbind(nus, data.frame(r_t = r_ts))
    
    
plot_inv.coefs <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,-0.4,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', xlab="", ylab="", 
        xlim=c(1500, 3200),
        ylim=c(-1.2, 1.2)) #c(-350, 350))
    
    maxabs <- max(abs(inv.coefs))
    with(nus, segments(nu.y, 0, nu.y, inv.coefs/maxabs, 
        col=col.pal[l+1], lwd=3))
    
    with(nus, segments(
        nu.y, 
        inv.coefs/maxabs + inv.coefs.sd/maxabs, 
        nu.y, 
        inv.coefs/maxabs - inv.coefs.sd/maxabs, 
        col=1, #adjustcolor(1, alpha.f=0.5), 
        lwd=0.5))#par()$lwd))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    if (jj == 1) {
        par(family="Helvetica")
        #text(x=1950, y=0.70, cex=0.9*text.cex, labels='2', col=col.pal[3])
        #text(x=2168, y=-0.7, cex=0.9*text.cex, labels='0', col=col.pal[1])
        #text(x=2580, y=-0.645, cex=0.9*text.cex, labels='3', col=col.pal[4])
        #text(x=2583, y=0.54, cex=0.9*text.cex, labels='1', col=col.pal[2])
        par(family=font)
    }
    legend('bottomright', bty='n', legend=bquote(r[0] == .(rs[jj])),
        lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(Frequency~nu/mu*Hz))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Inv.~coef.))
}

    make_plots(plot_inv.coefs, paste0('CygA_', jj), 
        filepath=file.path('plots', 'invcoefs'), 
        slides=F, make_png=F, 
        cex.paper=1.16) 
    
plot_inv.coefs_ltp <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,-0.4,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', xlab="", ylab="", 
        xlim=c(0, 0.4), 
        ylim=c(-1.2, 1.2)) #c(-350, 350))
    
    abline(v=rs[jj], lwd=par()$lwd, lty=2)
    
    maxabs <- max(abs(inv.coefs))
    with(nus, segments(r_t, 0, r_t, inv.coefs/maxabs, 
        col=col.pal[l+1], lwd=3))
    
    with(nus, segments(
        r_t, 
        inv.coefs/maxabs + inv.coefs.sd/maxabs, 
        r_t, 
        inv.coefs/maxabs - inv.coefs.sd/maxabs, 
        col=1, #adjustcolor(1, alpha.f=0.5), 
        lwd=0.5))#par()$lwd))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    box(lwd=par()$lwd)
    
    if (jj == 1) {
        par(family="Helvetica")
        #text(x=2223, y=-0.68, cex=0.9*text.cex, labels='0', col=col.pal[1])
        #text(x=0.035, y=0.55,   cex=0.9*text.cex, labels='1', col=col.pal[2])
        #text(x=0.078, y=0.65, cex=0.9*text.cex, labels='2', col=col.pal[3])
        #text(x=0.145, y=-0.2, cex=0.9*text.cex, labels='3', col=col.pal[4])
        par(family=font)
    }
    legend('bottomright', bty='n', legend=bquote(r[0] == .(rs[jj])),
        lty=NULL, pch=NULL, inset=c(-0.0075, -0.07), cex=text.cex)
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression(Lower~turning~point~r[t]/R))
    par(mgp=mgp+c(0.9, 0, 0))
    title(ylab=expression(Inv.~coef.))
}

    make_plots(plot_inv.coefs_ltp, paste0('CygA_ltp_', jj), 
        filepath=file.path('plots', 'invcoefs'), 
        slides=F, make_png=F, 
        cex.paper=1.16) 

}
