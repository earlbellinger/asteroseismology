#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('../../scripts/utils.R') 
library(Bolstad)
library(RColorBrewer)
library(pracma)
library(parallel)
library(parallelMap)

### MODELS
freqcols <- c('l', 'n', 'nu')
paths <- list(diffusion=file.path('diffusion_best', 'LOGS_MS'),
           no_diffusion=file.path('no_diffusion_best', 'LOGS_MS'),
                 modelS=file.path('..', 'modelS', 'ModelS_5b-freqs'),
                    MHD=file.path('..', 'modelS', 'JCD_MHD-freqs'),
                  muram=file.path('..', 'modelS', 'muram-freqs'),
                hl.diff=file.path('high_ell', 'diffusion'),
                hl.no_d=file.path('high_ell', 'no_diffusion'))

path <- paths$modelS
modelS <- list(name="Model S 5b", short='ModelS',
    kerns.dir=file.path(path),
    fgong=read.table(file.path(path, 'ModelS_5b.FGONG.dat'), header=1),
    freq=read.table(file.path(path, 'ModelS_5b-freqs.dat'), col.names=freqcols))

freqcols <- c('l', 'n', 'nu', 'E')
path <- paths$hl.diff
hl.diff <- list(name='MESA Diffusion', short='hlD',
    kerns.dir=file.path(path),
    #prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    freq=read.table(file.path(path, 'diffusion-freqs.dat'), header=1),
    fgong=read.table(file.path(path, 'diffusion.FGONG.dat'), header=1))

path <- paths$hl.no_d
hl.no_d <- list(name='MESA No Diffusion', short='hlnoD',
    kerns.dir=file.path(path),
    #prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    freq=read.table(file.path(path, 'no_diffusion-freqs.dat'), header=1),
    fgong=read.table(file.path(path, 'no_diffusion.FGONG.dat'), header=1))

path <- paths$diffusion
diffusion <- list(name='Diffusion', short='D',
    kerns.dir=file.path(path, 'profile1-freqs'),
    #prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    freq=read.table(file.path(path, 'profile1-freqs.dat'), col.names=freqcols),
    fgong=read.table(file.path(path, 'profile1-freqs', 
        'profile1.data.FGONG.dat'), header=1))

path <- paths$no_diffusion
no_diffusion <- list(name='No Diffusion', short='noD',
    kerns.dir=file.path(path, 'profile1-freqs'),
    prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    freq=read.table(file.path(path, 'profile1-freqs.dat'), col.names=freqcols),
    fgong=read.table(file.path(path, 'profile1-freqs', 
        'profile1.data.FGONG.dat'), header=1))

model.1 <- hl.diff #diffusion #hl.diff #no_diffusion_basu #no_diffusion
model.2 <- hl.no_d

### KERNEL PAIRS
c2_rho <- list(name='(c^2, rho)', f1='c2', f2='rho',
    f1.exp=bquote(c^2), f2.exp=bquote(rho))

u_Y <- list(name='(u, Y)', f1='u', f2='Y',
    f1.exp=bquote(u), f2.exp=bquote(Y))

u_Gamma1 <- list(name='(u, Gamma1)', f1='u', f2='Gamma1',
    f1.exp=bquote(u), f2.exp=bquote(Gamma[1]))

rho_Y <- list(name='(rho, Y)', f1='rho', f2='Y',
    f1.exp=bquote(rho), f2.exp=bquote(Y))

Gamma1_rho <- list(name='(Gamma1, rho)', f1='Gamma1', f2='rho',
    f1.exp=bquote(Gamma[1]), f2.exp=bquote(rho))

c2_Gamma1 <- list(name='(c^2, Gamma1)', f1='c2', f2='Gamma1',
    f1.exp=bquote(c^2), f2.exp=bquote(Gamma[1]))

k.pair <- c2_rho 
#k.pair <- rho_Y 
#k.pair <- u_Gamma1 
#k.pair <- u_Y 
#k.pair <- Gamma1_rho 
#k.pair <- c2_Gamma1 
for (k.pair in list(c2_rho, rho_Y, u_Gamma1, u_Y, Gamma1_rho, c2_Gamma1)) {


### LOAD MODEL
m.1 <- model.1$fgong
r <- rev(m.1$x)
m1.f1 <- rev(m.1[[k.pair$f1]])
m1.f2 <- rev(m.1[[k.pair$f2]])

if ('fgong' %in% names(model.2)) {
    m.2 <- model.2$fgong 
    
    m2.f1 <- splinefun(m.2$x, m.2[[k.pair$f1]])(r) 
    m2.f2 <- splinefun(m.2$x, m.2[[k.pair$f2]])(r) 
    
    # calculate relative structural differences 
    d.m1.f1 <- (m1.f1 - m2.f1) / (if (k.pair$f1 == 'Y') 1 else m1.f1)
    d.m1.f2 <- (m1.f2 - m2.f2) / (if (k.pair$f2 == 'Y') 1 else m1.f2)

    d.m2.f1 <- (m2.f1 - m1.f1) / (if (k.pair$f1 == 'Y') 1 else m2.f1)
    d.m2.f2 <- (m2.f2 - m1.f2) / (if (k.pair$f2 == 'Y') 1 else m2.f2)
}


### LOAD KERNELS
k1.fname <- paste0('E_K_', k.pair$f1, '-', k.pair$f2, '.dat')
k2.fname <- paste0('E_K_', k.pair$f2, '-', k.pair$f1, '.dat')
k.m1.f1 <- read.table(file.path(model.1$kerns.dir, k1.fname), header=1)
k.m1.f2 <- read.table(file.path(model.1$kerns.dir, k2.fname), header=1)
if (k.m1.f1$x[1] == 0) k.m1.f1 <- k.m1.f1[-1,] # trim 
if (k.m1.f2$x[1] == 0) k.m1.f2 <- k.m1.f2[-1,] # trim 


### LOAD FREQUENCIES
solar_data_dir <- file.path('..', '..', 'inverse', 'data')
sun <- rbind(read.table(file.path(solar_data_dir, 'SolarFreq_MDI.txt'), 
            col.names=c('l', 'n', 'nu.Sun', 'dnu')),
        read.table(file.path(solar_data_dir, 'Sun-freqs.dat'), header=1,
            col.names=c('n', 'l', 'nu.Sun', 'dnu')))
sun <- sun[order(sun$l, sun$n),]

#nus <- merge(sun, model.1$freq, by=c('l', 'n'))
ln <- c('l', 'n')
nus <- merge(merge(model.1$freq, model.2$freq, by=ln), sun, by=ln)
nus <- cbind(nus, data.frame(r.diff=with(nus, (nu.x-nu.y)/(nu.x)),
                         r.diff.Sun=with(nus, (nu.Sun-nu.x)/nu.Sun)))


ell_0 <- nus[nus$l==0,]
Q_0.x <- splinefun(ell_0$nu.x, ell_0$E.x)
Q_norm.x <- nus$E.x / Q_0.x(nus$nu.x)
nus <- cbind(nus, data.frame(m1.Q_norm=Q_norm.x))

nus <- nus[order(nus$l, nus$n),]

inertia <- nus$E.x #Q_norm.x #
nu <- nus$nu.x/5000
Xpinv <- ginv( matrix(c(nu**-2, nu**2)/inertia, ncol=2) )

a.r.1 <- Xpinv %*% nus$r.diff
m1.F_surf.r <- ( a.r.1[[1]]*nu**-2 + a.r.1[[2]]*nu**2 ) / inertia


plot_freq_diffs <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    low.l <- nus[nus$l <= 3,]
    corrected <- low.l$r.diff.Sun - low.l$m1.F_surf.r
    plot(low.l$nu.Sun, low.l$r.diff.Sun,
        axes=F, pch=low.l$l+1, col=blue, cex=0.5,
        ylim=range(0, low.l$r.diff.Sun, corrected),
        xlab=expression('Solar frequency' ~ nu/mu*Hz),
        ylab="")
    points(low.l$nu.Sun, corrected, col='darkred', pch=low.l$l+1, cex=0.5)
    abline(h=0, lty=2)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( 'Relative frequency difference' ~ delta*nu/nu ))
    legend("bottomleft", col=c(1,1,1,1, blue, 'darkred'), bty='n',
        inset=c(0.01, 0.01),
        pch=c(1:4, 20, 20), cex=text.cex, legend=as.expression(c(
            expression("\u2113"==0), expression("\u2113"==1), 
            expression("\u2113"==2), expression("\u2113"==3),
            bquote('Uncorrected'),
            bquote('Corrected'))))
}
plot_name <- paste0('solar_freq_diffs-', model.1$short)
print(plot_name)
make_plots(plot_freq_diffs, plot_name, filepath=file.path('plots', 'diffs'),
    paper_pdf_height=3/4*4.17309)


## calculate acoustic depth
int_0.r <- function(x, y) -as.numeric(cumtrapz(x, y)) # int_0^r 
int_r.R <- function(x, y) -trapz(x, y) - int_0.r(x, y) # int_r^R
tau <- int_r.R(m.1$r, 1/sqrt(m.1[['c2']])) # acoustic depth, Aerts et al eq 3.228
tau[tau<0] <- 0

num_nu_knots <- 30
num_r_knots <- 50

tau.points <- seq(max(tau), min(tau), length=num_r_knots-6)
r.knots.internal <- splinefun(tau, r)(tau.points)
r.knots <- c(r[1], r[1], r[1], r.knots.internal, rep(r[length(r)], 3))

nu.knots <- c(rep(min(nus$nu.x[nus$nu.x>0]), 3),
    seq(min(nus$nu.x), max(nus$nu.x), length=num_nu_knots-6),
    rep(max(nus$nu.x), 3))


## implement conservation of mass condition
im.mode <- 'l.-1_n.-1'
all.zeros <- 0 * 1:length(r)
k.m1.f1[[im.mode]] <- if (k.pair$f1 == 'rho') m1.f1 * r**2 else all.zeros
k.m1.f2[[im.mode]] <- if (k.pair$f2 == 'rho') m1.f2 * r**2 else all.zeros
nus <- plyr:::rbind.fill(nus, data.frame(l=-1, n=-1, nu.x=0, dnu=1, nu.y=0, 
    E.x=1, r.diff=0, m1.Q_norm=1, r.diff.Sun=0))

## get basis functions
parallelStartMulticore(12)#max(1,#as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
r.bfs <- do.call(rbind, 
    parallelMap(function(knot) sapply(r, function(x) B(x, knot, 3, r.knots)), 
        knot=1:(num_r_knots-4)))

## plot basis functions
plot_basis <- function(..., text.cex=1, mgp=utils.mgp, font="Times", 
        mar=utils.mar) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(0,0,2.2,0), cex.main=text.cex)
    plot(r, r.bfs[1,],
        xlim=c(0, 1),
        ylim=c(0, 1),
        axes=F, type='l', pch=20, 
        xlab=expression('Radius'~r/R),
        ylab=bquote( 'Cubic B-Spline Basis' ))
    for (ii in 2:nrow(r.bfs)) {
        lines(r, r.bfs[ii,], lty=ii)
    }
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0), 
            family=font, cex.axis=text.cex)
    #par(mgp=mgp+c(1.5, 0, 0))
    title(main=bquote( 'Acoustic depth' ~ tau/s ), line=2.02)
    axis(3, at=r.knots.internal, labels=F, cex.axis=text.cex, tcl=-0.1)
    indices <- c(which(r.knots.internal < 0.9)[c(T,F,F)], 
        length(r.knots.internal))
    axis(3, at=r.knots.internal[indices], labels=signif(tau.points[indices], 3),
        cex.axis=text.cex, tcl=-0.2)
}
plot_name <- paste0('basis_functions-', model.1$short)
print(plot_name)
make_plots(plot_basis, plot_name)



## calculate A matrix
get_row <- function(ell, nn, nu, inertia, Q_norm) {
    mode <- paste0('l.', ell, '_', 'n.', nn)
    if (!(mode %in% names(k.m1.f1) && mode %in% names(k.m1.f2))) return(NA)
    print(mode)
    
    #surfs <- do.call(cbind, 
    #    Map(function(knot) B(nu, knot, 3, nu.knots) / Q_norm, 
    #        knot=1:(num_nu_knots-4)))
    #names(surfs) <- paste0('a_', 1:length(surfs))
    
    surfs <- do.call(cbind,
         Map(function(pow) ifelse(nu != 0, (nu/5000)**pow/inertia, 0),
             pow=c(-2, 2)))
    names(surfs) <- paste0('a_', 1:length(surfs))
    
    #surfs <- data.frame(a_1=ifelse(nu != 0, (nu/5000)**-2/inertia, 0), 
    #                    a_2=(nu/5000)**2/inertia)
    
    f1.ints <- do.call(cbind, 
        Map(function(knot) 
                sintegral(r, r.bfs[knot,] * k.m1.f1[[mode]])$value,
                #trapz(r, r.bfs[knot,] * k.m1.f1[[mode]]), 
            knot=1:nrow(r.bfs)))
    names(f1.ints) <- paste0('b_', 1:length(f1.ints))
    
    f2.ints <- do.call(cbind, 
        Map(function(knot) 
                sintegral(r, r.bfs[knot,] * k.m1.f2[[mode]])$value, 
                #trapz(r, r.bfs[knot,] * k.m1.f2[[mode]]),
            knot=1:nrow(r.bfs)))
    names(f2.ints) <- paste0('c_', 1:length(f2.ints))
    
    data.frame(rbind(c(surfs, f1.ints, f2.ints)))
}

parallelStartMulticore(12)#max(1,#as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
ptm <- proc.time()
A.mat <- do.call(plyr:::rbind.fill, 
    parallelMap(get_row, ell=nus$l, nn=nus$n, 
        nu=nus$nu.x, inertia=nus$E.x, Q_norm=nus$m1.Q_norm))
proc.time() - ptm

## plot A matrix
plot_A_mat <- function(..., text.cex=1, mgp=utils.mgp, font="Times", 
        mar=utils.mar) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(0, 1, 0, 0))
    plot(NA, axes=F, 
        xlim=range(nus$nu.x),
        ylim=range(A.mat$b_34, A.mat$b_35, A.mat$b_36),
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="")
    abline(h=0, lty=2, col='gray', lwd=2)
    points( nus$nu.x, A.mat$b_34, pch=20, cex=0.3 )
    points( nus$nu.x, A.mat$b_35, pch=3, cex=0.3, col=red )
    points( nus$nu.x, A.mat$b_36, pch=4, cex=0.3, col=blue )
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1, 0, 0))
    title(ylab=bquote( integral( B[k](r) * K^(.(k.name)) * dr^"'" ) ))
}
k.name <- bquote(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))
plot_name <- paste0('Amat-', k.pair$f1, '_', k.pair$f2, '-', model.1$short)
print(plot_name)
make_plots(plot_A_mat, plot_name)


prof <- seq(min(r), max(r), length=100)
alpha <- 0.01
get_regularization_row <- function(r_j) {
    data.frame(alpha / sqrt(length(prof)) * do.call(cbind, 
        Map(function(knot) dB_dx(r_j, knot, 3, r.knots), 
            knot=1:(num_r_knots-4))))
}
r.block <- do.call(plyr:::rbind.fill, 
    parallelMap(get_regularization_row, r_j=prof))

nu.zeros <- zeros(length(prof), sum(grepl('a_', names(A.mat))))
r.zeros <- zeros(length(prof), num_r_knots-4)
f1.smooth <- cbind(nu.zeros, r.block, r.zeros)
names(f1.smooth) <- names(A.mat)
f2.smooth <- cbind(nu.zeros, r.zeros, r.block)
names(f2.smooth) <- names(A.mat)
R.mat <- rbind(f1.smooth, f2.smooth)


#svd.res <- svd(A.mat)# / nus$dnu)
A.mat.tot <- rbind(A.mat, R.mat)
d.vec.tot <- c(nus$r.diff, 0*1:nrow(R.mat))
svd.res <- svd(A.mat.tot)# / nus$dnu)
Sigma <- svd.res$d
Sigma.inv <- diag(length(svd.res$d)) * 1/Sigma
U <- svd.res$u
V <- svd.res$v
#res <- as.numeric(V %*% Sigma.inv %*% t(U) %*% (nus$r.diff))# / nus$dnu))
res <- as.numeric(V %*% Sigma.inv %*% t(U) %*% d.vec.tot)# / nus$dnu))


## calculate chi^2
resids <- rowSums(sweep(A.mat.tot, MARGIN=2, res, `*`))
chi2 <- sum((resids-d.vec.tot)**2)

plot_residuals <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(0, 1, 0, 0))
    
    plot(head(nus$nu.x, -1), head(nus$r.diff - resids[1:nrow(nus)], -1),
        cex=0.5, pch=20,
        #ylim=c(-0.05, 0.08),
        #xlim=c(0, 1),
        axes=F, 
        xlab=expression("Frequency"~nu/mu*Hz),
        ylab="")
    
    abline(h=0, lty=3)
    #points(prof, d.f1_d.r, col='darkred', pch=20, cex=0.5)
    
    #lines(r, d.m1.f1)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression('Residuals'~delta*nu/nu-Delta*nu))
}
plot_name <- paste0('residuals-', k.pair$f1, '_', k.pair$f2, '-', model.1$short)
print(plot_name)
make_plots(plot_residuals, plot_name)




## recover profiles
a.coefs <- res[grepl('a_', names(A.mat))]
b.coefs <- res[grepl('b_', names(A.mat))]
c.coefs <- res[grepl('b_', names(A.mat))]
d.f1_d.r <- sapply(prof, function(x) 
    sum(sapply(1:length(b.coefs), function(knot_i) 
        b.coefs[knot_i] * B(x, knot_i, 3, r.knots))))


plot_inversion <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    plot(prof, d.f1_d.r,
        type='l', lty=2, col='gray', lwd=2,
        ylim=c(-0.05, 0.08),
        xlim=c(0, 1),
        axes=F, 
        xlab="Radius r/R",
        ylab="")
    
    abline(h=0, lty=3)
    points(prof, d.f1_d.r, col='darkred', pch=20, cex=0.5)
    
    lines(r, d.m1.f1)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=expression('Relative sound speed difference' ~ delta*c^2/c^2))
}
plot_name <- paste0('inversion-', k.pair$f1, '_', k.pair$f2, '-', model.1$short)
print(plot_name)
make_plots(plot_inversion, plot_name)



plot(r.knots.internal, d.f1_d.r, ylim=c(-0.01, 0.01), type='l', lty=2, col='gray')
points(r.knots.internal, d.f1_d.r, pch=20, cex=0.5)

#A.mat <- NULL
#for (ell in unique(nus$l)) {
#    for (nn in sort(nus[nus$l == ell,]$n)) {
#        mode <- paste0('l.', ell, '_', 'n.', nn)
#        print(mode)
#        if (!(mode %in% names(k.m1.f1) && mode %in% names(k.m1.f2))) next
#        
#        k.1 <- splinefun(k.m1.f1$x, k.m1.f1[[mode]])(r)
#        k.2 <- splinefun(k.m1.f2$x, k.m1.f2[[mode]])(r)
#        
#        f1.ints <- c() 
#        f2.ints <- c() 
#        for (knot in 1:num_r_knots) {
#            spl <- sapply(r, function(x) N(x, knot, 3, r.knots))
#            f1.int <- k.1 * spl
#            f2.int <- k.2 * spl
#            f1.ints <- c(f1.ints, sintegral(r, f1.int)$value)
#            f2.ints <- c(f2.ints, sintegral(r, f2.int)$value)
#        }
#        
#        surfs <- c()
#        for (knot in 1:num_nu_knots) {
#            nu.row <- nus[nus$l==ell & nus$n==nn,]
#            spl <- N(nu.row$nu.x, knot, 3, nu.knots)
#            surfs <- c(surfs, spl)
#        }
#        
#        new.row <- data.frame(rbind(c(f1.ints, f2.ints, surfs)))
#        A.mat <- if (is.null(A.mat)) {
#            new.row 
#        } else plyr:::rbind.fill(A.mat, new.row)
#    }
#}

#P_1 <- function(x, ks, a) a[1] * ( x - ks[1] )**3
#P_2 <- function(x, ks, a) a[1] * ( ks[2] - ks[1] )**3 +
#        3 * a[1] * ( ks[2] - ks[1] )**2 * ( x - ks[2] ) +
#        6 * a[1] * ( ks[2] - ks[1] ) * ( x - ks[2] )**2 + 
#        a[3] * ( x - ks[2] )**3
#P_3 <- function(x, ks, a) a[2] * ( ks[4] - ks[5] )**3 +
#        3 * a[2] * ( ks[4] - ks[5] )**2 * ( x - ks[4] )**3 +
#        6 * a[2] * ( ks[4] - ks[5] ) * ( x - ks[4] )**2 +
#        a[4] * ( x - ks[4] )**3
#P_4 <- function(x, ks, a) a[2] * ( x - ks[5] )**3
#B_0 <- function(x, ks, a) {
#    if (x <= ks[1] || x >= ks[5]) 0
#    else if (x >= ks[1] && x <= ks[2]) P_1(x, ks, a)
#    else if (x >= ks[2] && x <= ks[3]) P_2(x, ks, a)
#    else if (x >= ks[3] && x <= ks[4]) P_3(x, ks, a)
#    else if (x >= ks[4] && x <= ks[5]) P_4(x, ks, a)
#}


}
