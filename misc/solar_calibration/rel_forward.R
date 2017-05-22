#### Calculate the "relative forward problem" in asteroseismology  
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('../../scripts/utils.R') 
library(Bolstad)


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

# path <- paths$MHD
# MHD <- list(name="MHD", short='MHD',
    # kerns.dir=file.path(path),
    # fgong=read.table(file.path(path, 'JCD_MHD.FGONG.dat'), header=1),
    # freq=read.table(file.path(path, 'JCD_MHD-freqs.dat'), col.names=freqcols))

# path <- paths$muram
# muram <- list(name="MURaM", short='MURaM',
    # kerns.dir=file.path(path),
    # fgong=read.table(file.path(path, 'muram.FGONG.dat'), header=1),
    # freq=read.table(file.path(path, 'muram-freqs.dat'), col.names=freqcols))

# path <- paths$no_diffusion
# no_diffusion_basu <- list(name='MESA No Diffusion (Basu)', short='BnoD',
    # kerns.dir=file.path(path, 'basu'),
    # prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    # freq=read.table(file.path(path, 'basu', 'model.freqs'), col.names=freqcols),
    # fgong=read.table(file.path(path, 'profile1-freqs', 
        # 'profile1.data.FGONG.dat'), header=1))

# path <- paths$no_diffusion
# no_diffusion_ikit <- list(name='MESA No Diffusion (iKit)', short='noD',
    # kerns.dir='/scratch/seismo/bellinger/asteroseismology/inversions/compare/ik',
    # prof=read.table(file.path(path, 'profile1.data'), skip=5, header=1),
    # freq=read.table(file.path(path, 'profile1-freqs.dat'), col.names=freqcols),
    # fgong=read.table(file.path(path, 'profile1-freqs', 
        # 'profile1.data.FGONG.dat'), header=1))

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

model.1 <- diffusion #hl.diff #no_diffusion_basu #no_diffusion
model.2 <- no_diffusion #hl.no_d #modelS #MHD #diffusion #no_diffusion


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

#k.pair <- c2_rho 
#k.pair <- rho_Y 
#k.pair <- u_Gamma1 
#k.pair <- u_Y 
#k.pair <- Gamma1_rho 
#k.pair <- c2_Gamma1 
for (k.pair in list(u_Y, c2_rho, rho_Y, u_Gamma1, Gamma1_rho, c2_Gamma1)) {

### LOAD MODELS
m.1 <- model.1$fgong 
m.2 <- model.2$fgong 

r <- m.1$x

m1.f1 <- m.1[[k.pair$f1]]
m1.f2 <- m.1[[k.pair$f2]]

m2.f1 <- splinefun(m.2$x, m.2[[k.pair$f1]])(r) 
m2.f2 <- splinefun(m.2$x, m.2[[k.pair$f2]])(r) 

# calculate relative structural differences 
d.m1.f1 <- (m1.f1 - m2.f1) / (if (k.pair$f1 == 'Y') 1 else m1.f1)
d.m1.f2 <- (m1.f2 - m2.f2) / (if (k.pair$f2 == 'Y') 1 else m1.f2)

d.m2.f1 <- (m2.f1 - m1.f1) / (if (k.pair$f1 == 'Y') 1 else m2.f1)
d.m2.f2 <- (m2.f2 - m1.f2) / (if (k.pair$f2 == 'Y') 1 else m2.f2)

plot_profs <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot(r, m1.f1, axes=F, cex=0.1, pch=20, #type='l', #log='xy',
        xlim=c(0.9995, max(r)),
        ylim=range(m1.f1[r>=0.9995 & r<=max(r)]),
        xlab=expression('Radius'~r/R),
        ylab=expression(c))
    #points(r, m2.f1, cex=0.1, col='darkred')
    points(m.2$x, m.2[[k.pair$f1]], cex=0.1, pch=20, col='darkred')
    #legend("topright", lty=1:2, col=c('black', 'darkred'), cex=text.cex,
    #    legend=c('Profile file', 'FGONG'))
    abline(h=0, lty=2)
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0),
            las=0, mgp=mgp, family=font, cex.axis=text.cex)
}
#make_plots(plot_profs, 'atmosphere', filepath=file.path('plots', 'diffs'))

plot_diffs <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font, 
        mar=utils.mar, short=F) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex)
    if (short) par(mar=mar+c(0,0.1,0,0))
    plot(r, diffs, axes=F, type='l',
        lty=2,
        xlim=c( round(min(r[abs(diffs)>0.00005]), 1), 1 ),
        xlab=expression('Radius'~r/R),
        ylab="")#bquote( 'Relative difference in' ~ .(d.name) ))
        #ylab=bquote((.(d.name)[1] - .(d.name)[2]) / .(d.name)[1]))
    points(r, diffs, pch=20, cex=0.25)
    abline(h=0, lty=3)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    if (short) par(mgp=mgp+c(1, 0, 0))
    title(ylab=bquote( 'Difference in' ~ .(d.name) ))
}
diffs <- d.m1.f1
d.name <- k.pair$f1.exp
plot_name <- paste0('d_', k.pair$f1, '-', model.1$short, '_', model.2$short)
print(plot_name)
#make_plots(plot_diffs, plot_name, filepath=file.path('plots', 'diffs'))
diffs <- d.m1.f2
d.name <- k.pair$f2.exp
plot_name <- paste0('d_', k.pair$f2, '-', model.1$short, '_', model.2$short)
print(plot_name)
#make_plots(plot_diffs, plot_name, filepath=file.path('plots', 'diffs'))


### LOAD KERNELS
k1.fname <- paste0('E_K_', k.pair$f1, '-', k.pair$f2, '.dat')
k2.fname <- paste0('E_K_', k.pair$f2, '-', k.pair$f1, '.dat')
k.m1.f1 <- read.table(file.path(model.1$kerns.dir, k1.fname), header=1)
k.m1.f2 <- read.table(file.path(model.1$kerns.dir, k2.fname), header=1)
k.m2.f1 <- read.table(file.path(model.2$kerns.dir, k1.fname), header=1)
k.m2.f2 <- read.table(file.path(model.2$kerns.dir, k2.fname), header=1)
#k.m2.f1$x[1] <- min(modelS$fgong$x)
#k.m2.f2$x[1] <- min(modelS$fgong$x)


plot_kernels <- function(..., text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex)
    mode <- 'l.0_n.5'
    #k.1 <- splinefun(k.m1.f1$x, k.m1.f1[[mode]])(r)
    #k.2 <- splinefun(k.m2.f1$x, k.m2.f1[[mode]])(r)
    #max_inside <- max(abs(k[[mode]][k$x<0.9]), abs(k[['l.1_n.5']][k$x<0.9]))
    plot(k$x, k[[mode]],
        #k.1/max(k.1)-k.2/max(k.2), 
        #xlim=if (max_inside < 0.005) c(0.9, 1) else c(0, 1),
        xlim=c(round(min(k$x[ abs(k[[mode]])>0.005 ],
                         k$x[ abs(k[['l.1_n.5']])>0.005 ],
                         k$x[ abs(k[['l.2_n.5']])>0.005 ]), 1), 1),
        #ylim=range( k[[mode]] ),
        ylim=range( k[['l.0_n.5']], k[['l.1_n.5']], k[['l.2_n.5']] ),
        #ylim=if(any(k[[mode]] < -0.005)) c(-1, 1) else c(0, 2),
        axes=F, type='l', pch=20, #type='l', #log='xy',
        xlab=expression('Radius'~r/R),
        #log=if (max_inside < 0.005) 'x' else '',
        ylab=bquote( 'Kernel'~K^(.(k.name)) ))
    #lines(k.m1.f2$x, k.m1.f2[[mode]], col=red)
    #points(r, k.2, cex=0.1, pch=20, col='darkred')
    #legend("topright", lty=1:2, col=c('black', 'darkred'), cex=text.cex,
    #    legend=c('Profile file', 'FGONG'))
    lines(k$x, k[['l.1_n.5']], lty=2, col=blue)
    lines(k$x, k[['l.2_n.5']], lty=3, col=red)
    abline(h=0, lty=2)
    legend("topleft", lty=1:3, col=c('black', blue, red), cex=text.cex, bty='n',
        legend=c( 
            expression("\u2113" == 0*','~ n == 5), 
            expression("\u2113" == 1*','~ n == 5), 
            expression("\u2113" == 2*','~ n == 5)
        ))
    #magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0),
    #        las=1, mgp=mgp, family=font, cex.axis=text.cex)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
}
k <- k.m1.f1
k.name <- bquote(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))
plot_name <- paste0('kernel-', k.pair$f1, '_', k.pair$f2, '-', model.1$short)
print(plot_name)
#make_plots(plot_kernels, plot_name, filepath=file.path('plots', 'diffs'))
k <- k.m1.f2
k.name <- bquote(.(k.pair$f2.exp)*','~.(k.pair$f1.exp))
plot_name <- paste0('kernel-', k.pair$f2, '_', k.pair$f1, '-', model.1$short)
print(plot_name)
#make_plots(plot_kernels, plot_name, filepath=file.path('plots', 'diffs'))

plot_kernels2 <- function(..., text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex)
    mode <- 'l.2_n.3'
    #k.1 <- splinefun(k.m1.f1$x, k.m1.f1[[mode]])(r)
    #k.2 <- splinefun(k.m2.f1$x, k.m2.f1[[mode]])(r)
    #max_inside <- max(abs(k[[mode]][k$x<0.9]), abs(k[['l.2_n.9']][k$x<0.9]))
    plot(k$x, k[[mode]],
        #k.1/max(k.1)-k.2/max(k.2), 
        xlim=c(round(min(k$x[ abs(k[[mode]])>0.005 ],
                         k$x[ abs(k[['l.2_n.6']])>0.005 ], 
                         k$x[ abs(k[['l.2_n.9']])>0.005 ]), 1), 1),
        ylim=range( k[['l.2_n.3']], k[['l.2_n.6']], k[['l.2_n.9']] ),
        #xlim=if (max_inside < 0.005) c(0.9, 1) else c(0, 1),
        #ylim=if(any(k[[mode]] < -0.005)) c(-1, 1) else c(0, 2),
        axes=F, type='l', pch=20, #type='l', #log='xy',
        xlab=expression('Radius'~r/R),
        #log=if (max_inside < 0.005) 'x' else '',
        ylab=bquote( 'Kernel'~K^(.(k.name)) ))
    #lines(k.m1.f2$x, k.m1.f2[[mode]], col=red)
    #points(r, k.2, cex=0.1, pch=20, col='darkred')
    #legend("topright", lty=1:2, col=c('black', 'darkred'), cex=text.cex,
    #    legend=c('Profile file', 'FGONG'))
    lines(k$x, k[['l.2_n.6']], lty=2, col=blue)
    lines(k$x, k[['l.2_n.9']], lty=3, col=red)
    abline(h=0, lty=2)
    legend("topleft", lty=1:3, col=c('black', blue, red), cex=text.cex, bty='n',
        legend=c( 
            expression("\u2113" == 2*','~ n == 3), 
            expression("\u2113" == 2*','~ n == 6), 
            expression("\u2113" == 2*','~ n == 9)
        ))
    #magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0),
    #        las=1, mgp=mgp, family=font, cex.axis=text.cex)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
}
k <- k.m1.f1
k.name <- bquote(.(k.pair$f1.exp)*','~.(k.pair$f2.exp))
plot_name <- paste0('kernel2-', k.pair$f1, '_', k.pair$f2, '-', model.1$short)
print(plot_name)
#make_plots(plot_kernels2, plot_name, filepath=file.path('plots', 'diffs'))
k <- k.m1.f2
k.name <- bquote(.(k.pair$f2.exp)*','~.(k.pair$f1.exp))
plot_name <- paste0('kernel2-', k.pair$f2, '_', k.pair$f1, '-', model.1$short)
print(plot_name)
#make_plots(plot_kernels2, plot_name, filepath=file.path('plots', 'diffs'))
#plot_name <- paste0('ker_diffs-', 
#   k.pair$f1,     '_', k.pair$f2, '-',
#   model.1$short, '_', model.2$short)
#make_plots(plot_kernels, plot_name, filepath=file.path('plots', 'diffs'))


### LOAD FREQUENCIES
solar_data_dir <- file.path('..', '..', 'inverse', 'data')
sun <- rbind(read.table(file.path(solar_data_dir, 'SolarFreq_MDI.txt'), 
            col.names=c('l', 'n', 'nu.Sun', 'dnu')),
        read.table(file.path(solar_data_dir, 'Sun-freqs.dat'), header=1,
            col.names=c('n', 'l', 'nu.Sun', 'dnu')))
sun <- sun[order(sun$l, sun$n),]

ln <- c('l', 'n')
nus <- merge(merge(model.1$freq, model.2$freq, by=ln), sun, by=ln)
nus <- nus[order(nus$l, nus$n),]

#nus <- merge(model.1$freq, model.2$freq, by=c('l', 'n'))
#nus <- nus[nus$nu.x >= 1000 & nus$nu.y >= 1000 & 
#           nus$nu.x <= 4000 & nus$nu.y <= 4000,]

# adjust frequencies for mass and radius 
m.1.M <- max(model.1$fgong$m)
m.2.M <- max(model.2$fgong$m)

m.1.R <- with(model.1$fgong, splinefun(x, r)(1))
m.2.R <- with(model.2$fgong, splinefun(x, r)(1))

m1.scaler <- sqrt( m.1.M / m.1.R^3 )
m2.scaler <- sqrt( m.2.M / m.2.R^3 )

nus$m1.nu <- nus$nu.x / m1.scaler
nus$m2.nu <- nus$nu.y / m2.scaler

# calculate relative frequency differences directly 
nus <- cbind(nus, data.frame(
    #m1.r.diff=with(nus, (m1.nu-m2.nu)/(m1.nu)),
    m1.r.diff=with(nus, (m1.nu-m2.nu)/(m1.nu)),
    m2.r.diff=with(nus, (m2.nu-m1.nu)/(m2.nu))))

# calculate relative frequency differences through the kernel functions 
k.diffs <- NULL
for (ell in unique(nus$l)) {
    for (nn in nus[nus$l == ell,]$n) {
        mode <- paste0('l.', ell, '_', 'n.', nn)
        #print(mode)
        new.row <- data.frame(l=ell, n=nn)
        
        min.R <- min(r[r>0])# 0 # modelS$fgong$x[nrow(modelS$fgong)-1]
        if (mode %in% names(k.m1.f1) && mode %in% names(k.m1.f2)) {
            k.1 <- splinefun(k.m1.f1$x, k.m1.f1[[mode]])(r)
            k.2 <- splinefun(k.m1.f2$x, k.m1.f2[[mode]])(r)
            integrand <- k.1*d.m1.f1 + k.2*d.m1.f2
            k1.integrand <- k.1*d.m1.f1
            k2.integrand <- k.2*d.m1.f2
            m1.diff <- sintegral(r[r>=min.R], integrand[r>=min.R])$value
            m1.k1.diff <- sintegral(r[r>=min.R], k1.integrand[r>=min.R])$value
            m1.k2.diff <- sintegral(r[r>=min.R], k2.integrand[r>=min.R])$value
            new.row <- cbind(new.row, data.frame(m1.K.diff=m1.diff,
                m1.k1.diff=m1.k1.diff, m1.k2.diff=m1.k2.diff))
        }
        
        if (mode %in% names(k.m2.f1) && mode %in% names(k.m2.f2)) {
            k.1 <- splinefun(k.m2.f1$x, k.m2.f1[[mode]])(r)
            k.2 <- splinefun(k.m2.f2$x, k.m2.f2[[mode]])(r)
            integrand <- k.1*d.m2.f1 + k.2*d.m2.f2
            k1.integrand <- k.1*d.m2.f1
            k2.integrand <- k.2*d.m2.f2
            m2.diff <- sintegral(r[r>=min.R], integrand[r>=min.R])$value
            m2.k1.diff <- sintegral(r[r>=min.R], k1.integrand[r>=min.R])$value
            m2.k2.diff <- sintegral(r[r>=min.R], k2.integrand[r>=min.R])$value
            new.row <- cbind(new.row, data.frame(m2.K.diff=m2.diff,
                m2.k1.diff=m2.k1.diff, m2.k2.diff=m2.k2.diff))
        }
        
        k.diffs <- if (is.null(k.diffs)) {
            new.row 
        } else plyr:::rbind.fill(k.diffs, new.row)
   }
}
nus <- merge(nus, k.diffs, by=c('l', 'n'))
#nus <- nus[complete.cases(nus),]

# calculate surface term
if ('E' %in% names(model.1$freq)) {
    ell_0 <- nus[nus$l==0,]
    Q_0.x <- splinefun(ell_0$m1.nu, ell_0$E.x)
    Q_norm.x <- nus$E.x / Q_0.x(nus$m1.nu)
    
    inertia <- nus$E.x #Q_norm.x #
    #nu <- nus$m1.nu#/(5000*m1.scaler)
    nu <- nus$nu.x / 5000
    Xpinv <- ginv( matrix(c(nu**-2, nu**2)/inertia, ncol=2) )
    
    a.r.1 <- Xpinv %*% ( nus$m1.K.diff - nus$m1.r.diff )
    m1.F_surf <- ( a.r.1[[1]]*nu**-2 + a.r.1[[2]]*nu**2 ) / inertia
    
    #a.K.1 <- Xpinv %*% nus$m1.K.diff
    #m1.F_surf.K <- ( a.K.1[[1]]*nu**-2 + a.K.1[[2]]*nu**2 ) / inertia
    
    nus <- cbind(nus, data.frame(
        m1.Q_norm=Q_norm.x,
        m1.F_surf=m1.F_surf,
        m1.diff=nus$m1.r.diff - nus$m1.K.diff + m1.F_surf
        #m1.F_surf.K=m1.F_surf.K,
        #m1.diff=(nus$m1.K.diff-m1.F_surf.K)-(nus$m1.r.diff-m1.F_surf.r)
    ))
}

if ('E' %in% names(model.2$freq)) {
    ell_0 <- nus[nus$l==0,]
    Q_0.y <- splinefun(ell_0$m2.nu, ell_0$E.y)
    Q_norm.y <- nus$E.y / Q_0.y(nus$m2.nu)
    
    inertia <- nus$E.y #Q_norm.y #
    nu <- nus$m1.nu/(5000*m1.scaler)
    Xpinv <- ginv( matrix(c(nu**-2, nu**2)/inertia, ncol=2) )
    
    a.r.2 <- Xpinv %*% nus$m1.r.diff
    m2.F_surf.r <- ( a.r.2[[1]]*nu**-2 + a.r.2[[2]]*nu**2 ) / inertia
    
    a.K.2 <- Xpinv %*% nus$m1.K.diff
    m2.F_surf.K <- ( a.K.2[[1]]*nu**-2 + a.K.2[[2]]*nu**2 ) / inertia
    
    nus <- cbind(nus, data.frame(
        m2.Q_norm=Q_norm.y,
        m2.F_surf.r=m2.F_surf.r,
        m2.F_surf.K=m2.F_surf.K,
        m2.diff=(nus$m2.K.diff-m2.F_surf.K)-(nus$m2.r.diff-m2.F_surf.r)
    ))
}

## save result
write.table(nus, file.path('nus', paste0('nus-', 
        k.pair$f1,     '_', k.pair$f2, '-',
        model.1$short, '_', model.2$short, '.dat')),
    quote=F, row.names=F)


### PLOT RESULT
plot_rel_diffs <- function(...,
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    #with(new.nus, 
    #plot(nu.x, (nu.x-nu.y)/(nu.x) - rdiff), axes=F,
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab=bquote( 'Relative frequency difference' ~ delta*nu/nu ), 
        xlim=range(nus$nu.x, nus$nu.y),
        #ylim=range(nus[,-1:-6]) * c(1.8, 1.3))
        ylim=c(-0.0075, 0.0075))
        #ylim=c(-2, 2))
        #ylim=c(-0.003, 0.003))
    abline(h=0, lty=2)
    points(nus$nu.x, nus$m1.K.diff, col=blue, cex=0.2, pch=20)#pch=nus$l+1)
    points(nus$nu.x, nus$m1.r.diff, col='darkred', cex=0.2, pch=20)# pch=nus$l+1)
    points(nus$nu.y, nus$m2.r.diff, col='black', cex=0.2, pch=20)# pch=nus$l+1)
    points(nus$nu.y, nus$m2.K.diff, col='purple', cex=0.2, pch=20)# pch=nus$l+1)
    legend('topleft', cex=text.cex, 
        legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    legend("bottomleft", col=c(blue, 'darkred', 'black', 'purple'), 
        pch=1, cex=text.cex, legend=as.expression(c(
            bquote(.(model.1$name)~'Kernels'),
            bquote('Differences w.r.t.'~.(model.1$name)), 
            bquote('Differences w.r.t.'~.(model.2$name)),
            bquote(.(model.2$name)~'Kernels'))))
    #legend("bottomright", pch=1:4, legend=c('l=0', 'l=1', 'l=2', 'l=3'),
    #    cex=text.cex)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=0, cex.axis=text.cex)
}

plot_name <- paste0('rel_diffs-', 
   k.pair$f1,     '_', k.pair$f2, '-',
   model.1$short, '_', model.2$short)
print(plot_name)
make_plots(plot_rel_diffs, plot_name)


### PLOT RESULT
plot_one_surfless <- function(...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    #with(new.nus, 
    #plot(nu.x, (nu.x-nu.y)/(nu.x) - rdiff), axes=F,
    Kdiff <- nus$m1.K.diff
    rdiff <- nus$m1.r.diff
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", 
        xlim=range(nus$nu.x, nus$nu.y),
        ylim=range(Kdiff, rdiff))
        #ylim=range(nus[,-1:-6]) * c(1.8, 1.3))
        #ylim=c(-0.0075, 0.0075))
        #ylim=c(-2, 2))
        #ylim=c(-0.003, 0.003))
    abline(h=0, lty=2)
    segments(nus$nu.x, Kdiff, nus$nu.x, rdiff, 
        col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    points(nus$nu.x, Kdiff, col=blue, cex=0.33, pch=pch)
    points(nus$nu.x, rdiff, col='darkred', cex=0.33, pch=pch)
    #legend('topleft', cex=text.cex, 
    #    legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    legend("topright", col=c(1,1,1,1,1, blue, 'darkred'), bty='n',
        inset=c(0.01, 0.01),
        pch=c(NA, 1:4, 20, 20), cex=text.cex, legend=as.expression(c(
            #bquote(.(model.1$name)~'Kernels'),
            bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            expression("\u2113"==0), expression("\u2113"==1), 
            expression("\u2113"==2), expression("\u2113"==3),
            bquote('Kernels'),
            bquote('Exact'))))
    #legend('topright', pch=c(NA, 1:4), cex=text.cex, legend=c(
    #    bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
    #    expression(l==0), expression(l==1), expression(l==2), expression(l==3)))
    #magaxis(side=1:4, family=font, tcl=-0.25, labels=c(1,1,0,0), mgp=mgp, 
    #    las=1, cex.axis=text.cex)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( 'Relative frequency difference' ~ delta*nu/nu ))
}

plot_name <- paste0('rel_diffs_one_surfless-', 
   k.pair$f1,     '_', k.pair$f2, '-',
   model.1$short, '_', model.2$short)
print(plot_name)
make_plots(plot_one_surfless, plot_name)



### PLOT RESULT
plot_one <- function(...,
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    #with(new.nus, 
    #plot(nu.x, (nu.x-nu.y)/(nu.x) - rdiff), axes=F,
    Kdiff <- nus$m1.K.diff - nus$m1.F_surf
    rdiff <- nus$m1.r.diff #- nus$m1.F_surf #- nus$m1.F_surf.r
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab="", #bquote( 'Relative frequency difference' ~ delta*nu/nu ), 
        xlim=range(nus$nu.x, nus$nu.y),
        ylim=range(Kdiff, rdiff))
        #ylim=range(nus[,-1:-6]) * c(1.8, 1.3))
        #ylim=c(-0.0075, 0.0075))
        #ylim=c(-2, 2))
        #ylim=c(-0.003, 0.003))
    abline(h=0, lty=2)
    segments(nus$nu.x, Kdiff, nus$nu.x, rdiff, 
        col=adjustcolor('black',alpha.f=0.25), lty=3)
    pch <- if (4 %in% nus$l) 20 else (nus$l+1)
    points(nus$nu.x, Kdiff, col=blue, cex=0.33, pch=pch)
    points(nus$nu.x, rdiff, col='darkred', cex=0.33, pch=pch)
    #legend('topleft', cex=text.cex, 
    #    legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    legend("topright", col=c(1,1,1,1,1, blue, 'darkred'), bty='n',
        inset=c(0.01, 0.01),
        pch=c(NA, 1:4, 20, 20), cex=text.cex, legend=as.expression(c(
            #bquote(.(model.1$name)~'Kernels'),
            bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
            expression("\u2113"==0), expression("\u2113"==1), 
            expression("\u2113"==2), expression("\u2113"==3),
            bquote('Kernels'),
            bquote('Exact'))))
    #legend("topright", col=c(1,1,1,1, blue, 'darkred'), bty='n',
    #    inset=c(0.01, 0.01),
    #    pch=c(1:4, 20, 20), cex=text.cex, legend=as.expression(c(
    #        #bquote(.(model.1$name)~'Kernels'),
    #        #bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))),
    #        expression("\u2113"==0), expression("\u2113"==1), 
    #        expression("\u2113"==2), expression("\u2113"==3),
    #        bquote('Kernels'),
    #        bquote('Exact'))))
    #legend('top', pch=NA, cex=text.cex, bty='n', inset=c(0.01, 0.01), 
    #    legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    #    expression(l==0), expression(l==1), expression(l==2), expression(l==3)))
    #magaxis(side=1:4, family=font, tcl=-0.25, labels=c(1,1,0,0), mgp=mgp, 
    #    las=0, cex.axis=text.cex)
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=1, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote( 'Relative frequency difference' ~ delta*nu/nu ))
}

plot_name <- paste0('rel_diffs_one-', 
   k.pair$f1,     '_', k.pair$f2, '-',
   model.1$short, '_', model.2$short)
print(plot_name)
make_plots(plot_one, plot_name)



### PLOT RESULT
plot_diff_diffs <- function(...,
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    Kdiff <- nus$m1.K.diff - nus$m1.F_surf
    rdiff <- nus$m1.r.diff #- nus$m1.F_surf.r
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab=bquote( 'Kernel Rel. Diff. Minus Actual Rel. Diff.' ), 
        xlim=range(nus$nu.x, nus$nu.y),
        ylim=range(Kdiff - rdiff))
    abline(h=0, lty=2)
    #pch <- if (4 %in% nus$l) 20 else nus$l+1
    points(nus$nu.x, Kdiff - rdiff, col=blue, cex=0.5, pch=20)
    for (ell_i in 0:4) {
        ell <- nus[nus$l == ell_i,]
        Kdiff <- ell$m1.K.diff - ell$m1.F_surf
        rdiff <- ell$m1.r.diff #- ell$m1.F_surf.r
        points(ell$nu.x, Kdiff - rdiff, col=ell_i+1, cex=0.5, pch=1)
    }
    legend('topleft', cex=text.cex, col=blue, pch=20,
        legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=0, cex.axis=text.cex)
}

plot_name <- paste0('diff_diffs-', 
   k.pair$f1,     '_', k.pair$f2, '-',
   model.1$short, '_', model.2$short)
print(plot_name)
make_plots(plot_diff_diffs, plot_name)




plot_components <- function(...,
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    #with(new.nus, 
    #plot(nu.x, (nu.x-nu.y)/(nu.x) - rdiff), axes=F,
    plot(NA, axes=F,
        xlab=expression('Frequency'~nu/mu*Hz),
        ylab=bquote( delta*nu/nu ), col='black', pch=nus$l+1,
        xlim=range(nus$nu.x, nus$nu.y),
        #ylim=c(-5, 0))
        ylim=range(nus[,-1:-6]) * c(1.8, 1.3))
        #ylim=c(-0.0075, 0.0075))
        #ylim=c(-2, 2))
        #ylim=c(-0.003, 0.003))
    abline(h=0, lty=2)
    
    points(nus$nu.x, nus$m1.k1.diff, col=blue, cex=0.2, pch=20)#, pch=nus$l+1)
    points(nus$nu.x, nus$m1.k2.diff, col='darkred', cex=0.2, pch=20)#, pch=nus$l+1)
    points(nus$nu.y, nus$m2.k1.diff, col='black', cex=0.2, pch=20)#, pch=nus$l+1)
    points(nus$nu.y, nus$m2.k2.diff, col='purple', cex=0.2, pch=20)#, pch=nus$l+1)
    
    legend("bottomleft", col=c(blue, 'darkred', 'black', 'purple'), 
        pch=1, cex=text.cex, legend=as.expression(c(
            bquote(.(model.1$name)~'u Kernel'),
            bquote(.(model.1$name)~'Y Kernel'), 
            bquote(.(model.2$name)~'u Kernel'),
            bquote(.(model.2$name)~'Y Kernel'))))
    legend('topleft', cex=text.cex, 
        legend=bquote((.(k.pair$f1.exp)*','~.(k.pair$f2.exp))))
    #legend("bottomright", pch=1:4, legend=c('l=0', 'l=1', 'l=2', 'l=3'),
    #    cex=text.cex)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp, 
        las=0, cex.axis=text.cex)
}

plot_name <- paste0('comps-', 
   k.pair$f1,     '_', k.pair$f2, '-',
   model.1$short, '_', model.2$short)
print(plot_name)
#make_plots(plot_components, plot_name)

#plot_rel_diffs()
#dev.off()


}

