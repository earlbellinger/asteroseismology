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
library(optimx)

### MODELS
freqcols <- c('l', 'n', 'nu')
paths <- list(diffusion=file.path('diffusion_best', 'LOGS_MS'),
           no_diffusion=file.path('no_diffusion_best', 'LOGS_MS'),
                 modelS=file.path('..', 'modelS', 'ModelS_5b-freqs'),
                    MHD=file.path('..', 'modelS', 'JCD_MHD-freqs'),
                  muram=file.path('..', 'modelS', 'muram-freqs'),
                hl.diff=file.path('high_ell', 'diffusion'),
                hl.no_d=file.path('high_ell', 'no_diffusion'))

freqcols <- c('l', 'n', 'nu', 'E')
path <- paths$hl.diff
hl.diff <- list(name='MESA Diffusion', short='hlD',
    kerns.dir=file.path(path),
    freq=read.table(file.path(path, 'diffusion-freqs.dat'), header=1),
    fgong=read.table(file.path(path, 'diffusion.FGONG.dat'), header=1))

path <- paths$hl.no_d
hl.no_d <- list(name='MESA No Diffusion', short='hlnoD',
    kerns.dir=file.path(path),
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

model.1 <- diffusion
model.2 <- no_diffusion

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

k.pair <- u_Y 

### LOAD MODEL
m.1 <- model.1$fgong[nrow(model.1$fgong):1,]
r <- m.1$x
m1.f1 <- m.1[[k.pair$f1]]
m1.f2 <- m.1[[k.pair$f2]]

if ('fgong' %in% names(model.2)) {
    m.2 <- model.2$fgong[nrow(model.1$fgong):1,]
    
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
#if (k.m1.f1$x[1] == 0) k.m1.f1 <- k.m1.f1[-1,] # trim 
#if (k.m1.f2$x[1] == 0) k.m1.f2 <- k.m1.f2[-1,] # trim 

### LOAD FREQUENCIES
solar_data_dir <- file.path('..', '..', 'inverse', 'data')
sun <- rbind(read.table(file.path(solar_data_dir, 'SolarFreq_MDI.txt'), 
            col.names=c('l', 'n', 'nu.Sun', 'dnu')),
        read.table(file.path(solar_data_dir, 'Sun-freqs.dat'), header=1,
            col.names=c('n', 'l', 'nu.Sun', 'dnu')))
sun <- sun[order(sun$l, sun$n),]

ln <- c('l', 'n')
nus <- merge(merge(model.1$freq, model.2$freq, by=ln), sun, by=ln)
nus$nu.y <- rnorm(nrow(nus), nus$nu.y, nus$dnu)
nus <- cbind(nus, data.frame(r.diff=with(nus, (nu.x-nu.y)/(nu.x)),
                         r.diff.Sun=with(nus, (nu.Sun-nu.x)/nu.Sun)))

ell_0 <- nus[nus$l==0,]
Q_0.x <- splinefun(ell_0$nu.x, ell_0$E.x)
Q_norm.x <- nus$E.x / Q_0.x(nus$nu.x)
nus <- cbind(nus, data.frame(m1.Q_norm=Q_norm.x))
nus <- nus[order(nus$l, nus$n),]

modes <- paste0('l.', nus$l, '_', 'n.', nus$n)

## implement conservation of mass condition
if (k.pair$f1 == 'rho' || k.pair$f2 == 'rho') {
    im.mode <- 'l.-1_n.-1'
    #all.zeros <- 
    k.m1.f1[[im.mode]] <- if (k.pair$f1 == 'rho')
        4*pi*splinefun(r, m1.f1)(k.m1.f1$x)*k.m1.f1$x**2 else 0*1:nrow(k.m1.f1)
    k.m1.f2[[im.mode]] <- if (k.pair$f2 == 'rho') 
        4*pi*splinefun(r, m1.f2)(k.m1.f2$x)*k.m1.f2$x**2 else 0*1:nrow(k.m1.f2)
    nus <- plyr:::rbind.fill(nus, data.frame(l=-1, n=-1, nu.x=0, dnu=10**-6, 
        nu.y=0, E.x=1, r.diff=0, m1.Q_norm=1, r.diff.Sun=0))
}

parallelStartMulticore(16)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

get_square_Ks <- function(K=k.m1.f1, x0=NA) {
    Ks <- do.call(plyr:::rbind.fill, parallelMap(function(mode_i) {
        data.frame(do.call(cbind, Map(function(mode_j) {
            K_ij.integrand <- K[[mode_i]] * K[[mode_j]]
            MOLA <- if (!is.na(x0)) (1-x0)**2 * K_ij else 1
            MOLA * sintegral(K$x, K_ij.integrand)$value
        }, mode_j=modes)))
    }, mode_i=modes))
    names(Ks) <- modes
    rownames(Ks) <- modes
    Ks
}

K_ijs <- get_square_Ks(k.m1.f1, NA)
C_ijs <- get_square_Ks(k.m1.f2, NA)
K.ints <- do.call(c, parallelMap(function(mode) {
    sintegral(k.m1.f1$x, k.m1.f1[[mode]])$value
}, mode=modes))

get_A.mat <- function(beta, mu, nu.knots=NULL, use.BG=T, x0=NA) {
    
    first.ij <- do.call(plyr:::rbind.fill, parallelMap(function(mode_i, dnu) {
        data.frame(do.call(cbind, Map(function(mode_j) {
            E_ij <- ifelse(mode_i == mode_j, dnu**2, 0)
            K_ijs[mode_i, mode_j] + beta * C_ijs[mode_i, mode_j] + mu * E_ij
        }, mode_j=modes)))
    }, mode_i=modes, dnu=nus$dnu))
    names(first.ij) <- paste0('A_', 1:nrow(nus))
    
    mat <- rbind(cbind(first.ij, K.ints), c(K.ints, 0))
    
    surf <- do.call(plyr:::rbind.fill, 
        parallelMap(function(ell, nn, nu, E_i, Q_norm) {
            mode <- paste0('l.', ell, '_', 'n.', nn)
            print(mode)
            
            if (!(mode %in% names(k.m1.f1) && mode %in% names(k.m1.f2))) 
                return(NA)
            
            if (use.BG) {
                BG_surfs <- data.frame(do.call(cbind, Map(function(pow) 
                        ifelse(nu != 0, (nu/5000)**pow / E_i, 0), 
                    pow=c(-2, 2))))
                names(BG_surfs) <- paste0('a_', 1:length(BG_surfs)) 
            } else BG_surfs <- NULL 
            
            if (!is.null(nu.knots)) {
                F_surfs <- data.frame(do.call(cbind, Map(function(knot) 
                        B(nu, knot, deg, nu.knots) / Q_norm, 
                    knot=1:(length(nu.knots)-(deg+1)))))
                names(F_surfs) <- paste0('F_', 1:length(F_surfs))
            } else F_surfs <- NULL
            
            if (!is.null(BG_surfs) && !is.null(F_surfs)) {
                data.frame(rbind(c(BG_surfs, F_surfs)))
            } else if (!is.null(BG_surfs) && is.null(F_surfs)) {
                BG_surfs
            } else if (is.null(BG_surfs) && !is.null(F_surfs)) {
                F_surfs
            } else {
                NULL
            }
        }, ell=nus$l, nn=nus$n, nu=nus$nu.x, E_i=nus$E.x, Q_norm=nus$m1.Q_norm))
    
    mat <- cbind(mat, rbind(surf, 0))
    tsurf <- cbind(t(surf), zeros(ncol(surf), ncol(surf)+1))
    colnames(tsurf) <- names(mat)
    rbind(mat, tsurf)
}





solve(mat, system="LDLt", tol=0)


