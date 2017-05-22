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

model.1 <- hl.diff
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

## calculate acoustic depth, Aerts et al eq 3.228 
int_0.r <- function(x, y) as.numeric(cumtrapz(x, y)) # int_0^r 
int_r.R <- function(x, y) trapz(x, y) - int_0.r(x, y) # int_r^R 
tau <- int_r.R(m.1$r, 1/sqrt(m.1[['c2']])) 
tau[tau<0] <- 0 

## implement conservation of mass condition
im.mode <- 'l.-1_n.-1'
#all.zeros <- 
k.m1.f1[[im.mode]] <- if (k.pair$f1 == 'rho')
    4*pi*splinefun(r, m1.f1)(k.m1.f1$x) * k.m1.f1$x**2 else 0*1:nrow(k.m1.f1)
k.m1.f2[[im.mode]] <- if (k.pair$f2 == 'rho') 
    4*pi*splinefun(r, m1.f2)(k.m1.f2$x) * k.m1.f2$x**2 else 0*1:nrow(k.m1.f2)
nus <- plyr:::rbind.fill(nus, data.frame(l=-1, n=-1, nu.x=0, dnu=10**-6, 
    nu.y=0, E.x=1, r.diff=0, m1.Q_norm=1, r.diff.Sun=0))


parallelStartMulticore(16)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

get_r.knots.tau <- function(r, tau, num_r_knots, degree=4) {
    #n <- num_r_knots-(2*degree+1)
    #inside <- tau[2:(length(tau)-1)]
    #taus <- sort(splinefun(tau, m.1$x)(seq(min(inside), max(inside), 
    #    length=n)))
    #c(rep(min(r), degree+1), taus, rep(max(taus), degree))
    n <- num_r_knots-(2*degree)
    #min.tau <- splinefun(m.1$x, tau)(1)
    #taus <- sort(splinefun(tau, m.1$x)(seq(min.tau, max(tau), length=n)))
    #min.tau <- splinefun(m.1$x, tau)(1)
    int.r <- r[2:(length(r)-1)]
    taus <- splinefun(m.1$x, tau)(int.r)
    taus <- sort(splinefun(taus, int.r)(seq(min(taus), max(taus), length=n)))
    c(rep(0, degree), taus, rep(r[length(r)-1], degree))
    #c(rep(0, degree), taus, rep(1, degree))
}

get_knots <- function(x, num_knots, degree=4) {
    x <- x[x>0]
    n <- num_knots-(2*degree)
    c(rep(min(x), degree), 
      seq(min(x), max(x), length=n),
      rep(max(x), degree))
}

get_r.bfs <- function(r, r.knots, degree=4) {
    do.call(rbind, parallelMap(function(knot_i) {
        sapply(r, function(x) B(x, knot_i, degree, r.knots))
    }, knot_i=1:(length(r.knots)-(degree+1))))
}

#get_basis_functions <- function(x, knots, degree=3) {
#    do.call(rbind, parallelMap(function(knot_i) {
#        sapply(x, function(x) B(x, knot_i, degree, knots))
#    }, knot_i=1:(length(knots)-(degree+1))))
#}

get_A.mat <- function(r, r.bfs, nus=nus, k.m1.f1=k.m1.f1, k.m1.f2=k.m1.f2,
        nu.knots=NULL, use.BG=T) {
    if (!is.null(nu.knots)) deg <- (length(nu.knots)-length(unique(nu.knots)))/2
    do.call(plyr:::rbind.fill, parallelMap(function(ell, nn, nu, E_i, Q_norm) {
        mode <- paste0('l.', ell, '_', 'n.', nn)
        print(mode)
        
        if (!(mode %in% names(k.m1.f1) && mode %in% names(k.m1.f2))) 
            return(NA)
        
        if (use.BG) {
            BG_surfs <- do.call(cbind, Map(function(pow) 
                    ifelse(nu != 0, (nu/5000)**pow / E_i, 0),
                pow=c(-2, 2)))
            names(BG_surfs) <- paste0('a_', 1:length(BG_surfs))
        } else BG_surfs <- NULL
        
        if (!is.null(nu.knots)) {
            F_surfs <- do.call(cbind, Map(function(knot) 
                    B(nu, knot, deg, nu.knots) / Q_norm, 
                knot=1:(length(nu.knots)-(deg+1))))
            names(F_surfs) <- paste0('F_', 1:length(F_surfs))
        } else F_surfs <- NULL
        
        f1.ints <- do.call(cbind, 
            Map(function(knot) 
                    sintegral(r, r.bfs[knot,] * k.m1.f1[[mode]])$value,
                knot=1:nrow(r.bfs)))
        names(f1.ints) <- paste0('b_', 1:length(f1.ints))
        
        f2.ints <- do.call(cbind, 
            Map(function(knot) 
                    sintegral(r, r.bfs[knot,] * k.m1.f2[[mode]])$value,
                knot=1:nrow(r.bfs)))
        names(f2.ints) <- paste0('c_', 1:length(f2.ints))
        
        if (!is.null(BG_surfs) && !is.null(F_surfs)) {
            data.frame(rbind(c(BG_surfs, f1.ints, f2.ints, F_surfs)))
        } else if (!is.null(BG_surfs) && is.null(F_surfs)) {
            data.frame(rbind(c(BG_surfs, f1.ints, f2.ints)))
        } else if (is.null(BG_surfs) && !is.null(F_surfs)) {
            data.frame(rbind(c(f1.ints, f2.ints, F_surfs))) 
        } else {
            data.frame(rbind(c(f1.ints, f2.ints)))
        }
    }, ell=nus$l, nn=nus$n, nu=nus$nu.x, E_i=nus$E.x, Q_norm=nus$m1.Q_norm))
}

get_R.mat <- function(alpha, prof, r.knots, coef.names, p.degree) {
    num_knots <- sum(grepl('b_', coef.names))
    r.block <- alpha / sqrt(length(prof)) * do.call(plyr:::rbind.fill, 
        parallelMap(function(r_j) {
            data.frame(do.call(cbind, 
                Map(function(knot_i) 
                    d.n_B_dx.n(r_j, knot_i, p.degree, r.knots, n=2), 
                knot_i=1:num_knots)))
        }, r_j=prof))
    
    nu.zeros <- zeros(length(prof), sum(grepl('a_', coef.names)))
    r.zeros <- zeros(length(prof), num_knots)
    
    f1.smooth <- cbind(nu.zeros, r.block, r.zeros)
    names(f1.smooth) <- coef.names
    f2.smooth <- cbind(nu.zeros, r.zeros, r.block)
    names(f2.smooth) <- coef.names
    rbind(f1.smooth, f2.smooth)
    
    #R.mat <- cbind(nu.zeros, r.block, r.block)
    #names(R.mat) <- coef.names
    #R.mat
}

get_R.p <- function(alpha.f1, alpha.f2, alpha.nu, 
        coef.names, r.knots, p.deg.r=2, p.deg.nu=1) {
    
    num_BG <- sum(grepl('a_', coef.names))
    num_r_knots <- sum(grepl('b_', coef.names))
    num_nu_knots <- sum(grepl('F_', coef.names))
    
    p.rows.r <- num_r_knots-p.deg.r
    blank <- data.frame(row.names=1:p.rows.r)
    BG.zeros <- if (num_BG) zeros(p.rows.r, num_BG) else blank
    nu.zeros <- if (num_nu_knots) zeros(p.rows.r, num_nu_knots) else blank
    r.zeros <- zeros(p.rows.r, num_r_knots)
    r.penalties <- diff(diag(num_r_knots), diff=p.deg.r)
    
    #smooth.prof <- p.deg.r / 
    #    (r.knots[(length(r.knots)-p.rows.r+1) : length(r.knots)] + .03)
    #r.penalties <- r.penalties * smooth.prof 
    
    f1.smooth <- as.data.frame(cbind(BG.zeros, r.penalties, r.zeros, nu.zeros))
    names(f1.smooth) <- coef.names
    f2.smooth <- as.data.frame(cbind(BG.zeros, r.zeros, r.penalties, nu.zeros))
    names(f2.smooth) <- coef.names
    
    nu.smooth <- if (num_nu_knots > 2) {
        p.rows.nu <- num_nu_knots-p.deg.nu
        blank <- data.frame(row.names=1:p.rows.nu)
        BG.zeros <- if (num_BG) zeros(p.rows.nu, num_BG) else blank
        nu.zeros <- if (num_nu_knots) zeros(p.rows.nu, num_nu_knots) else blank
        r.zeros <- zeros(p.rows.nu, num_r_knots)
        nu.penalties <- diff(diag(num_nu_knots), diff=p.deg.nu)
        nu.smooth <- as.data.frame(cbind(BG.zeros, 
            r.zeros, r.zeros, nu.penalties))
        names(nu.smooth) <- coef.names
        nu.smooth
    } else data.frame()
    
    
    #sweep(r.penalties, MARGIN=2, smooth.prof, `*`)
    
    rbind(sqrt(alpha.f1) * f1.smooth, 
          sqrt(alpha.f2) * f2.smooth, 
          sqrt(alpha.nu) * nu.smooth)
}

get_A.tot <- function(A, nus, R) {
    rbind(A / nus$dnu, R)
    #rbind(A, R)
}

get_d.tot <- function(nus, R) {
    c(nus$r.diff / nus$dnu, 0*1:nrow(R))
    #c(nus$r.diff, 0*1:nrow(R))
}

get_inversion_coefs <- function(A, d) {
    svd.res <- svd(A)
    Sigma.inv <- diag(length(svd.res$d)) * 1/svd.res$d
    x <- as.numeric(svd.res$v %*% Sigma.inv %*% t(svd.res$u) %*% d)
    #x
    e <- sapply(1:ncol(A), function(k) 
        sum( svd.res$v[k,]**2 * svd.res$u[k,]**2 
           / svd.res$d**2 )) 
    list(x=x, e=e)
}

boot.coefs <- function(tmp, indices) {
    print(length(unique(indices)))
    A.data <- A.mat[indices,]
    nus.data <- nus[indices,]
    A <- get_A.tot(A.data, nus, R.mat)
    d <- get_d.tot(nus.data, R.mat)
    svd.res <- svd(A)
    Sigma.inv <- diag(length(svd.res$d)) * 1/svd.res$d
    x <- as.numeric(svd.res$v %*% Sigma.inv %*% t(svd.res$u) %*% d)
    x
}

gcv <- function(alpha, A, nus, degree) {
    R <- get_R.p(alpha, names(A), degree)
    d <- get_d.tot(nus, R)
    A.tot <- as.matrix(get_A.tot(A, nus, R))
    svd.res <- svd(A.tot)
    Sigma.inv <- diag(length(svd.res$d)) * 1/svd.res$d
    xTx.inv.x <- svd.res$v %*% Sigma.inv %*% t(svd.res$u)
    res <- as.numeric(xTx.inv.x %*% d)
    y_hat <- rowSums(sweep(A.tot[1:nrow(nus),], MARGIN=2, res, `*`))
    rss <- sum((y_hat - d[1:nrow(nus)])**2)
    tr.H <- sum(diag( 1 - A.tot %*% xTx.inv.x ))
    ( rss / nrow(nus) ) / ( tr.H/nrow(nus) )**2
}

gcv2 <- function(alpha, A, nus, p.deg.r=1, p.deg.nu=1) {
    R <- get_R.p(alpha[1], alpha[2], alpha[3], names(A), 
        r.knots, p.deg.r, p.deg.nu)
    d <- get_d.tot(nus, R)
    A.tot <- as.matrix(get_A.tot(A, nus, R))
    svd.res <- svd(A.tot)
    Sigma.inv <- diag(length(svd.res$d)) * 1/svd.res$d
    xTx.inv.x <- svd.res$v %*% Sigma.inv %*% t(svd.res$u)
    res <- as.numeric(xTx.inv.x %*% d)
    y_hat <- rowSums(sweep(A.tot[1:nrow(nus),], MARGIN=2, res, `*`))
    rss <- sum((y_hat - d[1:nrow(nus)])**2)
    tr.H <- sum(diag( 1 - A.tot %*% xTx.inv.x ))
    ( rss / nrow(nus) ) / ( tr.H/nrow(nus) )**2
}

# gcv2 <- function(alpha, A, nus, p.deg.r=1, p.deg.nu=1) {
    # R <- get_R.p(alpha[1], alpha[2], alpha[3], names(A), p.deg.r, p.deg.nu)
    # #d <- get_d.tot(nus, R)
    # #A.tot <- as.matrix(get_A.tot(A, nus, R))
    # #svd.res <- svd(A.tot)
    # #Sigma.inv <- diag(length(svd.res$d)) * 1/svd.res$d
    # #xTx.inv.x <- svd.res$v %*% Sigma.inv %*% t(svd.res$u)
    # #res <- as.numeric(xTx.inv.x %*% d)
    # XtX <- t(A) %*% A + t(R) %*% R
    # a <- solve(XtX, t(A) %*% d)
    # #res <- 
    # s <- sum(  )
    # y_hat <- rowSums(sweep(A.tot[1:nrow(nus),], MARGIN=2, res, `*`))
    # rss <- sum((y_hat - d[1:nrow(nus)])**2)
    # tr.H <- sum(diag( 1 - A.tot %*% xTx.inv.x ))
    # ( rss / nrow(nus) ) / ( tr.H/nrow(nus) )**2
# }

#get_inversion_coefs <- function(A, R, d) {
#    as.numeric(A %*% solve(t(A) %*% A + t(R) %*% R, t(A) %*% d))
#}

cross_validate <- function(nus, A.mat, R.mat, n.trials=64, n.folds=100) {
    nus.nums <- cut(seq(1, nrow(A.mat)), breaks=n.folds, labels=F)
    chi2_dist <- sapply(parallelMap(function(trial) {
        nus.folds <- nus.nums[sample(nrow(A.mat))]
        sapply(1:n.folds, function(fold_i) {
            nus.test.indices <- which(nus.folds == fold_i, arr.ind=TRUE)
            
            A.train <- A.mat[-nus.test.indices,]
            nus.train <- nus[-nus.test.indices,]
            
            A.tot.train <- get_A.tot(A.train, nus.train, R.mat)
            d.tot.train <- get_d.tot(nus.train, R.mat)
            
            res <- get_inversion_coefs(A.tot.train, d.tot.train)
            
            A.test <- A.mat[nus.test.indices,]
            nus.test <- nus[nus.test.indices,]
            
            A.tot.test <- get_A.tot(A.test, nus.test, R.mat)
            d.tot.test <- get_d.tot(nus.test, R.mat)
            
            Delta.nu <- rowSums(sweep(A.tot.test, MARGIN=2, res, `*`))
            chi2 <- sum((Delta.nu - d.tot.test)**2)
            chi2
        })
    }, trial=1:n.trials), mean)
}

leave_one_out <- function(nus, A.mat, R.mat) {
    as.numeric(parallelMap(function(nu_i) {
        A.train <- A.mat[-nu_i,]
        nus.train <- nus[-nu_i,]
        
        A.tot.train <- get_A.tot(A.train, nus.train, R.mat)
        d.tot.train <- get_d.tot(nus.train, R.mat)
        
        res <- get_inversion_coefs(A.tot.train, d.tot.train)
        
        A.test <- A.mat[nu_i,]
        nus.test <- nus[nu_i,]
        
        A.tot.test <- get_A.tot(A.test, nus.test, R.mat)
        d.tot.test <- get_d.tot(nus.test, R.mat)
        
        Delta.nu <- rowSums(sweep(A.tot.test, MARGIN=2, res$x, `*`))
        chi2 <- sum((Delta.nu - d.tot.test)**2)
        chi2
    }, nu_i=1:(nrow(nus)-1)))
}

two_fold.cv <- function(nus, A.mat, R.mat) {
    #nus.nums <- cut(seq(1, nrow(A.mat)), breaks=n.folds, labels=F)
    nus.nums <- c(1:(nrow(A.mat)-1))%%2 + 1
    sapply(1:2, function(fold_i) {
            nus.test.indices <- which(nus.nums == fold_i, arr.ind=TRUE)
            
            A.train <- A.mat[-nus.test.indices,]
            nus.train <- nus[-nus.test.indices,]
            
            A.tot.train <- get_A.tot(A.train, nus.train, R.mat)
            d.tot.train <- get_d.tot(nus.train, R.mat)
            
            res <- get_inversion_coefs(A.tot.train, d.tot.train)
            
            A.test <- A.mat[nus.test.indices,]
            nus.test <- nus[nus.test.indices,]
            
            A.tot.test <- get_A.tot(A.test, nus.test, R.mat)
            d.tot.test <- get_d.tot(nus.test, R.mat)
            
            Delta.nu <- rowSums(sweep(A.tot.test, MARGIN=2, res$x, `*`))
            sum((Delta.nu - d.tot.test)**2)
            #sum((Delta.nu[1:nrow(nus.test)] - d.tot.test[1:nrow(nus.test)])**2)
    })
}

L_curve <- function() {
    alphas <- c(seq(0.01, 1, 0.01), seq(1.1, 8, 0.1))
    #omegas <- c()
    #thetas <- c()
    #for (alpha in alphas) {
    #alphas.1 <- runif(100, 0, 8)
    #alphas.2 <- runif(100, 0, 8)
    L.curve <- do.call(plyr:::rbind.fill, 
        #parallelMap(function(alpha.1, alpha.2) {
        parallelMap(function(alpha) {
        R.mat <- get_R.p(alpha, alpha, 0, coef.names, r.knots, p.deg.r, p.deg.nu)
        
        A.tot <- get_A.tot(A.mat, nus, R.mat)
        d.tot <- get_d.tot(nus, R.mat)
        res <- get_inversion_coefs(A.tot, d.tot)
        
        Delta.nu <- rowSums(sweep(A.tot, MARGIN=2, res, `*`))
        omega <- sum((Delta.nu[1:nrow(nus)] - d.tot[1:nrow(nus)])**2)
        
        #b.cols <- grepl('b_|c_', coef.names)
        #b.coefs <- res[b.cols]
        #Delta.nu <- rowSums(sweep(A.tot[,b.cols], MARGIN=2, res[b.cols], `*`))
        R.mat2 <- get_R.p(1, 1, 1, coef.names, r.knots, p.deg.r, p.deg.nu)
        Delta.nu2 <- rowSums(sweep(R.mat2, MARGIN=2, res, `*`))
        theta <- sum(Delta.nu2**2)
        #theta <- sum(Delta.nu[-(1:nrow(nus))]**2)
        #theta <- sum(res**2)
        
        data.frame(alpha, omega, theta)
        #omegas <- c(omegas, omega)
        #thetas <- c(thetas, theta)
    }, alpha=alphas))
    #}, alpha.1=alphas.1, alpha.2=alphas.2))
    #psis <- log(omegas)
    #phis <- log(thetas)
    
    
    dpsi <- splinefun(alphas, log(L.curve$omega))
    dphi <- splinefun(alphas, log(L.curve$theta))
    alp <- alphas#seq(min(alphas), max(alphas), length=10000)
    curvature <- (dpsi(alp, 1) * dphi(alp, 2) - dpsi(alp, 2) *
        dphi(alp, 1)) / (dpsi(alp, 1)**2 + dphi(alp, 1)**2)**(3/2)
    #curvature <- (dpsi(alphas, 1) * dphi(alphas, 2) - dpsi(alphas, 2) *
    #    dphi(alphas, 1)) / (dpsi(alphas, 1)**2 + dphi(alphas, 1)**2)**(3/2)
    alpha.L <- alphas[which.max(curvature)]
    
    plot(alp, curvature, type='l')
    #plot(alphas, curvature)
    
    plot(alp, splinefun(alp, sqrt(dpsi(alp, 1)**2 + dphi(alp, 1)**2))(alp, 1))
    
    plot(log(L.curve$omega), log(L.curve$theta))
    dev.off()
    
    plot(L.curve$alpha, L.curve$omega)
    
    Delta.nu <- rowSums(sweep(A.tot, MARGIN=2, res, `*`))
    omega <- sum((Delta.nu[1:nrow(nus)] - d.tot[1:nrow(nus)])**2)
    theta <- sum((Delta.nu[-(1:nrow(nus))] - d.tot[-(1:nrow(nus))])**2)
    
    #A <- as.matrix(A.mat / nus$dnu)
    #R <- as.matrix(R.mat)
    #d <- nus$r.diff / nus$dnu
    #res <- get_inversion_coefs( A.mat, R.mat, nus$r.diff / nus$dnu )
    #res <- as.numeric(A %*% solve(t(A) %*% A + t(R) %*% R, t(A) %*% d))
    
}

optimize_alpha <- function(nus, A.mat, coef.names, r.knots, p.deg.r, p.deg.nu) {
    #cv <- optim(0.1, fn=function(alpha.cv) {
    #cv <- optim(c(1,1,1), fn=function(alpha.cv) {
    cv <- optim(c(2.6, 2.6, 2.6), fn=function(alpha.cv) {
        if (any(alpha.cv < 0)) return(Inf)
        R.mat <- get_R.p(alpha.cv[1], alpha.cv[2], alpha.cv[3], 
            coef.names, r.knots, p.deg.r, p.deg.nu)
        #chi2_dist <- two_fold.cv(nus, A.mat, R.mat)
        #chi2_dist <- cross_validate(nus, A.mat, R.mat)
        #chi2_dist <- leave_one_out(nus, A.mat, R.mat)
        chi2_dist <- two_fold.cv(nus, A.mat, R.mat)
        #chi2_dist <- gcv2(alpha.cv, A.mat, nus, p.deg.r, p.deg.nu)
        chi2 <- mean(chi2_dist)
        #print(paste(num_r_knots, alpha.cv[1], alpha.cv[2], alpha.cv[3], chi2))
        print(paste(alpha.cv[1], alpha.cv[2], alpha.cv[3], chi2))
        chi2
    })#, lower=c(0, 0, 0), upper=c(100, 100, 100), method="Brent")
    #}, lower=0, upper=100, method="Brent")
    #alpha.cv <- cv$par
    alpha.cv <- cv$par
    score <- cv$value
    
    #cv <- optim(15, fn=function(alpha.cv) {
    ##cv <- DEoptim(fn=function(alpha.cv) {
    #    gen.cv <- gcv(alpha.cv, A.mat, nus, p.degree) #1)
    #    print(paste(alpha.cv, gen.cv))
    #    gen.cv
    ##}, lower=0, upper=100, control=list(parallelType=1))
    #}, lower=0, upper=1000, method="Brent")
    ##alpha.cv <- cv$bestmem
    #alpha.cv <- cv$par
    #score <- cv$value
    
    alpha.cv 
}

plot_inversion <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    plot(prof, df_dr,
        type='l', lty=2, col='gray', lwd=3,
        ylim=range(d.m1.f)*2, #c(-0.015, 0.02),
        xlim=c(0, 1),
        axes=F, 
        xlab="Radius r/R",
        ylab="")
    
    lines(r, d.m1.f, lty=1, lwd=2)
    
    abline(h=0, lty=2)
    points(prof, df_dr, col='darkred', pch=20, cex=0.5)
    
    points(r.knots, splinefun(prof, df_dr)(r.knots), col=1, pch=1, cex=0.5)
    
    lines(prof, df_dr+unc, col='darkred', lty=3)
    lines(prof, df_dr-unc, col='darkred', lty=3)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote('Relative difference' ~ delta*.(f.exp)/.(f.exp)))
}

get_profile <- function(xs, coefs, deg, r.knots) {
    sapply(xs, function(x) 
        sum(sapply(1:length(coefs), function(knot_i) 
            coefs[knot_i] * B(x, knot_i, deg, r.knots)) ))
}

plot_results <- function(best) {
    prof <<- seq(0, 1, length=200)
    r.knots <<- best$r.knots
    deg <- best$degree
    
    d.f1_d.r <- get_profile(prof, best$b, best$degree, best$r.knots)
    d.f2_d.r <- get_profile(prof, best$c, best$degree, best$r.knots)
    f1.unc <- get_profile(prof, best$b.e, best$degree, best$r.knots)
    f2.unc <- get_profile(prof, best$c.e, best$degree, best$r.knots)
    
    df_dr <<- d.f1_d.r
    unc <<- f1.unc
    d.m1.f <<- d.m1.f1
    f.exp <<- k.pair$f1.exp
    plot_name <- paste0('inversion-', k.pair$f1, '_', k.pair$f2, 
        '-', model.1$short, '_', model.2$short)
    print(plot_name)
    make_plots(plot_inversion, plot_name)
    
    df_dr <<- d.f2_d.r
    unc <<- f2.unc
    d.m1.f <<- d.m1.f2
    f.exp <<- k.pair$f2.exp
    plot_name <- paste0('inversion-', k.pair$f2, '_', k.pair$f1, 
        '-', model.1$short, '_', model.2$short)
    print(plot_name)
    make_plots(plot_inversion, plot_name)

}

invert <- function(num_r_knots=50, num_nu_knots=0, alpha=0.5, num_r_smooth=500,
        spl.degree=4) {
    r.knots <- get_r.knots.tau(k.m1.f1$x, tau, num_r_knots, spl.degree)
    r.bfs <- get_r.bfs(k.m1.f1$x, r.knots, spl.degree)
    
    nu.knots <- if (num_nu_knots <= 2) { NULL
    } else get_knots(nus$nu.x, num_nu_knots, nu.degree) 
    
    A.mat <- get_A.mat(k.m1.f1$x, r.bfs, nus, k.m1.f1, k.m1.f2, nu.knots, T)
    coef.names <- names(A.mat)
    
    prof <- seq(min(k.m1.f1$x), max(head(k.m1.f1$x, -1)), length=num_r_smooth)
    R.mat <- get_R.mat(alpha, prof, r.knots, coef.names, degree-1)
    
    chi2_dist <- cross_validate(nus, A.mat, R.mat)
    
    A.tot <- get_A.tot(A.mat, nus, R.mat)
    d.tot <- get_d.tot(nus, R.mat)
    res <- get_inversion_coefs(A.tot, d.tot)
    
    best <- list(#chi2=mean(chi2_dist), 
        a=res[grepl('a_', coef.names)],
        b=res[grepl('b_', coef.names)],
        c=res[grepl('c_', coef.names)],
        degree=degree,
        r.knots=r.knots,
        prof=prof)
    
    best
}

invert2 <- function(num_r_knots, num_r_smooth=200, degree=3) {
    r.knots <- get_r.knots(r, tau, num_r_knots, degree)
    r.bfs <- get_r.bfs(r, r.knots, degree)
    
    A.mat <- get_A.mat(r.bfs, nus, k.m1.f1, k.m1.f2)
    coef.names <- names(A.mat)
    
    prof <- seq(min(r), max(head(r, -1)), length=num_r_smooth)
    tried <<- list()
    cv <- optim(par=0.001, fn=function(alpha.cv) {
        alpha.cv <- signif(alpha.cv, 3)
        stringed <- toString(alpha.cv)
        if ( stringed %in% names(tried) ) return(tried[[stringed]])
        R.mat <- get_R.mat(alpha.cv, prof, r.knots, coef.names)
        chi2_dist <- cross_validate(nus, A.mat, R.mat)
        chi2 <- mean(chi2_dist)
        tried[[stringed]] <<- chi2
        print(paste(num_r_knots, alpha.cv, chi2))
        chi2
    }, lower=0, upper=1, method="Brent")
    alpha.cv <- cv$par
    score <- cv$value
    
    R.mat <- get_R.mat(alpha.cv, prof, r.knots, coef.names)
    #chi2_dist <- cross_validate(nus, A.mat, R.mat)
    
    A.tot <- get_A.tot(A.mat, nus, R.mat)
    d.tot <- get_d.tot(nus, R.mat)
    res <- get_inversion_coefs(A.tot, d.tot)
    
    list(chi2=score, 
        a=res[grepl('a_', coef.names)],
        b=res[grepl('b_', coef.names)],
        c=res[grepl('c_', coef.names)],
        r.knots=r.knots,
        prof=prof)
}

p.invert <- function(num_r_knots=80, num_nu_knots=65, 
        spl.degree=4, nu.degree=4, p.deg.r=2, p.deg.nu=2) {
    r.knots <- get_r.knots.tau(k.m1.f1$x, tau, num_r_knots, spl.degree)
    r.bfs <- get_r.bfs(k.m1.f1$x, r.knots, spl.degree)
    
    nu.knots <- if (num_nu_knots <= 2) { NULL
    } else get_knots(nus$nu.x, num_nu_knots, nu.degree) 
    
    A.fname <- paste0("A_mat_", num_r_knots, "_", spl.degree) 
    
    if (file.exists(A.fname)) {
        load(A.fname)
    } else {
        #A.mat <- get_A.mat(k.m1.f1$x, r.bfs, nus, k.m1.f1, k.m1.f2)
        A.mat <- get_A.mat(k.m1.f1$x, r.bfs, nus, k.m1.f1, k.m1.f2, nu.knots, F)
        save(A.mat, file=A.fname)
    }
    
    coef.names <- names(A.mat)
    
    alpha.cv <- optimize_alpha(nus, A.mat, coef.names, r.knots, 
        p.deg.r, p.deg.nu)
    #R.mat <- get_R.p(alpha.cv, coef.names, p.degree)
    R.mat <- get_R.p(alpha.cv[1], alpha.cv[2], alpha.cv[3], coef.names, r.knots,
        p.deg.r, p.deg.nu)
    
    A.tot <- get_A.tot(A.mat, nus, R.mat)
    d.tot <- get_d.tot(nus, R.mat)
    res <- get_inversion_coefs(A.tot, d.tot)
    
    #R.mat <- get_R.p(alpha.cv[1], alpha.cv[2], alpha.cv[3], coef.names, 
    #    p.deg.r, p.deg.nu)
    #low.l <- nus$l<4
    #A.tot <- get_A.tot(A.mat[low.l,], nus[low.l,], R.mat)
    #d.tot <- get_d.tot(nus[low.l,], R.mat)
    #res <- get_inversion_coefs(A.tot, d.tot)
    
    #cv <- optim(c(1,1,1), fn=function(alpha.cv) {
    #    if (any(alpha.cv < 0)) return(Inf)
    #    R.mat <- get_R.p(alpha.cv[1], alpha.cv[2], alpha.cv[3], coef.names,
    #        p.deg.r, p.deg.nu)
    #    #chi2_dist <- two_fold.cv(nus, A.mat, R.mat)
    #    chi2_dist <- cross_validate(nus, A.mat, R.mat)
    #    #chi2_dist <- gcv2(alpha.cv, A.mat, nus, p.deg.r, p.deg.nu)
    #    chi2 <- mean(chi2_dist)
    #    print(paste(num_r_knots, alpha.cv[1], alpha.cv[2], alpha.cv[3], chi2))
    #    chi2
    #})
    
    best <- list(#chi2=score, 
        degree=spl.degree,
        #a=res$x[grepl('a_', coef.names)],
        b=res$x[grepl('b_', coef.names)],
        b.e=res$e[grepl('b_', coef.names)],
        c=res$x[grepl('c_', coef.names)],
        c.e=res$e[grepl('c_', coef.names)],
        r.knots=r.knots)
    
    plot_results(best)
    
    best
}


### SEARCH
# iteration <<- 0
# ran <<- 0
# best_chisq <<- Inf
# best_param <<- c(45, 0.001)
# param_results <<- list()

# run <- function(params) {
    # params[1] <- as.integer(params[1])
    # print(params)
    # num_r_knots <- params[1]
    # alpha <- params[2]
    
    # iteration <<- iteration + 1 
    # if (alpha < 0 || alpha > 2 || num_r_knots < 20 || num_r_knots > 90) return(Inf) 
    
    # stringed <- toString(params)
    # if ( stringed %in% names(param_results) ) return(param_results[[stringed]])
    
    # ran <<- ran + 1 
    # cat(paste("**** iter:", iteration, "; ran:", ran, "\n")) 
    
    # result <- invert(num_r_knots, alpha)
    # chi_sq <- result$chi2 
    # cat(paste("**** chi_sq:", chi_sq, "\n")) 
    
    # if (chi_sq < best_chisq) { 
        # best_chisq <<- chi_sq 
        # best <<- result
        # cat(paste(c("******", params, "\n")))
        # best_param <<- params
        # cat("****** New record!\n")
    # }
    
    # param_results[[stringed]] <- chi_sq
    
    # chi_sq
# }

# result <- optim(par=best_param, fn=run, 
    # #lower=c(20, 0), upper=c(500, 3),
    # control=list(trace=999, parscale=c(10000, 1)))
# print(result)

#best <- invert(best_param[1], best_param[2])



#iteration <<- 0
#ran <<- 0
#best_chisq <<- Inf
#best_param <<- 50
#param_results <<- list()

#run <- function(params) {
#    params[1] <- as.integer(params[1])
#    print(params)
#    num_r_knots <- params[1]
#    
#    iteration <<- iteration + 1 
#    #if (num_r_knots < 20 || num_r_knots > 90) return(Inf) 
#    
#    stringed <- toString(params)
#    if ( stringed %in% names(param_results) ) return(param_results[[stringed]])
#    
#    ran <<- ran + 1 
#    cat(paste("**** iter:", iteration, "; ran:", ran, "\n")) 
#    
#    result <- invert2(num_r_knots)
#    chi_sq <- result$chi2 
#    cat(paste("**** chi_sq:", chi_sq, "\n")) 
#    
#    if (chi_sq < best_chisq) { 
#        best_chisq <<- chi_sq 
#        best <<- result
#        cat(paste(c("******", params, "\n")))
#        best_param <<- params
#        cat("****** New record!\n")
#    }
#    
#    param_results[[stringed]] <- chi_sq
#    
#    chi_sq
#}
#
#result <- optim(par=best_param, fn=run, lower=20, upper=90, method="Brent",
#    control=list(trace=999))
#print(result)

## recover profiles
best <- p.invert()

#r.bfs <- get_r.bfs(prof, r.knots)
#d.f1_d.r <- do.call(c, parallelMap(function(knot_j) {
#    sum(#sapply(1:length(b.coefs), function(knot_i) 
#        #b.coefs[knot_i] * r.bfs[knot_i, knot_j]))
#        b.coefs * r.bfs[, knot_j])#)
#}, knot_j=1:length(prof)))

#rowSums(sweep(r.bfs, MARGIN=2, b.coefs, `*`))


