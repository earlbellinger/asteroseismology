#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('kernels.R')
#source('invert_OLA.R')
library(parallel)
library(parallelMap)
library(Bolstad)

solar.mass <- 1.9892e+33
solar.radius <- 6.9598e+10

### MODELS
paths <- list(diffusion=file.path('models', 'diffusion', 'LOGS_MS'),
              no_diffusion=file.path('models', 'no_diffusion', 'LOGS_MS'),
              hl.diff=file.path('..', 'misc', 'solar_calibration', 
                                'high_ell', 'diffusion'),
              hl.no_d=file.path('..', 'misc', 'solar_calibration',
                                'high_ell', 'no_diffusion'),
              modmix=file.path('models', 'BPB2000'),
              CygAdiff=file.path('models', 'CygAdiff', 'LOGS_MS'),
              CygAno_diff=file.path('models', 'CygAno_diff', 'LOGS_MS'),
              CygAwball=file.path('models', 'wball', 'sample_0327'),
              CygBwball=file.path('models', 'CygBwball'),
              CygAbasu1=file.path('models', 'CygAbasu'),
              CygAbasu2=file.path('models', 'CygAbasu'),
              CygAbasu3=file.path('models', 'CygAbasu'),
              CygAbasu4=file.path('models', 'CygAbasu'))

path <- paths$hl.diff
hl.diff <- list(name='MESA Diffusion', short='hlD',
                kerns.dir=file.path(path),
                freq.path=file.path(path, 'diffusion-freqs.dat'),
                #'diffusion-freqs.dat'),
                freq.col.names=F,
                profile.path=file.path(path, 'diffusion.data'),
                fgong.path=file.path(path, 'diffusion.FGONG.dat'))
                #'diffusion.FGONG.dat'))

path <- paths$hl.no_d
hl.no_d <- list(name='MESA No Diffusion', short='hlnoD',
                kerns.dir=file.path(path),
                freq.path=file.path(path, 'no_diffusion-freqs.dat'),
                #'no_diffusion-freqs.dat'),
                freq.col.names=F,
                profile.path=file.path(path, 'no_diffusion.data'),
                fgong.path=file.path(path, 'no_diffusion.FGONG.dat'))
                #'no_diffusion.FGONG.dat'))

path <- paths$diffusion
diffusion <- list(name='Diffusion', short='D',
                  kerns.dir=file.path(path, 'profile1-freqs'),
                  freq.path=file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1-freqs', 
                                             'profile1.data.FGONG.dat'))

path <- paths$no_diffusion
no_diffusion <- list(name='No Diffusion', short='noD',
                     kerns.dir=file.path(path, 'profile1-freqs'),
                     freq.path=file.path(path, 'profile1-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'profile1.data'),
                     fgong.path=file.path(path, 'profile1-freqs', 
                                                'profile1.data.FGONG.dat'))

path <- paths$CygAdiff
CygAdiff <- list(name='16CygA Diffusion', short='CygAD',
                  kerns.dir=file.path(path, 'profile1-freqs'),
                  freq.path=file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1-freqs', 
                                             'profile1.data.FGONG.dat'))

path <- paths$CygAno_diff
CygAno_diff <- list(name='16CygA No Diffusion', short='CygAnoD',
                     kerns.dir=file.path(path, 'profile1-freqs'),
                     freq.path=file.path(path, 'profile1-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'profile1.data'),
                     fgong.path=file.path(path, 'profile1-freqs', 
                                                'profile1.data.FGONG.dat'))

path <- paths$CygAwball
CygAwball <- list(name='16CygA wball', short='CygAwball',
                     kerns.dir=file.path(path, 'sample_0327-freqs'),
                     freq.path=file.path(path, 'sample_0327-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'sample_0327.data'),
                     fgong.path=file.path(path, 'sample_0327-freqs', 
                                                'sample_0327.data.FGONG.dat'))

path <- paths$CygBwball
CygBwball <- list(name='16CygB wball', short='CygBwball',
                     kerns.dir=file.path(path, 'sample_0298-freqs'),
                     freq.path=file.path(path, 'sample_0298-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'sample_0298.data'),
                     fgong.path=file.path(path, 'sample_0298-freqs', 
                                                'sample_0298.FGONG.dat'))

path <- paths$modmix
modmix <- list(name='Modmix', short='Mm',
               freq.path=file.path(path, 'modmix-freqs.dat'), 
               freq.col.names=c('l', 'n', 'nu', 'E'), 
               fgong.path=file.path(path, 'modmix.dat'))

path <- paths$CygAbasu1 # matched to CygAwball
CygAbasu1 <- list(name='16 Cyg A Basu 1', short='CygAbasu1',
               freq.path=file.path(path, '16Cyg_1.freq'), 
               freq.col.names=FALSE, 
               fgong.path=file.path(path, '16Cyg_1.dat'))

path <- paths$CygAbasu2 # matched to CygAdiff
CygAbasu2 <- list(name='16 Cyg A Basu 2', short='CygAbasu2',
               freq.path=file.path(path, '16Cyg_2.freq'), 
               freq.col.names=FALSE, 
               fgong.path=file.path(path, '16Cyg_2.dat'))

path <- paths$CygAbasu3 # matched to CygAwball 
CygAbasu3 <- list(name='16 Cyg A Basu 3', short='CygAbasu3', 
               freq.path=file.path(path, '16Cyg_3.freq'), 
               freq.col.names=FALSE, 
               fgong.path=file.path(path, '16Cyg_3.dat'))

path <- paths$CygAbasu4 # matched to CygAdiff
CygAbasu4 <- list(name='16 Cyg A Basu 4', short='CygAbasu4',
               freq.path=file.path(path, '16Cyg_4.freq'), 
               freq.col.names=FALSE, 
               fgong.path=file.path(path, '16Cyg_4.dat'))

get_model_list <- function() {
    models <- list('hl.diff'=hl.diff,
        'hl.no_d'=hl.no_d,
        'diffusion'=diffusion,
        'no_diffusion'=no_diffusion,
        'modmix'=modmix,
        'CygAdiff'=CygAdiff,
        'CygAno_diff'=CygAno_diff,
        'CygAwball'=CygAwball,
        'CygBwball'=CygBwball,
        'CygAbasu1'=CygAbasu1,
        'CygAbasu2'=CygAbasu2,
        'CygAbasu3'=CygAbasu3,
        'CygAbasu4'=CygAbasu4)
    
    perturbed.model.names <<- c()
    for (R in c(0.98, 1, 1.02)) {
        for (M in c(0.984, 1, 1.016)) {
            path <- file.path('models', 'perturb', 
                paste0('M=', M, '_R=', R, '_D=0'), 'LOGS_MS')
            
            model.name.R <- if (R==0.98) 'lowR' else if (R==1.02) 'highR'
            model.name.M <- if (M==0.984) 'lowM' else if (M==1.016) 'highM'
            model.name <- if (length(model.name.R) < 1 && 
                              length(model.name.M) < 1)
                'Diffusion' else paste0(model.name.R, model.name.M)
            perturbed.model.names <<- c(perturbed.model.names, model.name)
            
            model <- list(name=model.name, short=model.name,
                          kerns.dir=file.path(path, 'profile1-freqs'),
                          freq.path=file.path(path, 'profile1-freqs.dat'),
                          freq.col.names=c('l', 'n', 'nu', 'E'),
                          fgong.path=file.path(path, 'profile1-freqs',
                                               'profile1.data.FGONG.dat'),
                          profile.path=file.path(path, 'profile1.data'))
            
            models[[model.name]] <- model
        }
    }
    
    perturbed.CygA.names <<- c()
    logRs <- c(0.07918125, 0.08635983, 0.09342169)
    Ms <- c(1.064, 1.08, 1.096)
    for (logR in logRs) {
        for (M in Ms) {
            path <- file.path('..', 'calibration', 'calibrate2', 
                paste0('M=', M, '_logR=', logR, '_age=6900000000'), 'LOGS_MS')
            
            model.name.R <- if (logR==logRs[1]) 'lowR' else 
                if (logR==logRs[3]) 'highR'
            model.name.M <- if (M==Ms[1]) 'lowM' else if (M==Ms[3]) 'highM'
            model.name <- if (length(model.name.R) < 1 && 
                              length(model.name.M) < 1)
                'Diffusion' else paste0(model.name.R, model.name.M)
            model.name <- paste0('CygA', model.name)
            perturbed.CygA.names <<- c(perturbed.CygA.names, model.name)
            
            model <- list(name=model.name, short=model.name,
                          kerns.dir=file.path(path, 'profile1-freqs'),
                          freq.path=file.path(path, 'profile1-freqs.dat'),
                          freq.col.names=c('l', 'n', 'nu', 'E'),
                          fgong.path=file.path(path, 'profile1-freqs',
                                               'profile1.data.FGONG.dat'),
                          profile.path=file.path(path, 'profile1.data'))
            
            models[[model.name]] <- model
        }
    }
    
    models
}

parse_freqs <- function(path, col.names=F) {
    if (length(col.names) > 1) { 
        read.table(path, col.names=col.names, stringsAsFactors=F) 
    } else read.table(path, header=1, stringsAsFactors=F) 
}

get_model <- function(freqs, model.name="diffusion", target.name=NULL, 
                      k.pair=u_Y, x0=NULL, square.Ks=F, 
                      fake.dnu=10**-6, match.nl=T) { 
    #model <- get(model.name) 
    model <- models[[model.name]]
    model$target.name <- target.name
    
    # parse model frequencies 
    nus <- parse_freqs(path=model$freq.path, col.names=model$freq.col.names)
    
    # calculate mode inertias
    ell_0 <- nus[nus$l==0,]
    Q_0 <- splinefun(ell_0$nu, ell_0$E)
    Q_norm <- nus$E / Q_0(nus$nu)
    nus <- cbind(nus, data.frame(Q_norm=Q_norm))
    
    # find relative differences with provided frequencies 
    #for (ell in unique(freqs$l)) {
    nus <- if (match.nl) merge(nus, freqs, by=c('l', 'n')) else {
        nus <- do.call(rbind, Map(function(ell) {
            ell.proxy <- freqs[freqs$l == ell,]
            ell.model <- nus[nus$l == ell,]
            closest <- find_closest(ell.proxy$nu, ell.model$nu)
            ell.nu <- ell.model[closest$y,]
            merge(ell.nu, 
                data.frame(n=ell.nu$n, l=ell.nu$l,
                    nu=ell.proxy$nu[closest$x], 
                    dnu=ell.proxy$dnu[closest$x]),
                by=c('l', 'n'))
        }, ell=unique(freqs$l)))
    }
    
    #nus <- merge(nus, freqs, by=c('l', 'n')) 
    diffs <- nus$nu.x - nus$nu.y 
    r.diff <- diffs / nus$nu.x
    d.r.diff <- (nus$dnu / abs(diffs) * abs(r.diff))
    weights <- 1/d.r.diff / sum(1/d.r.diff)
    w.mean <- weighted.mean(r.diff, weights) 
    w.std <- sqrt(sum( weights**2 * d.r.diff**2 ))
    nus <- cbind(nus, data.frame(
        #r.diff = r.diff, 
        #d.r.diff = d.r.diff
        r.diff = r.diff - w.mean, 
        d.r.diff = sqrt( d.r.diff**2 + w.std**2 )
    ))
    nus <- nus[order(nus$l, nus$n),]
    
    # parse model structure 
    model$fgong <- read.table(model$fgong.path, header=1) 
    model$fgong <- model$fgong[order(model$fgong$x),] 
    r <- model$fgong$x 
    model$r <- r 
    
    model$cs.spl <- splinefun(r, sqrt(model$fgong[['c2']])) 
    model$mass <- max(model$fgong$m)
    model$M <- model$mass / solar.mass # 1.9892e+33
    model$radius <- model$fgong$r[which.max(model$fgong$m)]
    model$R <- model$radius / solar.radius # 6.9598e+10
    model$Teff <- if ('profile.path' %in% names(model)) {
        read.table(model$profile.path, header=1, skip=1, nrow=1)$Teff
    } else 5777
    model$nu_ac <- model$M*model$R**-2*(model$Teff/5777)**(-1/2)*5000
    model$nu_max <- model$M*model$R**-2*(model$Teff/5777)**(-1/2)*3090 
    
    if (!is.null(k.pair)) {
        f1 <- model$fgong[[k.pair$f1]]
        f2 <- model$fgong[[k.pair$f2]]
        model$f1 <- f1
        model$f2 <- f2
        model$f1.spl <- splinefun(r, model$f1)
        model$f2.spl <- splinefun(r, model$f2)
        
        # parse model kernel functions 
        k1.fname <- paste0('E_K_', k.pair$f1, '-', k.pair$f2, '.dat')
        k2.fname <- paste0('E_K_', k.pair$f2, '-', k.pair$f1, '.dat')
        k1 <- read.table(file.path(model$kerns.dir, k1.fname), header=1)
        k2 <- read.table(file.path(model$kerns.dir, k2.fname), header=1)
        
        model$f1.exp <- k.pair$f1.exp
        model$f2.exp <- k.pair$f2.exp
        model$f1.name <- k.pair$f1.name
        model$f2.name <- k.pair$f2.name
    }
    
    if (target.name %in% names(models)) {
        target <- models[[target.name]]
        #get(target.name)
        
        m.2 <- read.table(target$fgong.path, header=1) 
        m.2 <- m.2[order(m.2$x),] 
        
        if (!'c2' %in% names(m.2)) m.2 <- cbind(m.2, data.frame(c2=m.2$c**2))
        if (!'u' %in% names(m.2)) m.2 <- cbind(m.2, data.frame(u=m.2$P/m.2$rho))
        
        if (!is.null(k.pair)) {
            m2.f1.spl <- splinefun(m.2$x, m.2[[k.pair$f1]])
            m2.f1 <- m2.f1.spl(r) 
            m2.f2.spl <- splinefun(m.2$x, m.2[[k.pair$f2]])
            m2.f2 <- m2.f2.spl(r) 
            
            model$m2.f1 <- m2.f1
            model$m2.f2 <- m2.f2
            model$m2.f1.spl <- m2.f1.spl
            model$m2.f2.spl <- m2.f2.spl
            
            # calculate relative structural differences 
            model$d.f1.true <- (f1 - m2.f1) / (if (k.pair$f1 == 'Y') 1 else f1)
            model$d.f2.true <- (f2 - m2.f2) / (if (k.pair$f2 == 'Y') 1 else f2)
            
            model$d.f1.spl <- splinefun(r, model$d.f1.true)
            
            k.diffs <- do.call(c, parallelMap(function(mode) {
                if (!mode %in% names(k1) | !mode %in% names(k2)) 
                    return(NA)
                k.1 <- splinefun(k1$x, k1[[mode]])(model$r)
                k.2 <- splinefun(k2$x, k2[[mode]])(model$r)
                integrand <- k.1*model$d.f1.true + k.2*model$d.f2.true
                sintegral(model$r, integrand)$value
            }, mode=paste0('l.', nus$l, '_', 'n.', nus$n)))
            nus <- cbind(nus, k.diffs) 
        }
    } else {
        cat(paste0("Target is real star ", target.name, "\n"))
    }
    
    #if (F) {
    # implement conservation of mass condition
    if (!is.null(k.pair) & (k.pair$f1 == 'rho' || k.pair$f2 == 'rho')) {
        im.mode <- 'l.-1_n.-1'
        
        k1[[im.mode]] <- if (k.pair$f1 == 'rho') {
            4*pi*splinefun(r, f1)(k1$x)*k1$x**2
        } else 0*1:nrow(k1)
        
        k2[[im.mode]] <- if (k.pair$f2 == 'rho') { 
            4*pi*splinefun(r, f2)(k2$x)*k2$x**2 
        } else 0*1:nrow(k2)
        
        nus <- plyr:::rbind.fill(nus, data.frame(l=-1, n=-1, nu.x=0, 
            dnu=fake.dnu, nu.y=0, E=1, E.x=1, r.diff=0, d.r.diff=fake.dnu, 
            Q_norm=1, m1.Q_norm=1))
    }
    #}
    
    model$modes <- paste0('l.', nus$l, '_', 'n.', nus$n)
    model$nus <- nus
    
    if (!is.null(k.pair)) {
        model$k1 <- k1
        model$k2 <- k2
        
        # get square K matrices for OLA inversions
        if (square.Ks) {
            model$K_ijs <- get_square_Ks(modes=model$modes, K=k1, x0=x0)
            model$C_ijs <- get_square_Ks(modes=model$modes, K=k2, x0=x0)
            model$K.ints <- get_K.ints(modes=model$modes, K=k1)
        }
    }
    
    # calculate lower turning points 
    rs <- seq(0, 1, 0.0001)
    model$r_ts <- sapply(model$modes, function(mode) {
        ell <- as.numeric(strsplit(strsplit(mode, '_')[[1]][1], '\\.')[[1]][2])
        nn <- as.numeric(strsplit(strsplit(mode, '_')[[1]][2], '\\.')[[1]][2])
        nu <- model$nus[with(model$nus, l==ell&n==nn),]$nu.x
        if (ell == 0) return(0)
        rs[which.min((
                model$cs.spl(rs)/(rs) - (2*pi*nu)**2/(ell*(ell+1))
            )**2)]
    })
    
    model
}

