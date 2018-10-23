#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

#source('kernels.R')
#source('invert_OLA.R')
source(file.path('..', 'scripts', 'utils.R'))
library(parallel)
library(parallelMap)
library(Bolstad)

#solar.mass <- 1.9892e+33
#solar.radius <- 6.9598e+10

### MODELS
paths <- list(diffusion=file.path('models', 'diffusion', 'LOGS_MS'),
              no_diffusion=file.path('models', 'no_diffusion', 'LOGS_MS'),
              diffusion_final=file.path('models', 'diffusion_final', 'LOGS_3MS'),
              no_diffusion_final=file.path('models', 'no_diffusion_final', 'LOGS_3MS'),
              hl.diff=file.path('..', 'misc', 'solar_calibration', 
                                'high_ell', 'diffusion'),
              hl.no_d=file.path('..', 'misc', 'solar_calibration',
                                'high_ell', 'no_diffusion'),
              modmix=file.path('models', 'BPB2000'),
              ModelS=file.path('models', 'ModelS'),
              CygAdiff=file.path('models', 'CygAdiff', 'LOGS_MS'),
              CygAno_diff=file.path('models', 'CygAno_diff', 'LOGS_MS'),
              CygAwball=file.path('models', 'wball', 'sample_0327'),
              CygBwball=file.path('models', 'CygBwball'),
              CygAbasu1=file.path('models', 'CygAbasu'),
              CygAbasu2=file.path('models', 'CygAbasu'),
              CygAbasu3=file.path('models', 'CygAbasu'),
              CygAbasu4=file.path('models', 'CygAbasu'),
              CygAjcd=file.path('models', 'CygAjcd'),
              CygAyoung=file.path('..', 'calibration', 'old', 'calibrate2',
                'M=1.08_logR=0.08635983_age=6000000000', 'LOGS_MS'),
              CygAyounger=file.path('..', 'calibration', 'old', 'calibrate2',
                'M=1.08_logR=0.08635983_age=5000000000', 'LOGS_MS'),
              calibrate=file.path('models', 'calibrate', 'LOGS_MS'),
              overshoot=file.path('models', 'overshoot', 'os', 'LOGS_3MS'),
              no_overshoot=file.path('models', 'overshoot', 'nos', 'LOGS_3MS'))

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

path <- paths$diffusion_final
diffusion_final <- list(name='Diffusion', short='D',
                  kerns.dir=file.path(path, 'profile1-freqs'),
                  freq.path=file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1-freqs', 
                                             'profile1.data.FGONG.dat'))

path <- paths$no_diffusion_final
no_diffusion_final <- list(name='No Diffusion', short='noD',
                     kerns.dir=file.path(path, 'profile1-freqs'),
                     freq.path=file.path(path, 'profile1-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'profile1.data'),
                     fgong.path=file.path(path, 'profile1-freqs', 
                                                'profile1.data.FGONG.dat'))


path <- paths$ModelS
ModelS <- list(name='Model S', short='ModelS',
                     kerns.dir=file.path(path, 'ModelS-freqs'),
                     freq.path=file.path(path, 'ModelS-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'ModelS.data'),
                     fgong.path=file.path(path, 'ModelS-freqs', 
                                                'ModelS.data.FGONG.dat'))

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

path <- paths$CygAjcd
CygAjcd <- list(name='16CygA JCD', short='CygAJCD',
                     kerns.dir=file.path(path, 'CygA_JCD'),
                     freq.path=file.path(path, 'CygA_JCD-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     #profile.path=file.path(path, 'CygA_JCD.data'),
                     fgong.path=file.path(path, 'CygA_JCD', 
                                                'CygA_JCD.FGONG.dat'))

path <- paths$modmix
modmix <- list(name='Modmix', short='Mm',
               freq.path=file.path(path, 'modmix-freqs.dat'), 
               freq.col.names=c('l', 'n', 'nu', 'E'), 
               fgong.path=file.path(path, 'modmix.dat'),
               M=1.989 * 10**30 / solar_mass,
               R=695.98*10**8 / solar_radius)

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

path <- paths$CygAyoung
CygAyoung <- list(name='16CygA low age', short='CygAyoung',
                     kerns.dir=file.path(path, 'profile1-freqs'),
                     freq.path=file.path(path, 'profile1-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'profile1.data'),
                     fgong.path=file.path(path, 'profile1-freqs', 
                                                'profile1.data.FGONG.dat'))

path <- paths$CygAyounger
CygAyounger <- list(name='16CygA lower age', short='CygAyounger',
                     kerns.dir=file.path(path, 'profile1-freqs'),
                     freq.path=file.path(path, 'profile1-freqs.dat'), 
                     freq.col.names=c('l', 'n', 'nu', 'E'),
                     profile.path=file.path(path, 'profile1.data'),
                     fgong.path=file.path(path, 'profile1-freqs', 
                                                'profile1.data.FGONG.dat'))

path <- paths$calibrate
calibrate <- list(name='Diffusion', short='D',
                  kerns.dir=file.path(path, 'profile1-freqs'),
                  freq.path=file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1-freqs', 
                                             'profile1.data.FGONG.dat'))

path <- paths$overshoot
overshoot <- list(name='Overshoot', short='ov',
                  kerns.dir=NULL,#file.path(path, 'profile1-freqs'),
                  freq.path=NULL,#file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=NA, #c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1.data.FGONG.dat'))

path <- paths$no_overshoot
no_overshoot <- list(name='No overshoot', short='nov',
                  kerns.dir=NULL,#file.path(path, 'profile1-freqs'),
                  freq.path=NULL,#file.path(path, 'profile1-freqs.dat'), 
                  freq.col.names=NA, #c('l', 'n', 'nu', 'E'),
                  profile.path=file.path(path, 'profile1.data'),
                  fgong.path=file.path(path, 'profile1.data.FGONG.dat'))

names.list <<- list()
get_model_list <- function() {
    models <- list('hl.diff'=hl.diff,
        'hl.no_d'=hl.no_d,
        'diffusion'=diffusion,
        'no_diffusion'=no_diffusion,
        'diffusion_final'=diffusion_final,
        'no_diffusion_final'=no_diffusion_final,
        'modmix'=modmix,
        'ModelS'=ModelS, 
        'CygAdiff'=CygAdiff,
        'CygAno_diff'=CygAno_diff,
        'CygAwball'=CygAwball,
        'CygBwball'=CygBwball,
        'CygAbasu1'=CygAbasu1,
        'CygAbasu2'=CygAbasu2,
        'CygAbasu3'=CygAbasu3,
        'CygAbasu4'=CygAbasu4,
        'CygAjcd'=CygAjcd,
        'CygAyoung'=CygAyoung,
        'CygAyounger'=CygAyounger,
        'calibrate'=calibrate,
        'overshoot'=overshoot,
        'no_overshoot'=no_overshoot)
    
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
    names.list[['Sun']] <<- perturbed.model.names
    
    perturbed.CygA.names <<- c()
    logRs <- c(0.07918125, 0.08635983, 0.09342169)
    Ms <- c(1.064, 1.08, 1.096)
    for (logR in logRs) {
        for (M in Ms) {
            path <- file.path('..', 'calibration', 'old', 'calibrate2', 
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
    names.list[['CygA']] <<- perturbed.CygA.names
    
    perturbed.CygB.names <<- c()
    logRs <- c(0.04139269, 0.04921802, 0.05690485)
    Ms <- c(1.015, 1.03, 1.045)
    for (logR in logRs) {
        for (M in Ms) {
            path <- file.path('..', 'calibration', 'old', 'calibrate2', 
                paste0('M=', M, '_logR=', logR, '_age=6800000000'), 'LOGS_MS')
            
            model.name.R <- if (logR==logRs[1]) 'lowR' else 
                if (logR==logRs[3]) 'highR'
            model.name.M <- if (M==Ms[1]) 'lowM' else if (M==Ms[3]) 'highM'
            model.name <- if (length(model.name.R) < 1 && 
                              length(model.name.M) < 1)
                'Diffusion' else paste0(model.name.R, model.name.M)
            model.name <- paste0('CygB', model.name)
            perturbed.CygB.names <<- c(perturbed.CygB.names, model.name)
            
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
    names.list[['CygB']] <<- perturbed.CygB.names
    
    add_models <- function(star.name, models) {
        #name.list <- paste0('perturbed.', star.name, '.names')
        model.names <- c()
        lmh <- c('low', 'mean', 'high')
        for (R in lmh) {
            for (M in lmh) {
                model.type <- paste0(R, 'R', M, 'M')
                model.name <- paste0(star.name, '_', model.type)
                model.names <- c(model.names, model.name)
                
                path <- file.path('models', star.name, model.type, 'LOGS_MS')
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
        
        names.list[[paste0(star.name)]] <<- model.names
        models 
    }
    
    models <- add_models('8006161', models)
    models <- add_models('6106415', models)
    models <- add_models('12258514', models)
    models <- add_models('6225718', models)
    models <- add_models('10068307', models)
    models <- add_models('6116048', models)
    models <- add_models('6106415', models)
    models <- add_models('3632418', models)
    models <- add_models('7510397', models)
    models <- add_models('8938364', models)
    models <- add_models('5774694', models)
    models <- add_models('8760414', models)
    models <- add_models('7970740', models)
    models <- add_models('5184732', models)
    models <- add_models('7940546', models)
    models <- add_models('10162436', models)
    models <- add_models('8379927', models)
    models <- add_models('8228742', models)
    models <- add_models('8694723', models)
    #models <- add_models('6933899', models)
    #models <- add_models('3656476', models)
    models <- add_models('9414417', models)
    models <- add_models('4914923', models)
    models <- add_models('12317678', models)
    models <- add_models('10454113', models)
    models <- add_models('9139163', models)
    models <- add_models('7680114', models)
    models <- add_models('10516096', models)
    models <- add_models('10963065', models)
    models <- add_models('12009504', models)
    models <- add_models('9098294', models)
    models <- add_models('5773345', models)
    models <- add_models('8394589', models)
    models <- add_models('6679371', models)
    models <- add_models('7103006', models)
    models <- add_models('12069449', models)
    models <- add_models('12069424', models)
    
    models
}

parse_freqs <- function(path, col.names=F) {
    if (length(col.names) > 1) { 
        read.table(path, col.names=col.names, stringsAsFactors=F) 
    } else read.table(path, header=1, stringsAsFactors=F) 
}

get_model <- function(model.name, freqs=NULL, target.name=NULL, 
                      k.pair=NULL, x0=NULL, square.Ks=F, 
                      fake.dnu=10**-6, match.nl=T, subtract.mean=F,
                      trim.ks=T) { 
    #model <- get(model.name) 
    if (exists('models') && model.name %in% names(models)) {
        model <- models[[model.name]]
    } else {
        cat(paste("Error: can't find model", model.name, '\n'))
        return(NULL)
    }
    model$target.name <- target.name
    
    # parse model frequencies 
    if (!is.null(model$freq.path) && !is.null(freqs)) {
        nus <- parse_freqs(path=model$freq.path, col.names=model$freq.col.names)
        nus <- nus[as.integer(rownames(unique(nus[,c(1,2)]))),]
        
        # calculate mode inertias
        ell_0 <- nus[nus$l==0,]
        Q_0 <- splinefun(ell_0$nu, ell_0$E)
        Q_norm <- nus$E / Q_0(nus$nu)
        nus <- cbind(nus, data.frame(Q_norm=Q_norm))
        
        # find relative differences with provided frequencies 
        #for (ell in unique(freqs$l)) {
        if (!is.null(freqs)) {
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
            
            missing.nus <- which(! freqs$nu %in% nus$nu.y)
            for (ii in missing.nus) 
                cat(paste0("missing mode n=", freqs[ii,]$n, 
                                      ", l=", freqs[ii,]$l, '\n'))
            
            #nus <- merge(nus, freqs, by=c('l', 'n')) 
            diffs <- nus$nu.x - nus$nu.y 
            r.diff <- diffs / nus$nu.x 
            #d.r.diff <- (nus$dnu / abs(diffs) * abs(r.diff))
            d.r.diff <- nus$dnu / abs(nus$nu.y) * abs(r.diff) 
            weights <- 1/d.r.diff / sum(1/d.r.diff) 
            w.mean <- weighted.mean(r.diff, weights) 
            w.std <- sqrt(sum( weights**2 * d.r.diff**2 )) 
            nus <- cbind(nus, data.frame(
                r.diff = if (subtract.mean) 
                    r.diff - w.mean else r.diff, 
                d.r.diff = if (subtract.mean) 
                    sqrt( d.r.diff**2 + w.std**2 ) else d.r.diff
                #r.diff = r.diff - w.mean, 
                #d.r.diff = sqrt( d.r.diff**2 + w.std**2 )
            ))
            nus <- nus[order(nus$l, nus$n),]
        }
    }
    
    # parse model structure 
    fgong <- read.table(model$fgong.path, header=1) 
    fgong <- fgong[order(fgong$x),] 
    model$fgong <- fgong
    r <- fgong$x 
    model$r <- r 
    
    #model$c2.spl <- splinefun(r, fgong[['c2']]) 
    model$cs.spl <- splinefun(r, sqrt(fgong[['c2']])) 
    model$mass <- max(fgong$m)
    model$M <- model$mass / solar_mass # 1.9892e+33
    model$radius <- fgong$r[which.max(fgong$m)]
    model$R <- model$radius / solar_radius # 6.9598e+10
    model$Teff <- if ('profile.path' %in% names(model)) {
        read.table(model$profile.path, header=1, skip=1, nrow=1)$Teff
    } else 5777
    model$nu_ac <- model$M*model$R**-2*(model$Teff/5777)**(-1/2)*5000
    model$nu_max <- model$M*model$R**-2*(model$Teff/5777)**(-1/2)*3090 
    
    if (!is.null(k.pair)) {
        model$k.pair <- k.pair 
        
        f1 <- fgong[[k.pair$f1]]
        f2 <- fgong[[k.pair$f2]]
        model$f1 <- f1
        model$f2 <- f2
        model$f1.spl <- splinefun(r, model$f1)
        model$f2.spl <- splinefun(r, model$f2)
        
        # parse model kernel functions 
        if (!is.null(model$kerns.dir)) {
            k1.fname <- paste0('E_K_', k.pair$f1, '-', k.pair$f2, '.dat')
            k2.fname <- paste0('E_K_', k.pair$f2, '-', k.pair$f1, '.dat')
            k1 <- read.table(file.path(model$kerns.dir, k1.fname), header=1)
            k2 <- read.table(file.path(model$kerns.dir, k2.fname), header=1)
        }
        
        model$f1.exp <- k.pair$f1.exp
        model$f2.exp <- k.pair$f2.exp
        model$f1.name <- k.pair$f1.name
        model$f2.name <- k.pair$f2.name
    }
    
    if (!is.null(target.name) && target.name %in% names(models)) {
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
            
            if ('m' %in% names(m.2)) {
                model$m2.mass <- max(m.2$m)
                model$m2.M <- model$m2.mass / solar_mass # 1.9892e+33
            } else model$m2.M <- target$M
            if ('r' %in% names(m.2)) {
                model$m2.radius <- m.2$r[which.max(m.2$m)]
                model$m2.R <- model$m2.radius / solar_radius # 6.9598e+10
            } else model$m2.R <- target$R
            
            model$dR <- (model$R-model$m2.R)#/model$R
            model$dM <- (model$M-model$m2.M)#/model$M
            
            model$m2.f1 <- m2.f1
            model$m2.f2 <- m2.f2
            model$m2.f1.spl <- m2.f1.spl
            model$m2.f2.spl <- m2.f2.spl
            
            # calculate relative structural differences 
            model$d.f1.true <- (f1 - m2.f1) / (if (k.pair$f1 == 'Y') 1 else f1)
            model$d.f2.true <- (f2 - m2.f2) / (if (k.pair$f2 == 'Y') 1 else f2)
            
            if (k.pair$f1 == 'u') {
                model$d.f1.nondim <- model$d.f1.true + 
                    model$dR/model$R - model$dM/model$M
            }
            if (k.pair$f2 == 'u') {
                model$d.f2.nondim <- model$d.f2.true + 
                    model$dR/model$R - model$dM/model$M
            }
            
            model$d.f1.spl <- splinefun(r, model$d.f1.true)
            model$d.f2.spl <- splinefun(r, model$d.f2.true)
            
            if (!is.null(model$kerns.dir) && !is.null(model$freq.path)
                    && !is.null(freqs)) {
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
        }
        model$m.2 <- m.2
    } else if (!is.null(target.name)) {
        cat(paste0("Target is real star ", target.name, "\n"))
    }
    
    #if (F) {
    # implement conservation of mass condition
    if (!is.null(k.pair) && (k.pair$f1 == 'rho' || k.pair$f2 == 'rho') &&
            !is.null(model$kerns.dir) && !is.null(model$freq.path) &&
            !is.null(freqs)) {
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
    
    if (!is.null(model$freq.path) && !is.null(freqs)) {
        model$modes <- paste0('l.', nus$l, '_', 'n.', nus$n)
        model$nus <- nus
    }
    
    if (!is.null(k.pair) && !is.null(model$kerns.dir)) {
        include.k1 <- if (trim.ks) names(k1) %in% model$modes else T
        include.k2 <- if (trim.ks) names(k2) %in% model$modes else T
        model$k1 <- cbind(data.frame(x=k1$x), k1[,include.k1])
        model$k2 <- cbind(data.frame(x=k2$x), k2[,include.k2])
        
        # get square K matrices for OLA inversions
        if (square.Ks) {
            model$K.ints <- get_K.ints(modes=model$modes, K=k1)
            model$K_ijs <- get_square_Ks(modes=model$modes, K=k1, x0=NULL)
            model$C_ijs <- get_square_Ks(modes=model$modes, K=k2, x0=NULL)
            if (!is.null(x0)) {
                model$MOLA.K_ijs <- parallelMap(function(r)
                        get_square_Ks(modes=model$modes, K=k1, x0=r), r=x0)
                model$MOLA.C_ijs <- parallelMap(function(r)
                        get_square_Ks(modes=model$modes, K=k2, x0=r), r=x0)
                names(model$MOLA.K_ijs) <- x0
                names(model$MOLA.C_ijs) <- x0
            }
        }
    }
    
    # calculate lower turning points 
    #rs <- seq(0.0001, 1, 0.0001)
    if (!is.null(model$freq.path) && !is.null(freqs)) {
        rs <- model$r * model$R * solar_radius
        model$r_ts <- sapply(model$modes, function(mode) {
            ell <- as.numeric(strsplit(strsplit(mode, '_')[[1]][1], '\\.')[[1]][2])
            nn <- as.numeric(strsplit(strsplit(mode, '_')[[1]][2], '\\.')[[1]][2])
            nu <- model$nus[with(model$nus, l==ell&n==nn),]$nu.x
            if (ell == 0) return(0)
            model$r[which.min((
                    model$cs.spl(rs / (model$R * solar_radius))**2 / (rs**2) 
                    - (10**-6*nu*(2*pi))**2/(ell*(ell+1))
                )**2)]
        })
    }
    
    model
}

