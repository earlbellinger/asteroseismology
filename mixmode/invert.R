#### Invert Gemma models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Department of Astronomy, Yale University;
#### Stellar Ages & Galactic Evolution Group,
#### Max-Planck-Institut fur Sonnensystemforschung

library(parallel)
library(parallelMap)
library(data.table)
library(pracma)
library(akima)
library(Bolstad)
library(RColorBrewer)
source(file.path('..', 'inversion', 'kernels.R'))
source(file.path('..', 'inversion', 'invert_OLA.R'))
source(file.path('..', 'inversion', 'models.R'))
source(file.path('..', 'inversion', 'OLA_plots.R'))
#source(file.path('..', 'inversion', 'frequencies.R'))
#source(file.path('/', 'scratch', 'seismo', 'bellinger',
#    'asteroseismology', 'scripts', 'seismology.R'))
source(file.path('/', 'scratch', 'seismo', 'bellinger',
    'asteroseismology', 'scripts', 'utils.R'))

parallelStartMulticore(16)

k.pair <- u_Y #rho_c2 #u_Y #c2_rho #
ells <- c(10, 11, 9, 3)
dnus <- c(0.18, 0.32, 0.33, 0.12, 0.18, 0.16, 0.16, 0.15, 0.16, 0.15, 
          0.26, 0.11, 0.24, 0.17, 0.19, 0.19, 0.11, 0.11, 0.21, 0.25, 0.23, 
          0.19, 0.22, 0.23, 0.15, 0.18, 0.18, 0.12, 0.22, 0.28, 
          0.47, 0.37, 0.44)

get_closest <- function(freqs, nu_max, ells) {
    new.freqs <- data.frame()
    for (ell in 1:length(ells)) {
        freqs.ell <- freqs[freqs$l==ell-1,]
        for (ii in 1:ells[ell]) {
            closest <- find_closest2(freqs.ell$nu, nu_max)
            new.freqs <- plyr:::rbind.fill(new.freqs, freqs.ell[closest,])
            freqs.ell <- freqs.ell[-closest,]
        }
    }
    new.freqs[order(new.freqs$l, new.freqs$n),]
}

get_freqs <- function(model_number, log_dir, freq_files, mdl_nums) {
    idx <- which(mdl_nums == model_number)
    freq_file <- freq_files[idx]
    read.table(file.path(log_dir, freq_file), header=F,
        col.names=c('l', 'n', 'nu', 'E'))#, gyre=T)
}

collapse_freqs <- function(model_number, log_dir, ev.DF, freq_files, mdl_nums,
        closest=T) {
    #print(model_number)
    DF <- get_freqs(model_number, log_dir, freq_files, mdl_nums)
    if (nrow(DF) <= 0) return(data.frame(NA))
    if (closest) {
        nu_max <- ev.DF[ev.DF$model_number==model_number,]$nu_max
        DF <- get_closest(freqs=DF, nu_max=nu_max, ells=ells)
    }
    freqs <- rbind(DF$nu)
    colnames(freqs) <- paste0('l', DF$l, '.', 'n', DF$n)
    as.data.frame(freqs)
}

tracks <- Map(function(dir.name) {
    print(dir.name)
    log_dir <- file.path(dir.name, 'LOGS')
    ev.DF <- read.table(file.path(log_dir, 'history.data'), header=TRUE, skip=5)
    logs <- list.files(log_dir)
    
    freq_files <- logs[grep('profile.+-freqs.dat$', logs)]
    prof_files <- sub('-freqs.dat', '.data', freq_files)
    
    mdl_nums <- as.numeric(parallelMap(function(pro_file)
            read.table(pro_file, header=TRUE, nrows=1, skip=1)$model_number,
        pro_file=file.path(log_dir, prof_files)))
    ages <- as.numeric(parallelMap(function(pro_file)
            read.table(pro_file, header=TRUE, nrows=1, skip=1)$star_age,
        pro_file=file.path(log_dir, prof_files))) / 10**9
    min.age <- min(ages)
    ages <- ages[order(mdl_nums)] - min.age
    mdl_nums. <- sort(mdl_nums)
    
    track_freqs <- do.call(plyr:::rbind.fill, 
        parallelMap(function(model_number) 
            collapse_freqs(model_number, log_dir, ev.DF, freq_files, mdl_nums), 
        model_number=mdl_nums.))
    track_freqs <- track_freqs[,order(names(track_freqs))]
    
    all_freqs <- do.call(plyr:::rbind.fill, 
        parallelMap(function(model_number) 
            collapse_freqs(model_number, log_dir, ev.DF, freq_files, mdl_nums,
                closest=F), 
        model_number=mdl_nums.))
    all_freqs <- all_freqs[,order(names(all_freqs))]
    
    list(ev.DF=ev.DF, 
         mdl.nums.=mdl_nums., 
         ages=ages, 
         track_freqs=track_freqs,
         all_freqs=all_freqs,
         path=log_dir)
}, dir.name=c('gemma', 'gemma0'))

attach(tracks[['gemma']])
mdls <- c(25, 100, 136, 143, 180, 230)
cat(paste("Mass", "Radius", "Teff", "L/L_Sun", "Age\n", sep=' '))
for (model_number in mdls+min(mdl.nums.)) {
    selector <- ev.DF$model_number == model_number
    model <- ev.DF[selector,]
    properties <- c(model$star_mass, 
        10**model$log_R, 
        10**model$log_Teff, 
        10**model$log_L, 
        model$star_age/10**9)
    cat(signif(properties, 4))
    cat('\n')
}

make.df <- function(freqs) {
    data.frame(
        l=as.numeric(sapply(names(freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][1], 'l')[[1]][2])),
        n=as.numeric(sapply(names(freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][2], 'n')[[1]][2])),
        nu=as.numeric(freqs), 
        dnu=dnus) 
}

for (model_number in mdls) {
    freqs <- track_freqs[model_number,]
    freqs. <- freqs[,!is.na(freqs)]
    
    gemma.freqs <- make.df(freqs.)
    
    freqs.0 <- track_freqs[model_number-1,]
    freqs.0 <- freqs.0[,!is.na(freqs.0)]
    
    proxy_mdl <- model_number-1
    proxy_freqs <- freqs.0
    gemma0.freqs <- make.df(proxy_freqs)
    
    path <- tracks[['gemma']]$path
    Gemma0 <- list(name='Gemma0', short='Gemma0', 
        kerns.dir=file.path(path, paste0('profile', proxy_mdl, '-freqs')), 
        freq.path=file.path(path, paste0('profile', proxy_mdl, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'), 
        profile.path=file.path(path, paste0('profile', proxy_mdl, '.data')), 
        fgong.path=file.path(path, paste0('profile', proxy_mdl, '-freqs'), 
            paste0('profile', proxy_mdl, '.data.FGONG.dat')))
    
    
    
    
    if (F) {
    
    
    #gemma0.freqs$nu <- with(gemma0.freqs, rnorm(nrow(gemma0.freqs), nu, dnu))
    
    path <- tracks[['gemma0']]$path
    Gemma0 <- list(name='Gemma0', short='Gemma0', 
        kerns.dir=file.path(path, paste0('profile', proxy_mdl, '-freqs')), 
        freq.path=file.path(path, paste0('profile', proxy_mdl, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'), 
        profile.path=file.path(path, paste0('profile', proxy_mdl, '.data')), 
        fgong.path=file.path(path, paste0('profile', proxy_mdl, '-freqs'), 
            paste0('profile', proxy_mdl, '.data.FGONG.dat')))
    
    
    gemma0.freqs <- data.frame(
        l=as.numeric(sapply(names(proxy_freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][1], 'l')[[1]][2])),
        n=as.numeric(sapply(names(proxy_freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][2], 'n')[[1]][2])),
        nu=as.numeric(proxy_freqs),
        dnu=dnus)
    #gemma0.freqs$nu <- with(gemma0.freqs, rnorm(nrow(gemma0.freqs), nu, dnu))
    
    path <- tracks[['gemma0']]$path
    Gemma0 <- list(name='Gemma0', short='Gemma0', 
        kerns.dir=file.path(path, paste0('profile', proxy_mdl, '-freqs')), 
        freq.path=file.path(path, paste0('profile', proxy_mdl, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'), 
        profile.path=file.path(path, paste0('profile', proxy_mdl, '.data')), 
        fgong.path=file.path(path, paste0('profile', proxy_mdl, '-freqs'), 
            paste0('profile', proxy_mdl, '.data.FGONG.dat')))
    
    freqs.0 <- tracks$gemma0$all_freqs
    freqs.0 <- freqs.0[,names(freqs.)]
    freqs.0 <- freqs.0[complete.cases(freqs.0),]
    
    chi2s <- parallelMap(function(ii) {
        gemma0.freqs <- freqs.0[ii,]
        
        # remove surf term 
        ##nu <- freqs. - gemma0.freqs
        #r.diff <- as.numeric(freqs. - gemma0.freqs)
        #inertia <- m1$nus$Q_norm #m1$nu_ac * #* m1$nus$d.r.diff
        #Xpinv <- ginv( matrix(c(freqs.**-1, freqs.**3) / inertia, ncol=2) )
        
        #a.r.1 <- Xpinv %*% ( r.diff )#/ m1$nus$d.r.diff )
        #F_surf <- ( a.r.1[[1]]*freqs.**-1 + a.r.1[[2]]*freqs.**3 ) / inertia
        
        ##K.diffs.surf <- k.diffs + F_surf #r.diff - 
        sum( ((freqs. - gemma0.freqs)/dnus)**2 )
        #sum( ((freqs. + F_surf - gemma0.freqs)/dnus)**2 )
    }, ii=1:nrow(freqs.0))
    
    proxy_freqs <- freqs.0[order(as.numeric(chi2s))[1],]
    #freqs.0[which.min(chi2s),]
    proxy_mdl <- as.numeric(rownames(proxy_freqs))
    
    gemma0.freqs <- data.frame(
        l=as.numeric(sapply(names(proxy_freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][1], 'l')[[1]][2])),
        n=as.numeric(sapply(names(proxy_freqs), 
          function(mode) strsplit(strsplit(mode, '\\.')[[1]][2], 'n')[[1]][2])),
        nu=as.numeric(proxy_freqs),
        dnu=dnus)
    gemma0.freqs$nu <- with(gemma0.freqs, rnorm(nrow(gemma0.freqs), nu, dnu))
    
    path <- tracks[['gemma0']]$path
    Gemma0 <- list(name='Gemma0', short='Gemma0', 
        kerns.dir=file.path(path, paste0('profile', proxy_mdl, '-freqs')), 
        freq.path=file.path(path, paste0('profile', proxy_mdl, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'), 
        profile.path=file.path(path, paste0('profile', proxy_mdl, '.data')), 
        fgong.path=file.path(path, paste0('profile', proxy_mdl, '-freqs'), 
            paste0('profile', proxy_mdl, '.data.FGONG.dat')))
      
    }
    
    
    
    
    path <- tracks[['gemma']]$path
    Gemma <- list(name='Gemma', short='Gemma',
        kerns.dir=file.path(path, paste0('profile', model_number, '-freqs')),
        freq.path=file.path(path, paste0('profile', model_number, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'),
        profile.path=file.path(path, paste0('profile', model_number, '.data')),
        fgong.path=file.path(path, paste0('profile', model_number, '-freqs'), 
            paste0('profile', model_number, '.data.FGONG.dat')))
    
    models <- list(Gemma=Gemma, Gemma0=Gemma0)
    
    m1 <- get_model(freqs=gemma0.freqs, model.name='Gemma', 
        target.name='Gemma0', k.pair=k.pair, square.Ks=T) 
    
    m1 <- get_model(freqs=gemma.freqs, model.name='Gemma0', 
        target.name='Gemma', k.pair=k.pair, square.Ks=T) 
    
    k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
        '_r-', m1$short, model_number)
    
    make_plots(plot_one_surfless, paste0('kernel-comp', k.str),
        model=m1, k.pair=k.pair, legend.spot='bottomright')
    make_plots(plot_kernel_diffs, paste0('kernel-diffs', k.str),
        model=m1, k.pair=k.pair, legend.spot='left')
    make_plots(plot_kernel_diffs_surf, paste0('kernel-diffs-surf', k.str),
        model=m1, k.pair=k.pair, legend.spot='left')
    
    #rs <- c(0.25) #c(0.1, 0.2, 0.3)
    #rs <- c(0.05)#
    #rs <- c(0.001, 0.01, 0.1)
    #rs <- seq(0.001, 0.3, 0.001)
    rs <- c(0.002, seq(0.01, 0.3, 0.01))
    
    #m1.inversion <- invert.OLA(model=m1, rs=rs, cross.term=1, error.sup=0,
    #    width=0.01, use.BG=T)
    #m1.inversion <- invert.OLA(model=m1, rs=rs, cross.term=10**5, error.sup=1)#0.001)
    
    #m1.inversion <- minimize_dist(model=m1, rs=rs, initial_params=c(10, 1))#, 0.01))
    #    targ.kern.type='mod_sinc')
    #plot_inversion(m1, m1.inversion, log='x', xlim=c(min(rs)/2, 1))
    
    rs <- 10**seq(-4, log10(0.4), 0.1)
    
    m1.inversion <- minimize_dist_individual(model=m1, rs=rs, 
        initial_params=c(100, 100, 0.01, 10))#, targ.kern.type='mod_sinc')
    
    m1.inversion <- minimize_dist(model=m1, rs=rs, 
        initial_params=c(100, 100, 0.01))
    
    plot_inversion(model=m1, inversion=m1.inversion, k.pair=k.pair, log='x', ylim=c(-0.25, 0.2)); dev.off()
    
    make_plots(plot_inversion, paste0(model_number),
        model=m1, inversion=m1.inversion, legend.spot='topright')
    
    make_plots(plot_inversion, paste0(model_number, 'log'),
        model=m1, inversion=m1.inversion, legend.spot='topright',
        log='x', xlim=c(min(rs)/1.5, min(1, max(rs+.1))))#1))
    
    #par(mfrow=c(2,1))
    #plot_inversion(m1, m1.inversion, log='x', xlim=c(min(rs)/2, 1))
    #plot_kernels(m1, m1.inversion$avg_kerns, rs, log='x', xlim=c(min(rs)/2, 1))
    #dev.off()
    ##make_plots(plot_one_surfless, paste0(model_number, '-forward'),
    ##    model=m1, k.pair=k.pair)
    sidebyside <- function(model, k.pair, ...,
            text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font) {
        par(mfrow=c(1,2))
        plot_one_surfless(model=model, k.pair=k.pair, text.cex=text.cex, 
            mgp=mgp, mar=mar, font=font, legend.spot='center')
        plot_kernel_diffs(model=model, k.pair=k.pair, text.cex=text.cex, 
            mgp=mgp, mar=mar, font=font, legend.spot='topright')
    }
    make_plots(sidebyside, paste0("sidebyside-", model_number),
        model=m1, k.pair=k.pair)
}






