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
source(file.path('/', 'scratch', 'seismo', 'bellinger',
    'asteroseismology', 'scripts', 'utils.R'))
source(file.path('..', 'inversion', 'kernels.R'))
source(file.path('..', 'inversion', 'OLA_invert.R'))
source(file.path('..', 'inversion', 'models.R'))
source(file.path('..', 'inversion', 'OLA_plots.R'))
#source(file.path('..', 'inversion', 'frequencies.R'))
#source(file.path('/', 'scratch', 'seismo', 'bellinger',
#    'asteroseismology', 'scripts', 'seismology.R'))

parallelStartMulticore(16)

k.pair <- u_Y #rho_c2 #u_Y #c2_rho #
if (F) {
ells <- c(10, 11, 9, 3)
dnus <- c(0.18, 0.32, 0.33, 0.12, 0.18, 0.16, 0.16, 0.15, 0.16, 0.15, 
          0.26, 0.11, 0.24, 0.17, 0.19, 0.19, 0.11, 0.11, 0.21, 0.25, 0.23, 
          0.19, 0.22, 0.23, 0.15, 0.18, 0.18, 0.12, 0.22, 0.28, 
          0.47, 0.37, 0.44)
} else {
ells <- c(15, 15, 15, 11)
dnus <- c(0.0869020395071, 0.0966240629865, 0.0717355776071, 0.0602866261678, 0.0569259456774, 0.0482559959075, 0.0425987693245, 0.0405089288348, 0.0449376026152, 0.0584052965449, 0.100838859713, 0.151057946903, 0.419247010572, 1.04037503944, 1.90332442698, 0.101064009232, 0.117989028239, 0.0705683160588, 0.0565860281984, 0.05624610449, 0.0396667975797, 0.0410012932975, 0.041718087988, 0.0479737201351, 0.0631514715634, 0.0972306091726, 0.154562070684, 0.290389283301, 0.643039480533, 0.767873085246, 0.312873839805, 0.272978980918, 0.159850943565, 0.13704824559, 0.102349870847, 0.0831586146055, 0.0608233749943, 0.0598264219113, 0.0662369972022, 0.0793571905516, 0.123912247239, 0.235129409308, 0.507109031236, 1.53587077979, 2.74260008225, 5.13722504036, 1.4983748475, 1.22370713213, 0.373924294046, 0.299447049239, 0.146571555671, 0.168103089458, 0.276270137214, 0.435047563273, 1.12826860282, 1.41549015796)
}

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

for (mdl_num in mdls) {
    proxy_freqs <- track_freqs[mdl_num,]
    proxy_freqs <- make.df(proxy_freqs[,!is.na(proxy_freqs)])
    ref_mdl_num <- mdl_num + 1 
    
    path <- tracks[['gemma']]$path
    ref_mod <- list(name='Gemma0', short='Gemma0', 
        kerns.dir=file.path(path, paste0('profile', ref_mdl_num, '-freqs')), 
        freq.path=file.path(path, paste0('profile', ref_mdl_num, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'), 
        profile.path=file.path(path, paste0('profile', ref_mdl_num, '.data')), 
        fgong.path=file.path(path, paste0('profile', ref_mdl_num, '-freqs'), 
            paste0('profile', ref_mdl_num, '.data.FGONG.dat')))
    proxy_star <- list(name='Gemma', short='Gemma',
        kerns.dir=file.path(path, paste0('profile', mdl_num, '-freqs')),
        freq.path=file.path(path, paste0('profile', mdl_num, '-freqs.dat')), 
        freq.col.names=c('l', 'n', 'nu', 'E'),
        profile.path=file.path(path, paste0('profile', mdl_num, '.data')),
        fgong.path=file.path(path, paste0('profile', mdl_num, '-freqs'), 
            paste0('profile', mdl_num, '.data.FGONG.dat')))
    
    models <- list(proxy_star=proxy_star, ref_mod=ref_mod)
    
    rs <- 0.01 #c(0.01)
    m1 <- get_model(freqs=proxy_freqs, model.name='ref_mod', 
        target.name='proxy_star', k.pair=k.pair, square.Ks=T, x0=rs) 
    
    #k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    #    '_r-', m1$short, mdl_num)
    #rs <- seq(0.15, 0.25, 0.025)
    
    m1.inversion <- minimize_dist(model=m1, rs=rs, initial_params=c(10e-10, 10))
    
    m1.inversion <- invert.OLA(model=m1, rs=rs, cross.term=0, 
        error.sup=10**3, use.BG=T, dM=m1$dM, dR=m1$dR, 
        subtract.mean=F, num_realizations=1, perturb=F) 
    make_plots_inversion_all(inversion=m1.inversion, model=m1, xlim=c(0, 0.45),
        plot_nondim=F, log='x', avg_kern_ylim=NULL, use.cairo=F, font='Times') 
    
    
    m1.inversion <- invert.OLA(model=m1, rs=rs, cross.term=0, 
        error.sup=1e1, use.BG=T, dM=m1$dM, dR=m1$dR, 
        subtract.mean=F, num_realizations=1, perturb=F)#, BG_pows=-2) 
    make_plots_inversion_all(inversion=m1.inversion, model=m1, xlim=c(0, 1),
        plot_nondim=F, log='x', avg_kern_ylim=NULL, use.cairo=F, font='Times',
        make.inset=F) 
    
    
    #make_plots_inversion_all(inversion=m1.inversion, model=m1, k.str=k.str, 
    #    xlim=c(0, 0.45), k.pair=k.pair, mode.set='Gemma', plot_nondim=F, 
    #    log='x', avg_kern_ylim=NULL)
    
    #m1.inversion <- minimize_dist(model=m1, rs=rs, 
    #    initial_params=c(100, 100, 0.01), dM=m1$dM, dR=m1$dR)
}






