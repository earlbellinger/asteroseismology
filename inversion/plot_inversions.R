#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)
debug <- F

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')

parallelStartMulticore(16)

models <- get_model_list()
k.pair <- u_Y
perturb <- T


# add non-dimensional u

target.name <- 'CygAwball'
ref.mod <- 'CygAno_diff'
mode.set <- 'BiSON'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
#m1.inversion <- invert.OLA(model=m1, rs=seq(0.05, 0.25, 0.05), 
#    cross.term=5*10**-2, error.sup=2*10**1, width=0.045)
m1.inversion <- invert.OLA(model=m1, rs=seq(0.05, 0.3, 0.05), 
    cross.term=1000, error.sup=1000, width=0.025)
k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m1$short, '_m-', mode.set)

if (debug) plot_inversion_all(model=m1, 
    inversion=m1.inversion, 
    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
    legend.spot='topright', ylim=c(-.08, .08), xlim=c(0, 0.5)); dev.off() 

plot_inversion_all2(model=m1, inversion=m1.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot='topright', 
    ylim=c(-.05, .05), xlim=c(0, 0.5)) 




target.name <- 'CygAwball'
ref.mod <- 'CygAno_diff'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m2 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
m2.inversion <- invert.OLA(model=m2, rs=seq(0.05, 0.25, 0.05), 
    cross.term=1000, error.sup=100, width=0.04)

k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m2$short, '_m-', mode.set)

if (debug) plot_inversion_all(model=m2, inversion=m2.inversion,
    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
    legend.spot='topright', #ylim=c(-.08, .08), 
    xlim=c(0, 1)); dev.off()

plot_inversion_all2(model=m2, inversion=m2.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot='topright', 
    #ylim=c(-.08, .08), 
    xlim=c(0, 0.5))




target.name <- 'CygA'
ref.mod <- 'CygAwball'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=F) 
m3 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m3$short, '_m-', mode.set)
m3.inversion <- invert.OLA(model=m3, rs=seq(0.1, 0.3, 0.05), 
    cross.term=10**4, error.sup=10**3, width=0.025)

if (debug) plot_inversion_all(model=m3, inversion=m3.inversion,
    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
    legend.spot='topright', ylim=c(-.08, .08), xlim=c(0, 0.4)); dev.off()

plot_inversion_all2(model=m3, inversion=m3.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot='topright', 
    ylim=c(-.08, .08), xlim=c(0, 0.4))

#plot_inversion_all(model=m3, inversion=invert.OLA(model=m3, 
#    rs=seq(0.15, 0.25, 0.05), 
#    cross.term=10**4, error.sup=5*10**2, width=0.01),
#    #, targ.kern.type='mod_sinc'), 
#    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
#    legend.spot='topright', ylim=c(-.08, .08), xlim=c(0, 0.4))



target.name <- 'CygB'
ref.mod <- 'CygBwball'
mode.set <- 'CygB'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=F) 
m4 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m4$short, '_m-', mode.set)
m4.inversion <- invert.OLA(model=m4, rs=seq(0.1, 0.3, 0.05), 
    cross.term=10**4, error.sup=5*10**2, width=0.01)

if (debug) plot_inversion_all(model=m4, inversion=m4.inversion,
    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
    legend.spot='topright', ylim=c(-.08, .08), xlim=c(0, 0.4)); dev.off()

plot_inversion_all2(model=m4, inversion=m4.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot='topright', 
    ylim=c(-.08, .08), xlim=c(0, 0.4))


target.name <- 'modmix'
ref.mod <- 'highRlowM'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m5 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
rs <- seq(0.08, 0.3, 0.02)
m5.inversion <- invert.OLA(model=m5, rs=rs, #seq(0.05, 0.3, 0.01), 
    cross.term=10**3, error.sup=10**2, width=0.0225)

k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m5$short, '_m-', mode.set)

sampler <- c(F, F, T, F, F, T, F, T, F, T, F, F) #

if (debug) plot_inversion_all(model=m5, 
    inversion=m5.inversion, 
    k.pair=k.pair, k.str=k.str, mode.set=mode.set, 
    xlim=c(0, 0.45),
    ylim=c(-0.15, 0.15),
    legend.spot=NULL, #ylim=c(-0.04, 0.1), 
    sampler=sampler)
#, ylim=c(-.2, .2), xlim=c(0, 0.4)); 
    dev.off() 

plot_inversion_all2(model=m5, inversion=m5.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot=NULL, sampler=sampler,
    ylim=c(-.15, .15), xlim=c(0, 0.45))





target.name <- 'CygAwball'
ref.mod <- 'CygAno_diff'
mode.set <- 'CygA'
error.set <- 'BiSON'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=T,
    error.set=error.set) 
m6 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=T) 
m6.inversion <- invert.OLA(model=m6, rs=seq(0.05, 0.25, 0.025), 
    cross.term=1*10**1, error.sup=2*10**2, width=0.035)

k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m6$short, '_m-', mode.set)

if (debug) plot_inversion_all(model=m6, inversion=m6.inversion,
    k.pair=k.pair, k.str=k.str, mode.set=mode.set,
    legend.spot='topright', #ylim=c(-.08, .08), 
    xlim=c(0, 1)); dev.off()

plot_inversion_all2(model=m6, inversion=m6.inversion, k.pair=k.pair, 
    k.str=k.str, mode.set=mode.set, legend.spot='topright', 
    #ylim=c(-.08, .08), 
    xlim=c(0, 0.5))



