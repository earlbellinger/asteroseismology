#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
parallelStartMulticore(16)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

models <- get_model_list()

k.pair <- u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2#c2_rho # 
target.name <- 'CygAwball'#'CygAwball'#'hl.no_d'#'CygAbasu3'#'no_diffusion'#'CygAdiff'#
ref.mod <- 'CygAhighRhighM'#'CygAno_diff'#'CygAdiff'#'diffusion'#'hl.diff'#'CygAwball'#
mode.set <- 'CygA'#'BiSON'#'BiSON.MDI'#
perturb <- F#T
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T, subtract.mean=F)# T) 

k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m1$short, '_m-', mode.set)


rs <- seq(0.05, 0.3, 0.05)

#m1.inversion <- minimize_dist(model=m1, rs=rs, initial_params=c(1000,1000,0.03),
#    num_realizations=1, dM=m1$dM, dR=m1$dR)

m1.inversion <- invert.OLA(model=m1, rs=rs, 
    cross.term=866.94, error.sup=1951.9366, width=0.035, dM=m1$dM, dR=m1$dR, 
    subtract.mean=T, num_realizations=256) 

#dev.off();dev.off();dev.off()

m1.inversion2 <- invert.OLA(model=m1, rs=rs, 
    cross.term=866.94, error.sup=1951.9366, width=0.035, #dM=m1$dM, dR=m1$dR,
    subtract.mean=T, num_realizations=256) 


source('OLA_plots.R')

plot_inversion_all2(model=m1, inversion=m1.inversion, 
                    k.pair=k.pair, k.str=k.str, mode.set=mode.set, 
                    legend.spot=NULL, cross.inset="topright",
                    xlim=c(0, 0.45), ylim=c(-0.1, 0.1), plot_nondim=F) 

plot_inversion_all2(model=m1, inversion=m1.inversion2, 
                    k.pair=k.pair, k.str=paste0(k.str, 2), mode.set=mode.set,
                    legend.spot=NULL, cross.inset="topright",
                    xlim=c(0, 0.45), ylim=c(-0.1, 0.1), plot_nondim=F) 



#m1.inversion2 <- minimize_dist(model=m1, rs=rs, 
#    initial_params=c(866.94, 1951.9366, 0.035),
#    num_realizations=1, d.f1.true=F)






