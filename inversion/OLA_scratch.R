#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')
parallelStartMulticore(8)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

models <- get_model_list()

k.pair <- u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2#c2_rho # 
target.name <- 'CygAwball'#'modmix'#'hl.no_d'#'CygAbasu3'#'no_diffusion'#'CygAdiff'#
ref.mod <- 'CygAdiff'#'diffusion'#'hl.diff'#'CygAwball'#'lowRhighM'#
mode.set <- 'CygA'#'BiSON'#'BiSON.MDI'#
perturb <- F
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=perturb) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T) 
#m2 <- get_model(freqs=freqs, model.name='lowRlowM', target.name=target.name, 
#                k.pair=k.pair, square.Ks=T) 

k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
    '-p_', target.name, '_r-', m1$short, '_m-', mode.set)

#save(m1, file='m1.hl.diff')


#plot_kernel_diffs(m1, k.pair, legend.spot='topleft')
make_plots(plot_one_surfless, paste0('kernel-comp', k.str),
    model=m1, k.pair=k.pair, legend.spot='right')
make_plots(plot_kernel_diffs, paste0('kernel-diffs', k.str),
    model=m1, k.pair=k.pair, legend.spot='left')
make_plots(plot_kernel_diffs_surf, paste0('kernel-surf', k.str),
    model=m1, k.pair=k.pair, legend.spot='bottomright')

make_plots(plot_inversion, filename=paste0("inversion", k.str),
    model=m1, inversion=m1.inversion, k.pair=k.pair, legend.spot='bottomright')


make_plots(plot_kernels, filename=paste0("kernels", k.str),
           model=m1, inversion=m1.inversion, xlim=c(0, 1.3))
           #kernels=m1.inversion$avg_kerns, 
           #target_radii=m1.inversion$result$rs, xlim=c(0, 1.3))

#rs <- c(0.13, 0.2, 0.24, 0.3)
#rs <- seq(0.1, 0.25, 0.05)
#rs <- seq(0.11, 0.26, 0.03)
#rs <- seq(0.05, 0.5, 0.06)
#rs <- seq(0, 1, 0.05)
rs <- seq(0.01, 0.35, 0.02)

target.name <- 'modmix'
ref.mod <- 'diffusion'#'CygAdiff'
#k.pair <- rho_Gamma1
#for (target.name in c('CygAbasu1', 'CygAbasu3')) {
for (k.pair in k.pairs) {#list(rho_c2, u_Y)) {
    m1 <- NULL
    inversions <- Map(function(mode.set) {
        freqs <- get_freqs(target.name=target.name, 
            mode.set=mode.set, perturb=perturb) 
        m1 <<- get_model(freqs=freqs, model.name=ref.mod, 
            target.name=target.name, k.pair=k.pair, square.Ks=T) 
        m1.inv <- minimize_dist_individual(model=m1, rs=rs, 
            initial_params=c(1,1))$result
    }, mode.set=c('BiSON'))#, 'CygA'))
    k.str <- paste0('-k_', k.pair$f1, k.pair$f2,
        '-p_', target.name, '_r-', m1$short)
    make_plots(plot_inversions, filename=paste0("inversions", k.str), 
        model=m1, inversions=inversions, k.pair=k.pair, 
        legend.spot='bottomright', xlim=c(0, 0.4))
}
}




m1.inversion <- minimize_dist_individual(model=m1, rs=rs, initial_params=c(1,1,0.05,30))
m1.inversion <- minimize_dist(model=m1, rs=rs, initial_params=c(100,100,0.01,22), use.BG=F)

plot_inversion(model=m1, inversion=m1.inversion, k.pair=k.pair, legend.spot='topleft', ylim=c(-0.1, 0.1)); dev.off()
dev.off()

make_plots(plot_inversion, filename=paste0("inversion", k.str),
    model=m1, inversion=m1.inversion, legend.spot='bottomright')



num_realizations <- 3
inv.params.mat <- parallelMap(function(ii) {
    freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=T) 
    model.list <- list()
    inversion.list <- list()
    for (perturbed.model.name in perturbed.model.names) {
        model <- get_model(freqs=freqs, 
            model.name=perturbed.model.name, 
            target.name=target.name, 
            k.pair=k.pair, square.Ks=T) 
        model.list[[perturbed.model.name]] <- model
        inversion.list[[perturbed.model.name]] <- minimize_dist(model=model, rs=rs)
    }
    sapply(inversion.list, function(inversion) as.numeric(inversion$params))
}, ii=1:num_realizations)

print(inv.params.mat)

inv.params <- Reduce('+', inv.params.mat)/length(inv.params.mat)

print(inv.params)





model.list <- list()
inversion.list <- list()
for (perturbed.model.name in perturbed.model.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=T) 
    model.list[[perturbed.model.name]] <- model
    inversion.list[[perturbed.model.name]] <- minimize_dist(model=model, rs=rs) 
        #invert.OLA(model=model, 
        #rs=rs, cross.term=0, error.sup=1)
}

targ.modes <- paste0('-', target.name, '-', mode.set)
make_plots(plot_freq_diffs, filename="freq-diffs-modmix", model.list=model.list)
make_plots(plot_dimless_freq_diffs, 
    filename=paste0("freq-diffs-dimless", targ.modes), 
    model.list=model.list)
make_plots(plot_function_differences, 
    filename=paste0("funct-diffs", targ.modes), 
    model.list=model.list, inversion.list=inversion.list)
make_plots(plot_inv_diffs, 
    filename=paste0("inv-diffs", targ.modes), 
    model.list=model.list, inversion.list=inversion.list)
make_plots(plot_inv_diffs_lines, 
    filename=paste0("inv-diffs-lines", targ.modes), 
    model.list=model.list, inversion.list=inversion.list)

make_plots(plot_inversion, 
    filename="inversion-modmix-R098M0984", 
    model=model.list[[1]], inversion=inversion.list[[1]])
make_plots(plot_inversion, 
    filename="inversion-modmix-R102M1016", 
    model=model.list[[9]], inversion=inversion.list[[9]])


m1.inversion <- minimize_dist_individual(model=m1, rs=rs)#, targ.kern.type='mod_sinc') 
m2.inversion <- minimize_dist(model=m2, rs=rs) 


plot_inversion(m1, m1.inversion, k.pair) 
plot_inversion(m2, m2.inversion) 


m1.inversion <- minimize_dist(model=m1, rs=rs, initial_params=c(1,1)) 
m1.inversion <- invert.OLA(model=m1, rs=rs, cross.term=1, error.sup=10)
m2.inversion <- invert.OLA(model=m2, rs=rs, 
                           cross.term=1, error.sup=1)

m1.inversion <- invert.OLA(model=m1, 
                           rs=rs, 
                           cross.term=m1.inversion$params[1], #
                           error.sup=m1.inversion$params[2],#, #
                           width=m1.inversion$params[3][[1]])#

m2.inversion <- invert.OLA(model=m2, 
                           rs=rs, 
                           cross.term=m2.inversion$params[1], #257.7427,#
                           error.sup=18.418,#m1.inversion$params[2], #74.98201,#
                           width=m1.inversion$params[3][[1]])#0.04083236)#

make_plots(plot_inversion, filename="SOLA_inversion-no_diffusion_modmix-BiSON-diffparams",
           model=m1, inversion=m1.inversion)
make_plots(plot_kernels, filename="SOLA_inversion-diffusion_modmix-16CygA-avg",
           model=m1, kernels=m1.inversion$avg_kerns, 
           target_radii=m1.inversion$result$rs)
make_plots(plot_kernels, filename="SOLA_inversion-diffusion_modmix-16CygA-cross",
           model=m1, kernels=m1.inversion$cross_kerns, cross=T,
           target_radii=m1.inversion$result$rs)
#make_plots(plot_fiducial_kernel, filename="SOLA_inversion-diffusion_modmix-16CygA-fkern",
#           model=m1, inversion=m1.inversion, make_png=F, short=F, thin=F)
make_plots(plot_convolution, filename="SOLA_inversion-diffusion_modmix-16CygA-convolution",
           model=m1, inversion=m1.inversion, make_png=F, short=F, thin=F)
make_plots(plot_response, filename="SOLA_inversion-diffusion_modmix-16CygA-response",
           model=m1, inversion=m1.inversion, make_png=F, short=F, thin=F)
make_plots(plot_error_corr, filename="SOLA_inversion-diffusion_modmix-16CygA-err_corr",
           model=m1, inversion=m1.inversion, make_png=F, short=F, thin=F)
make_fiducial_kernels_plots(m1, m1.inversion)
make_sensitivities_plots(m1, m1.inversion)
make_MOLA_sensitivity_plots(m1, m1.inversion)

plot_response(m1, m1.inversion)

plot_inversion(m1, invert.OLA(model=m1, rs=rs, 
                              cross.term=m1.inversion$params[1], 
                              error.sup=m1.inversion$params[2], 
                              width=m1.inversion$params[3][[1]])$result)

#plot_kernels(m1$k1$x, m1.inversion$avg_kerns, k.pair$f1.exp, k.pair$f2.exp)
plot_kernels(m1, m1.inversion$avg_kerns, rs)
plot_kernels(m1, m1.inversion$cross_kerns, rs, cross=T)
#plot_kernels(m1$k2$x, m1.inversion$cross_kerns, k.pair$f2.exp, k.pair$f1.exp)

initial_values <- do.call(c, c(m1.inversion$params, m2.inversion$params))
initial_values <- c(3000, 10, 0.01, 3000, 10, 0.01)
    #c(3000, 1, 0.015, 3000, 1, 0.015)


### VARY WIDTHS
best_params <- optim(log10(initial_values), fn=function(inversion_params) {
    inversion_params <- 10**inversion_params
    
    cat(paste("Trying", c("mu:", "width:", "beta:"), 
        inversion_params,'\n'))
    
    cross.term.m1 <- inversion_params[1]
    error.sup.m1 <- inversion_params[2]
    width.m1 <- inversion_params[3]
    m1.inversion <- invert.OLA(model=m1, cross.term=cross.term.m1, 
                               error.sup=error.sup.m1, width=width.m1, rs=rs)
    
    cross.term.m2 <- inversion_params[4]
    error.sup.m2 <- inversion_params[5]
    width.m2 <- inversion_params[6]
    m2.inversion <- invert.OLA(model=m2, cross.term=cross.term.m2, 
                               error.sup=error.sup.m1, width=width.m2, rs=rs)
    
    sum((m1.inversion$result$f - m2.inversion$result$f)**2)
})#, control=list(trace=999, parscale=c(1000, 0.1, 100)))

m1.inversion <- invert.OLA(model=m1, 
                           rs=rs, 
                           cross.term=10**best_params$par[1], #257.7427,#
                           error.sup=10**best_params$par[2], #74.98201,#
                           width=10**best_params$par[3][[1]])#0.04083236)#

m2.inversion <- invert.OLA(model=m2, 
                           rs=rs, 
                           cross.term=10**best_params$par[4], #257.7427,#
                           error.sup=10**best_params$par[5], #74.98201,#
                           width=10**best_params$par[6][[1]])#0.04083236)#


### SAME WIDTH PARAMETER
best_params <- optim(log10(initial_values[1:5]), fn=function(inversion_params) {
    inversion_params <- 10**inversion_params
    
    cat(paste("Trying", c("mu:", "beta:", "width:", "mu2:", "beta2:"), 
              inversion_params,'\n'))
    
    cross.term.m1 <- inversion_params[1]
    error.sup.m1 <- inversion_params[2]
    width.m1 <- inversion_params[3]
    m1.inversion <- invert.OLA(model=m1, cross.term=cross.term.m1, 
                               error.sup=error.sup.m1, width=width.m1, rs=rs)
    
    cross.term.m2 <- inversion_params[4]
    error.sup.m2 <- inversion_params[5]
    width.m2 <- inversion_params[3]
    m2.inversion <- invert.OLA(model=m2, cross.term=cross.term.m2, 
                               error.sup=error.sup.m1, width=width.m2, rs=rs)
    
    sum(abs(m1.inversion$result$f - m2.inversion$result$f))
})

m1.inversion <- invert.OLA(model=m1, 
                           rs=rs, 
                           cross.term=10**best_params$par[1], #257.7427,#
                           error.sup=10**best_params$par[2], #74.98201,#
                           width=10**best_params$par[3][[1]])#0.04083236)#

m2.inversion <- invert.OLA(model=m2, 
                           rs=rs, 
                           cross.term=10**best_params$par[4], #257.7427,#
                           error.sup=10**best_params$par[5], #74.98201,#
                           width=10**best_params$par[3][[1]])#0.04083236)#



### MOLA
best_params <- optim(log10(c(1,1,1,1)), fn=function(inversion_params) {
    inversion_params <- 10**inversion_params
    
    cat(paste("Trying", c("beta:", "mu:", "beta2:", "mu2:"), 
              inversion_params,'\n'))
    
    cross.term.m1 <- inversion_params[1]
    error.sup.m1 <- inversion_params[2]
    m1.inversion <- invert.OLA(model=m1, rs=rs,
                               cross.term=cross.term.m1, 
                               error.sup=error.sup.m1)
    
    cross.term.m2 <- inversion_params[3]
    error.sup.m2 <- inversion_params[4]
    m2.inversion <- invert.OLA(model=m2, cross.term=cross.term.m2, 
                               error.sup=error.sup.m1, rs=rs)
    #distances <- do.call(c, parallelMap(function(avg_kern_i) {
    #    m1.kern <- m1.inversion$avg_kerns[,avg_kern_i]
    #    m2.kern <- m2.inversion$avg_kerns[,avg_kern_i]
    #    sintegral(model$k1$x, (m1.kern - m2.kern)**2)$value
    #}, avg_kern_i=1:ncol(m1.inversion$avg_kerns)))
    
    #sum(distances)
    sum(abs(m1.inversion$result$f - m2.inversion$result$f)**2)
})

m1.inversion <- invert.OLA(model=m1, 
                           rs=rs, 
                           cross.term=10**best_params$par[1], #257.7427,#
                           error.sup=10**best_params$par[2])#0.04083236)#

m2.inversion <- invert.OLA(model=m2, 
                           rs=rs, 
                           cross.term=10**best_params$par[3], #257.7427,#
                           error.sup=10**best_params$par[4])#0.04083236)#

