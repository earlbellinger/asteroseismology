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
source('OLA_plots.R')
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

### CONSTANTS 
k.pair  = u_Y
rs      = seq(0.05, 0.3, 0.05) 
sampler = c(T)
models  = get_model_list() 
kern.interp.xs = seq(0, 1, 0.001)

targ.kern.type <- 'mod_Gauss'

for (star in c('CygA', 'CygB')) {

    mode.set     <- star
    error.set    <- star
    target.name  <- star
    ref.mod.name <- paste0(star, 'wball')
    
    targ.mode <- paste0(
        '-p_', target.name, 
        '-m_', mode.set, 
        '-e_', error.set, 
        '-r_', ref.mod.name, 
        paste0("-", targ.kern.type))
    
    load(file.path('save', paste0('inv.lists',       targ.mode)))
    load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
    load(file.path('save', paste0('cross.lists',     targ.mode)))
    
    ref.mod.name <- paste0(star, 'Basu')
    model.prof <- read.table(file.path('models', 'CygniBasu', 
            paste0(star, '.dat')), 
        header=1)
    model <- list(k.pair=k.pair, 
        f1.spl=with(model.prof, splinefun(r, P/rho)))
    
    inversion <- lists_to_inversion(model=model, rs=rs, 
        inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
        cross.lists=cross.lists, inv.params=NULL,
        kern.interp.xs=kern.interp.xs)
    
    make_plots_inversion_all(model, inversion, kern.interp.xs=kern.interp.xs,
        k.str=paste0("-Basu", star), cross.inset="bottomright",
        cross_kern_ylim=c(-0.8, 0.3))

}

