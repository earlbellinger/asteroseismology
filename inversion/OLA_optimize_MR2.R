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

### LISTS 
inv.lists       <- list()
avg.kerns.lists <- list()
cross.lists     <- list()
star.Ms         <- c()
star.Rs         <- c()
cross.terms     <- c()
error.sups      <- c()
widths          <- c()
dir.create('save', showWarnings = FALSE)

### COMMAND LINE ARGUMENTS
args <- commandArgs(TRUE)
mode.set        <- if (length(args)>0)   args[1] else 'CygA'
error.set       <- if (length(args)>1)   args[2] else 'CygA'
target.name     <- if (length(args)>2)   args[3] else 'CygAwball'
ref.mod.name    <- if (length(args)>3)   args[4] else 'diffusion'
model.list.name <- if (length(args)>4)   args[5] else 'perturbed.CygA.names'
targ.kern.type  <- if (length(args)>5)   args[6] else 'mod_Gauss'
initial.M       <- if (length(args)>6)   as.numeric(args[7])  else 1.08
initial.R       <- if (length(args)>7)   as.numeric(args[8])  else 1.22
sigma.M         <- if (length(args)>8)   as.numeric(args[9])  else 0.016
sigma.R         <- if (length(args)>9)   as.numeric(args[10]) else 0.02
n_trials        <- if (length(args)>10)  as.numeric(args[11]) else 128
perturb         <- if (length(args)>11)  as.logical(args[12]) else T

### PREPARE MODELS 
targ.mode <- paste0(
    '-p_', target.name, 
    '-m_', mode.set, 
    '-e_', error.set, 
    '-r_', ref.mod.name, 
    paste0("-", targ.kern.type))

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
    error.set=error.set, perturb=perturb) 

ref.mod <- get_model(freqs=freqs, model.name=ref.mod.name, 
    target.name=target.name, k.pair=k.pair, square.Ks=F)

model.names <- get(model.list.name) 
model.list <- parallelMap(function(model.name) 
        get_model(freqs=freqs, model.name=model.name, 
            target.name=target.name, k.pair=k.pair, square.Ks=T), 
    model.name=model.names)
names(model.list) <- model.names

### RUN TRIALS 
for (trial_i in 1:n_trials) {
    
    print(paste("TRIAL_I =", trial_i))
    
    # perturb quantities 
    nu.star <- rnorm(length(freqs$nu), freqs$nu, freqs$dnu)
    repeat {
        trial.M <- rnorm(1, initial.M, sigma.M)
        trial.R <- rnorm(1, initial.R, sigma.R)
        if (trial.M > initial.M - 10*sigma.M &&
            trial.M < initial.M + 10*sigma.M &&
            trial.R > initial.R - 10*sigma.R && 
            trial.R < initial.R + 10*sigma.R) break
    }
    for (perturbed.model.name in model.names) {
        model <- model.list[[perturbed.model.name]]
        model$nus$nu.y <- nu.star
        model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
            nu_ac=model$nu_ac)
        model.list[[perturbed.model.name]] <- model
    }
    
    # set up optimization 
    initial_params <- c(1000, 1000, 0.02, trial.M, trial.R)
    parscale <- c(0.1, 0.1, 0.01, sigma.M, sigma.R)
    
    iteration <<- 0
    best_result <<- Inf
    best_params <<- initial_params
    
    SOLA_optim <- function(log10_inversion_params) {
        inversion_params <- 10**log10_inversion_params
        
        cross.term <- inversion_params[1]
        error.sup  <- inversion_params[2]
        width      <- inversion_params[3]
        star.M     <- inversion_params[4]
        star.R     <- inversion_params[5]
        if (star.M < initial.M-10*sigma.M || star.M > initial.M+10*sigma.M || 
            star.R < initial.R-10*sigma.R || star.R > initial.R+10*sigma.R ||
            width < 0.001 || width > 0.1)
          return(Inf)
        
        iteration <<- iteration + 1
        cat(paste("**** iter:", iteration, "\n"))
        
        cat('Trying inversion params ')
        cat(inversion_params)
        cat('\n')
        
        inv.list <- list()
        inv.list <- parallelMap(function(model_i) {
            model <- model.list[[model_i]]
            invert.OLA(model=model, rs=rs, 
                cross.term=cross.term, error.sup=error.sup, width=width, 
                targ.kern.type=targ.kern.type, 
                get_cross_kerns=F, get_avg_kerns=F, 
                subtract.mean=F, dM=model$M-star.M, dR=model$R-star.R, 
                perturb=F, num_realizations=1, F_surf=model$F_surf)$result
        }, model_i=1:length(model.list))
        names(inv.list) <- model.names
        
        f.means <- sapply(inv.list, function(inversion) inversion$f) 
        result <- sum(apply(f.means, 1, var)) 
        result <- log(result) - 
            dnorm(star.M, trial.M, sigma.M, log=T) - 
            dnorm(star.R, trial.R, sigma.R, log=T) 
        cat(paste("Result:", format(result, digits=8), '\n')) 
        
        if (result < best_result) {
            best_result <<- result
            best_params <<- inversion_params
            cat(paste(c("*****", inversion_params, "\n")))
            cat("***** New record!\n")
        }
        
        result 
    }
    
    best_params <- optim(log10(initial_params), fn=SOLA_optim, 
        control=list(trace=999, parscale=parscale, maxit=512))
    print(best_params)
    
    inversion_params <- 10**best_params$par
    cross.term <- inversion_params[1]
    error.sup  <- inversion_params[2]
    width      <- inversion_params[3]
    star.M     <- inversion_params[4]
    star.R     <- inversion_params[5]
    
    inv.list <- parallelMap(function(model_i) {
        model <- model.list[[model_i]]
        invert.OLA(model=model, rs=rs, 
            cross.term=cross.term, error.sup=error.sup, width=width, 
            targ.kern.type=targ.kern.type, 
            get_cross_kerns=T, get_avg_kerns=T, 
            kern.interp.xs=kern.interp.xs,
            subtract.mean=F, dM=model$M-star.M, dR=model$R-star.R, 
            perturb=F, num_realizations=1, F_surf=model$F_surf)
    }, model_i=1:length(model.list))
    names(inv.list) <- model.names
    
    # save results 
    inv.lists[[trial_i]] <- Map(function(inv) inv$result, inv=inv.list)
    avg.kerns.lists[[trial_i]] <- Map(function(inv) inv$avg_kern, inv=inv.list)
    cross.lists[[trial_i]] <- Map(function(inv) inv$cross_kern, inv=inv.list)
    
    star.Ms <- c(star.Ms, star.M)
    star.Rs <- c(star.Rs, star.R)
    MRs <- list(Ms=star.Ms, Rs=star.Rs)
    
    cross.terms <- c(cross.terms, cross.term)
    error.sups  <- c(error.sups,  error.sup)
    widths      <- c(widths,      width)
    inv.params  <- list(
        cross.terms=cross.terms, 
        error.sups=error.sups,
        widths=widths)
    
    save(ref.mod, 
        file=file.path('save', paste0('ref.mod',         targ.mode)))
    save(inv.lists, 
        file=file.path('save', paste0('inv.lists',       targ.mode)))
    save(avg.kerns.lists, 
        file=file.path('save', paste0('avg.kerns.lists', targ.mode)))
    save(cross.lists, 
        file=file.path('save', paste0('cross.lists',     targ.mode)))
    save(MRs, 
        file=file.path('save', paste0('MRs',             targ.mode)))
    save(inv.params, 
        file=file.path('save', paste0('inv.params',      targ.mode)))
    
    inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
        inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
        cross.lists=cross.lists, inv.params=inv.params,
        kern.interp.xs=kern.interp.xs)
    
    make_plots_inversion_all(ref.mod, inversion, kern.interp.xs=kern.interp.xs,
        k.str=targ.mode, cross.inset="bottomright",
        #inversion_ylim=c(-0.15, 0.05),
        #col.pal="#F46D43", sampler=c(F,F,T,F,F,F),
        #inversion_ylim=c(0, 0.15),
        cross_kern_ylim=c(-0.8, 0.3))
}

print_latex_table(inversion)

if (F) {

    load(file.path('save', paste0('ref.mod',         targ.mode)))
    load(file.path('save', paste0('inv.lists',       targ.mode)))
    load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
    load(file.path('save', paste0('cross.lists',     targ.mode)))
    load(file.path('save', paste0('MRs',             targ.mode)))
    load(file.path('save', paste0('inv.params',      targ.mode)))
    
    
    
    make_plots(plot_ref_mods, paste0('ref_mods', targ.mode),
        model.list=model.list, inversion=inversion, 
        xlim=c(0, 0.35), ylim=c(0.7, 1.7)) 
    
}

