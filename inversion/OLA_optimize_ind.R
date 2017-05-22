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
library(matrixStats)
parallelStartMulticore(8)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))#9

args <- commandArgs(TRUE)
MOLA <- if (length(args)>0) as.logical(as.numeric(args[1])) else T
mode.set <- if (length(args)>1) args[2] else 'CygA'
error.set <- if (length(args)>2) args[3] else 'CygA'
target.name <- if (length(args)>3) args[4] else 'modmix'
targ.kern.type <- if (length(args)>4) args[5] else 'mod_Gauss'
n_trials <- if(length(args)>5) args[6] else 32
half <- F
k.pair <- u_Y
rs <- seq(0.05, 0.29, 0.06)

targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set,
    '-n_', n_trials, 
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type))

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
        error.set=error.set, perturb=T) 
perturbed.model.name <- perturbed.CygA.names[5]
ref.mod <- get_model(freqs=freqs, 
    model.name=perturbed.model.name, 
    target.name=target.name, 
    k.pair=k.pair, square.Ks=T) 
model.names <- c(perturbed.CygA.names[1], perturbed.CygA.names[1])
#perturbed.model.names[9])

invert <- function(model_i, model.list, inversion_params, r, targ.kern.type,
        MOLA, MOLA.K_ijs=NA, MOLA.C_ijs=NA) {
    cross.term <- inversion_params[model_i]
    error.sup <- inversion_params[length(model.list) + model_i]
    width <- if (MOLA) NULL else 
        inversion_params[2*length(model.list) + model_i]
    
    model <- model.list[[model_i]]
    cat(paste("Inverting model", names(model.list)[model_i], 
        "at r =", r,
        "with params:",
        if (MOLA) paste(c("beta =", "mu ="), 
            format(c(cross.term, error.sup), digits=4, nsmall=3)) else
                  paste(c("beta =", "mu =", "width ="),
            format(c(cross.term, error.sup, width), digits=4, nsmall=3)),
        #"beta =", format(cross.term, digits=4, nsmall=3), 
        #"; mu =", format(error.sup, digits=4, nsmall=3), 
        "\n"))
    invert.OLA(model=model, rs=r,
        cross.term=cross.term, 
        error.sup=error.sup,
        width=width,
        targ.kern.type=targ.kern.type,
        MOLA.K_ijs=MOLA.K_ijs,
        MOLA.C_ijs=MOLA.C_ijs)
}

pair_optim <- function(log10_inversion_params, model.list, 
                       r, targ.kern.type, MOLA, 
                       MOLA.K.list=list(), MOLA.C.list=list()) {
    inversion_params <- 10**log10_inversion_params
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    cat('Trying inversion params ')
    cat(inversion_params)
    cat('\n')
    
    #inversion.list <- list()
    inversion.list <- parallelMap(function(model_i) {
        invert(model_i, model.list, inversion_params, r, targ.kern.type, MOLA,
            MOLA.K_ijs=if(MOLA) MOLA.K.list[[model_i]] else NA,
            MOLA.C_ijs=if(MOLA) MOLA.C.list[[model_i]] else NA)
    }, model_i=1:length(model.list))
    
    if (with(inversion.list[[1]]$result, fwhm.left > rs | fwhm.right < rs |
                                         r.first_q > rs | r.third_q  < rs) | 
        with(inversion.list[[2]]$result, fwhm.left > rs | fwhm.right < rs |
                                         r.first_q > rs | r.third_q  < rs)) {
        result <- Inf
    } else {
        result <- (inversion.list[[1]]$result$f - 
                   inversion.list[[2]]$result$f)**2 #+ 
                  #(with(inversion.list[[1]]$result, fwhm.mid - rs)**2 +
                  # with(inversion.list[[2]]$result, fwhm.mid - rs)**2)
    }
    cat(paste("Result:", format(result, digits=4), '\n'))
    
    if (result < best_result) {
        best_result <<- result
        best_params <<- inversion_params
        cat(paste(c("*****", inversion_params, "\n")))
        cat("***** New record!\n")
    }
    
    result
}

trial <- function(r, model.names, target.name, mode.set, 
                  error.set, k.pair, targ.kern.type, MOLA) {
    model.list <- list()
    MOLA.K.list <- list()
    MOLA.C.list <- list()
    for (ii in 1:length(model.names)) {
        freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
            error.set=error.set, perturb=T) 
        model <- get_model(freqs=freqs, 
            model.name=model.names[[ii]], 
            target.name=target.name, 
            k.pair=k.pair, square.Ks=T) 
        model.list[[ii]] <- model
        if (MOLA) MOLA.K.list[[ii]] <- get_square_Ks(model$modes, model$k1, r)
        if (MOLA) MOLA.C.list[[ii]] <- get_square_Ks(model$modes, model$k2, r)
    }
    
    inversion_params <- rep(10, 2*length(model.list))
    parscale <- rep(100, 2*length(model.list))
    if (!MOLA) {
        inversion_params <- c(inversion_params, rep(0.02, length(model.list)))
        parscale <- c(parscale, rep(1, length(model.list)))
    }
    
    iteration <<- 0
    best_result <<- Inf
    best_params <<- inversion_params
    
    #MOLA.K_ijs <- if (MOLA) get_square_Ks(model$modes, model$k1, r) else NA
    #MOLA.C_ijs <- if (MOLA) get_square_Ks(model$modes, model$k2, r) else NA
    
    best_params <- optim(log10(inversion_params), fn=pair_optim, 
        model.list=model.list, r=r, targ.kern.type=targ.kern.type, MOLA=MOLA,
        MOLA.K.list=MOLA.K.list, MOLA.C.list=MOLA.C.list,
        control=list(trace=999, parscale=parscale))
    print(best_params)
    
    inversion_params <- 10**best_params$par
    do.call(rbind, parallelMap(function(model_i) {
        invert(model_i, model.list, inversion_params, r, targ.kern.type, MOLA,
            MOLA.K_ijs=if(MOLA) MOLA.K.list[[model_i]] else NA,
            MOLA.C_ijs=if(MOLA) MOLA.C.list[[model_i]] else NA)$result
    }, model_i <- 1:length(model.list)))
}

results.lists <- list()
for (r in rs) {
    results.list <- do.call(rbind, parallelMap(function(trial_i) {
        cat(paste("TRIAL_I =", trial_i, '\n'))
        trial(r, model.names, target.name, mode.set, error.set, k.pair,
              targ.kern.type, MOLA)
    }, trial_i=1:n_trials))
    results.lists[[paste(r)]] <- results.list
}

save(ref.mod, file=paste0("ref.mod", targ.mode))
save(results.lists, file=paste0("results.lists", targ.mode))

#make_plots(plot_results_lists_mean_points,
#    filename=paste0("res-lists-pair-optim", targ.mode), 
#    model=ref.mod, results.lists=results.lists)

make_plots(plot_results_lists,
    filename=paste0("res-lists-pair-optim-mean", targ.mode), 
    model=ref.mod, results.lists=results.lists)

