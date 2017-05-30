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
parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))#9

args <- commandArgs(TRUE)
MOLA <- if (length(args)>0) as.logical(as.numeric(args[1])) else T
mode.set <- if (length(args)>1) args[2] else 'CygA'
error.set <- if (length(args)>2) args[3] else 'CygA'
target.name <- if (length(args)>3) args[4] else 'CygAwball'
targ.kern.type <- if (length(args)>4) args[5] else 'mod_Gauss'
n_trials <- if(length(args)>5) args[6] else 32
half <- F
k.pair <- u_Y
rs <- seq(0.05, 0.29, 0.06)

targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set,
    '-n_', n_trials, 
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type))

print(targ.mode)

models <- get_model_list()
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
        error.set=error.set, perturb=T) 
perturbed.model.name <- perturbed.CygA.names[3]
ref.mod <- get_model(freqs=freqs, 
    model.name=perturbed.model.name, 
    target.name=target.name, 
    k.pair=k.pair, square.Ks=F) 
model.names <- perturbed.CygA.names
#c(perturbed.CygA.names[1], perturbed.CygA.names[1])
#perturbed.model.names[9])

#cat <- function(x) if (F) cat(x)

print_params <- function(model_i, r, cross.term, error.sup, width, MOLA) {
    cat(paste("Inverting model", model_i, 
        "at r =", r,
        "with params: "))
    cat(if (MOLA) paste(c("beta =", "mu ="), 
            format(c(cross.term, error.sup), digits=4, nsmall=3)) else
                  paste(c("beta =", "mu =", "width ="),
            format(c(cross.term, error.sup, width), digits=4, nsmall=3)),
        "\n")
}

optim_fn <- function(log10_inversion_params, model.list, 
                       r, targ.kern.type, MOLA, 
                       MOLA.K.list=list(), MOLA.C.list=list()) {
    inversion_params <- 10**log10_inversion_params
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    cat('Trying inversion params ')
    cat(inversion_params)
    cat('\n')
    
    inversion.list <- parallelMap(function(model_i) {
        cross.term <- inversion_params[model_i]
        error.sup <- inversion_params[length(model.list) + model_i]
        width <- if (MOLA) NULL else 
            inversion_params[2*length(model.list) + model_i]
        
        memo <- paste(model_i, r, cross.term, error.sup, width, sep='_')
        if (memo %in% names(result.list)) {
            cat(paste0("Returning memoized value for ", memo, "\n"))
            return(result.list[[memo]])
        }
        model <- model.list[[model_i]]
        
        print_params(model_i, r, cross.term, error.sup, width, MOLA)
        inv.res <- invert.OLA(model=model, rs=r,
            cross.term=cross.term, 
            error.sup=error.sup,
            width=width,
            targ.kern.type=targ.kern.type,
            MOLA.K_ijs=model$K_ijs,#MOLA.K.list[[model_i]],
            MOLA.C_ijs=model$C_ijs,#MOLA.C.list[[model_i]],
            F_surf=model$F_surf)$result
        
        result.list[[memo]] <<- inv.res
        inv.res
    }, model_i=1:length(model.list))
    
    for (ii in 1:length(inversion.list)) {
        if (with(inversion.list[[ii]],
                fwhm.left > rs | fwhm.right < rs)) {
            cat(paste0(ii, " out of bounds\n"))
            return(Inf)
        }
    }
    
    result <- if (length(inversion.list) == 2) {
        (inversion.list[[1]]$f - inversion.list[[2]]$f)**2
    } else {
        f.stds <- sapply(inversion.list, function(inv.) inv.$f.err)
        f.means <- sapply(inversion.list, function(inv.) inv.$f)
        mean(replicate(10000, sd(rnorm(length(f.means), f.means, f.stds))))
        #sd(sapply(inversion.list, function(inv.) inv.$f))
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

trial <- function(r, model.names, target.name, mode.set, error.set, k.pair, 
        targ.kern.type, MOLA, MOLA.K.list, MOLA.C.list, K.ints.list) {
    model.list <- parallelMap(function(ii) {
        model <- get_model(freqs=get_freqs(target.name=target.name, 
                mode.set=mode.set, error.set=error.set, perturb=T), 
            model.name=model.names[ii], 
            target.name=target.name, 
            k.pair=k.pair, square.Ks=F)
        if (MOLA) model$K_ijs <- MOLA.K.list[[ii]]
        if (MOLA) model$C_ijs <- MOLA.C.list[[ii]]
        model$K.ints <- K.ints.list[[ii]]
        model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
            nu_ac=model$nu_ac)
        model
    }, ii=1:length(model.names))
    
    inversion_params <- rep(10, 2*length(model.list))
    parscale <- rep(100, 2*length(model.list))
    if (!MOLA) {
        inversion_params <- c(inversion_params, rep(0.02, length(model.list)))
        parscale <- c(parscale, rep(1, length(model.list)))
    }
    
    iteration <<- 0
    best_result <<- Inf
    best_params <<- inversion_params
    result.list <<- list()
    
    best_params <- optim(log10(inversion_params), fn=optim_fn, 
        model.list=model.list, r=r, targ.kern.type=targ.kern.type, MOLA=MOLA,
        MOLA.K.list=MOLA.K.list, MOLA.C.list=MOLA.C.list,
        control=list(trace=999, parscale=parscale, maxit=10000, 
            abstol=model.list[[1]]$f1.spl(r)/10000))
    print(best_params)
    
    inversion_params <- 10**best_params$par
    parallelMap(function(model_i) {
        model <- model.list[[model_i]]
        cross.term <- inversion_params[model_i]
        error.sup <- inversion_params[length(model.list) + model_i]
        width <- if (MOLA) NULL else 
            inversion_params[2*length(model.list) + model_i]
        print_params(model_i, r, cross.term, error.sup, width, MOLA)
        inv.res <- invert.OLA(model=model, rs=r,
            cross.term=cross.term, error.sup=error.sup, width=width,
            targ.kern.type=targ.kern.type, 
            MOLA.K_ijs=model$K_ijs,
            MOLA.C_ijs=model$C_ijs,
            F_surf=model$F_surf)
    }, model_i <- 1:length(model.list))
}

inversion.lists <- list()
avg.kerns.lists <- list()
cross.kerns.lists <- list()

for (r in rs) {
    lists <- parallelMap(function(ii) {
        model <- get_model(freqs=freqs, 
            model.name=model.names[[ii]], 
            target.name=target.name, 
            k.pair=k.pair, square.Ks=T, x0=r) 
        list(MOLA.K=model$K_ijs, MOLA.C=model$C_ijs, K.ints=model$K.ints)
    }, ii=1:length(model.names))
    MOLA.K.list <- Map(function(list.) list.[[1]], list.=lists) 
    MOLA.C.list <- Map(function(list.) list.[[2]], list.=lists) 
    K.ints.list <- Map(function(list.) list.[[3]], list.=lists) 
    
    trials <- parallelMap(function(trial_i) {
        cat(paste("TRIAL_I =", trial_i, '\n'))
        trial(r, model.names, target.name, mode.set, 
            error.set, k.pair, targ.kern.type, MOLA, 
            MOLA.K.list, MOLA.C.list, K.ints.list)
    }, trial_i=1:n_trials)
    
    avg.kerns.list <- Map(function(invs) do.call(cbind, 
        Map(function(inv.) inv.$avg_kern, inv.=invs)), invs=trials)
    cross.kerns.list <- Map(function(invs) do.call(cbind, 
        Map(function(inv.) inv.$cross_kern, inv.=invs)), invs=trials)
    results.list <- Map(function(invs) do.call(rbind, 
        Map(function(inv.) inv.$result, inv.=invs)), invs=trials)
    
    inversion.lists[[paste(r)]] <- results.list
    avg.kerns.lists[[paste(r)]] <- avg.kerns.list
    cross.kerns.lists[[paste(r)]] <- cross.kerns.list
}

save(inversion.lists, file=paste0("save/ind_results.lists", targ.mode))
save(avg.kerns.lists, file=paste0("save/ind_avg.kerns.lists", targ.mode))
save(cross.kerns.lists, file=paste0("save/ind_cross.kerns.lists", targ.mode))

make_plots(plot_inversion_lists_mean,
    filename=paste0("inv-lists-optim-ind", targ.mode), 
    model=ref.mod, inversion.lists=inversion.lists, legend.spot="topleft",
    #ylim=c(0, 0.05), 
    #xlim=c(0.05, 0.35), 
    k.pair=k.pair, sampler=c(T)) 

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim-ind-kerns", targ.mode), 
    kernel.lists=avg.kerns.list, cross=F, legend.spot="topleft",
    model=ref.mod, k.pair=k.pair) 

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim-ind-cross", targ.mode), 
    kernel.lists=cross.kerns.list, cross=T, legend.spot="topleft",
    model=ref.mod, k.pair=k.pair) 


#make_plots(plot_results_lists_mean_points,
#    filename=paste0("res-lists-pair-optim", targ.mode), 
#    model=ref.mod, inversion.lists=inversion.lists)

#make_plots(plot_results_lists,
#    filename=paste0("res-lists-pair-optim-mean", targ.mode), 
#    model=ref.mod, inversion.lists=inversion.lists)

