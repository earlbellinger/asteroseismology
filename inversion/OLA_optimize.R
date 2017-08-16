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
num_procs <- max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

args <- commandArgs(TRUE)
MOLA            <- if (length(args)>0) as.logical(as.numeric(args[1])) else F
mode.set        <- if (length(args)>1)   args[2] else 'CygA'
error.set       <- if (length(args)>2)   args[3] else 'CygA'
target.name     <- if (length(args)>3)   args[4] else 'modmix'
targ.kern.type  <- if (length(args)>4)   args[5] else 'mod_Gauss'
n_trials        <- if (length(args)>5)   as.numeric(args[6]) else 128
model.list.name <- if (length(args)>6)   args[7] else 'perturbed.model.names'
#'perturbed.model.names'

half <- F
k.pair <- u_Y
rs <- seq(0.08, 0.3, 0.02)
targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set,
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type)) 

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
    error.set=error.set, perturb=T) 

inversion.lists <- list()
avg.kerns.lists <- list()
cross.kerns.lists <- list()

models <- get_model_list() 
model.names <- get(model.list.name) 

model.list <- list()
for (perturbed.model.name in model.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=T) 
    #model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
    #    nu_ac=model$nu_ac)
    model.list[[perturbed.model.name]] <- model
}

ref.mod <- model.list[[3]] 

ref.mod2.name <- 
    if (target.name=='CygA') 'CygAwball' else 
    if (target.name=='CygB') 'CygBwball' else 
    if (target.name=='BiSON') 'diffusion' 

ref.mod2 <- if (!(target.name %in% c('CygA', 'CygB', 'BiSON'))) ref.mod else
    get_model(freqs=freqs, model.name=ref.mod2.name, target.name=target.name, 
            k.pair=k.pair, square.Ks=F) 


for (trial_i in 1:n_trials) {
    
    print(paste("TRIAL_I =", trial_i))
    
    for (perturbed.model.name in model.names) {
        model <- model.list[[perturbed.model.name]]
        model$nus$nu.y <- rnorm(length(freqs$nu), freqs$nu, freqs$dnu)
        model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
            nu_ac=model$nu_ac)
        model.list[[perturbed.model.name]] <- model
    }
    
    if (MOLA) {
        #######################
        ### MOLA ##############
        #######################
        
        initial_params <- rep(10, 2*length(model.list))
        parscale <- rep(10, 2*length(model.list))
        
        iteration <<- 0
        best_result <<- Inf
        best_params <<- initial_params
        
        best_params <- optim(log10(initial_params), 
                fn=function(inversion_params) {
            inversion_params <- 10**inversion_params
            
            iteration <<- iteration + 1
            cat(paste("**** iter:", iteration, "\n"))
            
            cat('Trying inversion params ')
            cat(inversion_params)
            cat('\n')
            
            inversion.list <- list()
            inversion.list <- parallelMap(function(model_i) {
                cross.term <- inversion_params[model_i]
                error.sup <- inversion_params[length(model.list) + model_i]
                
                model <- model.list[[model_i]]
                cat(paste("Inverting model", model.names[model_i], 
                    "with params:",
                    "beta =", format(cross.term, digits=4, nsmall=3), 
                    "; mu =", format(error.sup, digits=4, nsmall=3), 
                    "\n"))
                inversion <- invert.OLA(model=model, rs=rs,
                           cross.term=cross.term, error.sup=error.sup)
                
                inversion.list[[model.names[model_i]]] <- inversion
            }, model_i <- 1:length(model.list))
            names(inversion.list) <- model.names
            
            f.stds <- sapply(inversion.list, 
                function(inversion) inversion$result$f.err)
            f.means <- sapply(inversion.list, 
                function(inversion) inversion$result$f)
            
            w.f.stds <- apply(replicate(100, {
                    f.means2 <- sapply(1:ncol(f.means), function(ii) 
                        rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii]))
                    apply(f.means2, 1, sd)
                }), 1, mean)
            
            result <- sum(w.f.stds)
            cat(paste("Result:", format(result, digits=4), '\n'))
            
            if (result < best_result) {
                best_result <<- result
                best_params <<- inversion_params
                cat(paste(c("*****", inversion_params, "\n")))
                cat("***** New record!\n")
            }
            
            result
        }, control=list(trace=999, parscale=parscale))

        print(best_params)

        inversion_params <- 10**best_params$par
        inversion.list <- list()
        inversion.list <- parallelMap(function(model_i) {
            cross.term <- inversion_params[model_i]
            error.sup <- inversion_params[length(model.list) + model_i]
            
            model <- model.list[[model_i]]
            cat(paste("Inverting model", model.names[model_i], 
                "with params:",
                "beta =", format(cross.term, digits=4, nsmall=3), 
                "; mu =", format(error.sup, digits=4, nsmall=3), 
                "\n"))
            inversion <- invert.OLA(model=model, rs=rs,
                       cross.term=cross.term, error.sup=error.sup)
            
            inversion.list[[model.names[model_i]]] <- inversion
        }, model_i <- 1:length(model.list))
        names(inversion.list) <- model.names
        
        
        
        
        
    } else {
        #######################
        ### SOLA ##############
        #######################
        
        initial_params <- c(rep(1000, 2*length(model.list)), 
                            rep(0.02,   length(model.list)))
        parscale <- c(rep(10, 2*length(model.list)), 
                      rep(1,    length(model.list)))
        
        iteration <<- 0
        best_result <<- Inf
        best_params <<- initial_params
        
        SOLA_optim <- function(inversion_params) {
            inversion_params <- 10**inversion_params
            
            if (any(inversion_params[-1:-(2*length(model.list))] > 0.08)) 
                return(Inf)
            
            iteration <<- iteration + 1
            cat(paste("**** iter:", iteration, "\n"))
            
            cat('Trying inversion params ')
            cat(inversion_params)
            cat('\n')
            
            inversion.list <- list()
            inversion.list <- parallelMap(function(model_i) {
                cross.term <- inversion_params[model_i]
                error.sup <- inversion_params[length(model.list) + model_i]
                width <- inversion_params[2*length(model.list) + model_i]
                
                model <- model.list[[model_i]]
                cat(paste("Inverting model", model.names[model_i], 
                    "with params:",
                    "beta =", format(cross.term, digits=4, nsmall=3), 
                    "; mu =", format(error.sup, digits=4, nsmall=3), 
                    "; width =", format(width, digits=4, nsmall=3), "\n"))
                invert.OLA(model=model, rs=rs, 
                    cross.term=cross.term, error.sup=error.sup, width=width,
                    targ.kern.type=targ.kern.type, 
                    get_cross_kerns=F, get_avg_kerns=F, 
                    perturb=F, num_realizations=1, F_surf=model$F_surf)$result
            }, model_i=1:length(model.list))
            names(inversion.list) <- model.names
            
            f.means <- sapply(inversion.list, function(inversion) 
                inversion$f)
            result <- sum(apply(f.means, 1, sd))
            cat(paste("Result:", format(result, digits=4), '\n'))
            
            if (result < best_result) {
                best_result <<- result
                best_params <<- inversion_params
                cat(paste(c("*****", inversion_params, "\n")))
                cat("***** New record!\n")
            }
            
            result
        }
        
        best_params <- optim(log10(initial_params), fn=SOLA_optim, 
            control=list(trace=999, parscale=parscale, 
            maxit=512))#, abstol=sum(model.list[[1]]$f1.spl(rs)/1000)))
        
        print(best_params)
        
        inversion_params <- 10**best_params$par
        inversion.list <- parallelMap(function(model_i) {
            cross.term <- inversion_params[model_i]
            error.sup <- inversion_params[length(model.list) + model_i]
            width <- inversion_params[2*length(model.list) + model_i]
            
            model <- model.list[[model_i]]
            cat(paste("Inverting model", model.names[model_i], 
                "with params:",
                "beta =", format(cross.term, digits=4, nsmall=3), 
                "; mu =", format(error.sup, digits=4, nsmall=3), 
                "; width =", format(width, digits=4, nsmall=3), "\n"))
            invert.OLA(model=model, rs=rs, 
                cross.term=cross.term, error.sup=error.sup, width=width,
                targ.kern.type=targ.kern.type, 
                get_cross_kerns=T, get_avg_kerns=T, 
                perturb=F, num_realizations=1, F_surf=model$F_surf)
        }, model_i=1:length(model.list))
        names(inversion.list) <- model.names
        
        #model.lists[[trial_i]] <- model.list 
        inversion.lists[[trial_i]] <- Map(function(inversion) 
            inversion$result, inversion=inversion.list)
        avg.kerns.lists[[trial_i]] <- Map(function(inversion) 
            inversion$avg_kern, inversion=inversion.list)
        cross.kerns.lists[[trial_i]] <- Map(function(inversion) 
            inversion$cross_kern, inversion=inversion.list)
        #inversion.list 
        #inversion.lists[[trial_i]] <- inversion.list
        
        save(ref.mod, file=paste0("save/ref.mod", targ.mode))
        save(inversion.lists, file=paste0("save/inversion.lists", targ.mode))
        save(avg.kerns.lists, file=paste0("save/avg.kerns.lists", targ.mode))
        save(cross.kerns.lists, file=paste0("save/cross.kerns.lists", 
            targ.mode))

    }
}

#save(model.lists, file="model.lists")

#sampler <- c(F, T)#c(F, T, F, T, F, T, F, T, F)

sampler <- T #c(F, F, T, F, F, T, F, T, F, T, F, F) #
#sampler <- c(F, F, F, F, F, T, T, T, F, F, F, F)#c(F, T)

#plot_inversion_lists_mean(model=ref.mod, inversion.lists=inversion.lists);dev.off;
make_plots(plot_inversion_lists_mean,
    filename=paste0("inv-lists-optim", targ.mode), 
    model=ref.mod2, inversion.lists=inversion.lists, legend.spot=NULL,
    #ylim=c(-0.03, 0.03), #c(-0.03, 0.015),#
    xlim=c(0, 0.45), 
    #ylim=c(-0.04, 0.1),
    k.pair=k.pair, sampler=sampler)#c(T))#

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim-kerns", targ.mode), 
    kernel.lists=avg.kerns.lists, cross=F, legend.spot="topleft",
    xlim=c(0, 0.45),
    ylim=c(-7, 25),
    model=ref.mod, k.pair=k.pair, sampler=sampler)#c(T))#

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim-cross", targ.mode), 
    kernel.lists=cross.kerns.lists, cross=T, legend.spot="topleft",
    xlim=c(0, 0.45),
    ylim=c(-0.2, 1),
    model=ref.mod, k.pair=k.pair, sampler=sampler)#c(T))#

#for (ref.mod in model.list) plot_inversion_lists_mean(ref.mod=ref.mod, 
#    inversion.lists=inversion.lists)
#plot_inv_diffs_mean(model.list, inversion.list)
#dev.off()

#make_plots(plot_inv_diffs_mean,
#    filename=paste0("inv-diffs-mean-optim", targ.mode), 
#    model.list=model.list, inversion.list=inversion.list)

#make_plots(plot_inversion_mean,
#    filename=paste0("inv-mean-optim", targ.mode), 
#    model=model.list[[5]], inversion.list=inversion.list)

if (F) {

    load(paste0("save/ref.mod", targ.mode))
    load(paste0("save/inversion.lists", targ.mode))
    load(paste0("save/avg.kerns.lists", targ.mode))
    load(paste0("save/cross.kerns.lists", targ.mode))

}

