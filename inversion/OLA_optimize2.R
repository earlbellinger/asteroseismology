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
MOLA <- if (length(args)>0) as.logical(as.numeric(args[1])) else F
mode.set <- if (length(args)>1) args[2] else 'CygA'
error.set <- if (length(args)>2) args[3] else 'CygA'
target.name <- if (length(args)>3) args[4] else 'CygAwball'
targ.kern.type <- if (length(args)>4) args[5] else 'mod_Gauss'
n_trials <- if (length(args)>5) args[6] else 1
model.list.name <- if (length(args)>6) args[7] else 'perturbed.CygA.names'
perturb <- if (length(args)>7) as.logical(as.numeric(args[8])) else F
k.pair <- u_Y
rs <- seq(0.08, 0.35, 0.01)
targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set,
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type))


#model.lists <- list()
inversion.lists <- list()
avg.kerns.lists <- list()
cross.kerns.lists <- list()

ref.mod <- NULL
models <- get_model_list()
model.names <- get(model.list.name)

for (trial_i in 1:n_trials) {

print(paste("TRIAL_I =", trial_i))

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
    error.set=error.set, perturb=perturb) 

model.list <- list()
for (perturbed.model.name in model.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=T) 
    model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
        nu_ac=model$nu_ac)
    model.list[[perturbed.model.name]] <- model
}

if (trial_i == 1) ref.mod <- model.list[[3]]#5]] 



if (MOLA) {
#######################
### MOLA ##############
#######################

initial_params <- rep(10, 2*length(model.list))
parscale <- rep(10, 2*length(model.list))

iteration <<- 0
best_result <<- Inf
best_params <<- initial_params

best_params <- optim(log10(initial_params), fn=function(inversion_params) {
    inversion_params <- 10**inversion_params
    
    #width <- inversion_params[length(inversion_params)]
    #if (width > min(diff(rs))) return(Inf) 
    #if (any(inversion_params[-1:-(2*length(model.list))] > 0.05)) return(Inf)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    cat('Trying inversion params ')
    cat(inversion_params)
    cat('\n')
    
    inversion.list <- list()
    #for (model_i in 1:length(model.list)) {
    inversion.list <- parallelMap(function(model_i) {
        cross.term <- inversion_params[model_i]
        error.sup <- inversion_params[length(model.list) + model_i]
        #width <- inversion_params[2*length(model.list) + model_i]
        
        model <- model.list[[model_i]]
        cat(paste("Inverting model", model.names[model_i], 
            "with params:",
            "beta =", format(cross.term, digits=4, nsmall=3), 
            "; mu =", format(error.sup, digits=4, nsmall=3), 
            #"; width =", format(width, digits=4, nsmall=3), 
            "\n"))
        inversion <- invert.OLA(model=model, rs=rs,
                   cross.term=cross.term, error.sup=error.sup)
        
        inversion.list[[model.names[model_i]]] <- inversion
    }, model_i <- 1:length(model.list))
    names(inversion.list) <- model.names
    
    f.stds <- sapply(inversion.list, function(inversion) inversion$result$f.err)
    f.means <- sapply(inversion.list, function(inversion) inversion$result$f)
    
    w.f.stds <- apply(replicate(100, {
            f.means2 <- sapply(1:ncol(f.means), function(ii) 
                rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii]))
            #f.means2 <- apply(f.means, 2, function(f.mean) 
            #    rnorm(length(f.mean), f.mean, f.std))
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
#width <- inversion_params[length(inversion_params)]
inversion.list <- list()
inversion.list <- parallelMap(function(model_i) {
    cross.term <- inversion_params[model_i]
    error.sup <- inversion_params[length(model.list) + model_i]
    
    model <- model.list[[model_i]]
    cat(paste("Inverting model", model.names[model_i], 
        "with params:",
        "beta =", format(cross.term, digits=4, nsmall=3), 
        "; mu =", format(error.sup, digits=4, nsmall=3), 
        #"; width =", format(width, digits=4, nsmall=3), 
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

initial_params <- c(rep(10,   2*length(model.list)), 
                    rep(0.02, length(model.list)))
parscale <- c(rep(10, 2*length(model.list)), 
              rep(1,    length(model.list)))

iteration <<- 0
best_result <<- Inf
best_params <<- initial_params

SOLA_optim <- function(inversion_params) {
    inversion_params <- 10**inversion_params
    
    if (any(inversion_params[-1:-(2*length(model.list))] > 0.05)) return(Inf)
    
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
                   targ.kern.type=targ.kern.type, get_cross_kerns=F,
                   F_surf=model$F_surf)$result
    }, model_i=1:length(model.list))
    names(inversion.list) <- model.names
    
    f.stds <- sapply(inversion.list, function(inversion) inversion$f.err)
    f.means <- sapply(inversion.list, function(inversion) inversion$f)
    w.f.stds <- apply(do.call(cbind, parallelMap(function(ii) {
        replicate(floor(10**5/num_procs),
            apply(sapply(1:ncol(f.means), function(ii) 
                rnorm(length(f.means[,ii]), f.means[,ii], f.stds[,ii])), 1, sd)
        )}, ii=1:num_procs)), 1, mean)
    result <- sum(w.f.stds)
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
    maxit=10000)) 

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
               targ.kern.type=targ.kern.type, F_surf=model$F_surf)
}, model_i <- 1:length(model.list))
names(inversion.list) <- model.names


inversion.lists[[trial_i]] <- Map(function(inversion) inversion$result, 
    inversion=inversion.list)
avg.kerns.lists[[trial_i]] <- Map(function(inversion) inversion$avg_kern, 
    inversion=inversion.list)
cross.kerns.lists[[trial_i]] <- Map(function(inversion) inversion$cross_kern, 
    inversion=inversion.list)

save(ref.mod, file=paste0("save3/ref.mod", targ.mode))
save(inversion.lists, file=paste0("save3/inversion.lists", targ.mode))
save(avg.kerns.lists, file=paste0("save3/avg.kerns.lists", targ.mode))
save(cross.kerns.lists, file=paste0("save3/cross.kerns.lists", targ.mode))

}
}

ref.mod2 <- if (target.name=='CygA') get_model(freqs=freqs, 
    model.name='CygAwball', 
    target.name='CygA', 
    k.pair=k.pair, square.Ks=F) else ref.mod 

sampler <- c(T) 

make_plots(plot_inversion_lists_mean,
    filename=paste0("inv-lists-optim3", targ.mode), 
    model=ref.mod2, inversion.lists=inversion.lists, legend.spot="topleft",
    k.pair=k.pair, sampler=sampler)

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim3-kerns", targ.mode), 
    kernel.lists=avg.kerns.lists, cross=F, legend.spot="topleft",
    model=ref.mod, k.pair=k.pair, sampler=sampler)

make_plots(plot_kernel_lists,
    filename=paste0("inv-lists-optim3-cross", targ.mode), 
    kernel.lists=cross.kerns.lists, cross=T, legend.spot="topleft",
    model=ref.mod, k.pair=k.pair, sampler=sampler)


