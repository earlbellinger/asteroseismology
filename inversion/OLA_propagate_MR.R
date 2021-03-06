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
MOLA            <- if (length(args)>0)  as.logical(as.numeric(args[1])) else F
mode.set        <- if (length(args)>1)  args[2] else 'CygA'
error.set       <- if (length(args)>2)  args[3] else 'CygA'
ref.mod.name    <- if (length(args)>3)  args[4] else 'CygADiffusion'
target.name     <- if (length(args)>4)  args[5] else 'CygAlowRhighM'
targ.kern.type  <- if (length(args)>5)  args[6] else 'mod_Gauss'
n_trials        <- if (length(args)>6)  as.numeric( args[7]) else 128
model.list.name <- if (length(args)>7)  args[8] else 'perturbed.model.names'
initial.M       <- if (length(args)>8)  as.numeric( args[9]) else 1
initial.R       <- if (length(args)>9)  as.numeric(args[10]) else 1
sigma.M         <- if (length(args)>10) as.numeric(args[11]) else 0.016
sigma.R         <- if (length(args)>11) as.numeric(args[12]) else 0.02

perturb <- F 
half <- F 
k.pair <- u_Y 
rs <- seq(0.05, 0.3, 0.05)
sampler <- T 
targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set, 
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type)) 

freqs <- get_freqs(target.name=target.name, mode.set=mode.set, 
    error.set=error.set, perturb=perturb) 

inv.lists <- list()
avg.kerns.lists <- list()
cross.lists <- list()
star.Ms <- c()
star.Rs <- c()
dir.create("save", showWarnings = FALSE)

models <- get_model_list() 
model.names <- get(model.list.name) 

model.list <- list()
for (perturbed.model.name in model.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=T) 
    model.list[[perturbed.model.name]] <- model
}

ref.mod <- model.list[[5]] 


for (trial_i in 1:n_trials) {
    
    print(paste("TRIAL_I =", trial_i))
    nu.star <- rnorm(length(freqs$nu), freqs$nu, freqs$dnu)
    trial.M <- rnorm(1, initial.M, sigma.M)
    trial.R <- rnorm(1, initial.R, sigma.R)
    
    for (perturbed.model.name in model.names) {
        model <- model.list[[perturbed.model.name]]
        model$nus$nu.y <- nu.star
        model$F_surf <- get_F_surf(model$nus, num.knots=0, use.BG=T, 
            nu_ac=model$nu_ac)
        model.list[[perturbed.model.name]] <- model
    }
    
    initial_params <- c(1000, 1000, 0.02, trial.M, trial.R)
    parscale <- c(0.1, 0.1, 0.01, sigma.M, sigma.R)
    
    iteration <<- 0
    best_result <<- Inf
    best_params <<- initial_params
    
    SOLA_optim <- function(inversion_params) {
        inversion_params <- 10**inversion_params
        
        cross.term <- inversion_params[1]
        error.sup <- inversion_params[2]
        width <- inversion_params[3]
        star.M <- inversion_params[4]
        star.R <- inversion_params[5]
        if (star.M < initial.M-10*sigma.M || 
            star.M > initial.M+10*sigma.M || 
            star.R < initial.R-10*sigma.R || 
            star.R > initial.R+10*sigma.R ||
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
        result <- sum(apply(f.means, 1, sd)) 
        result <- log(result) - 
            dnorm(star.M, initial.M, sigma.M, log=T) - 
            dnorm(star.R, initial.R, sigma.R, log=T) 
            #dnorm(star.M, initial.M, sigma.M, log=T) - 
            #dnorm(star.R, initial.R, sigma.R, log=T) 
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
        control=list(trace=999, parscale=parscale, 
        maxit=512))#, abstol=sum(model.list[[1]]$f1.spl(rs)/1000)))
    
    print(best_params)
    
    inversion_params <- 10**best_params$par
    cross.term <- inversion_params[1]
    error.sup <- inversion_params[2]
    width <- inversion_params[3]
    star.M <- inversion_params[4]
    star.R <- inversion_params[5]
    inv.list <- parallelMap(function(model_i) {
        model <- model.list[[model_i]]
        invert.OLA(model=model, rs=rs, 
            cross.term=cross.term, error.sup=error.sup, width=width, 
            targ.kern.type=targ.kern.type, 
            get_cross_kerns=T, get_avg_kerns=T, 
            subtract.mean=F, dM=model$M-star.M, dR=model$R-star.R, 
            kern.interp=seq(0, 1, 0.01),
            perturb=F, num_realizations=1, F_surf=model$F_surf)
    }, model_i=1:length(model.list))
    names(inv.list) <- model.names
    
    inv.lists[[trial_i]] <- Map(function(inv) inv$result, inv=inv.list)
    avg.kerns.lists[[trial_i]] <- Map(function(inv) inv$avg_kern, inv=inv.list)
    cross.lists[[trial_i]] <- Map(function(inv) inv$cross_kern, inv=inv.list)
    star.Ms <- c(star.Ms, star.M)
    star.Rs <- c(star.Rs, star.R)
    
    save(ref.mod, file=paste0("save/ref.mod", targ.mode))
    save(inv.lists, file=paste0("save/inv.lists", targ.mode))
    save(avg.kerns.lists, file=paste0("save/avg.kerns.lists", targ.mode))
    save(cross.lists, file=paste0("save/cross.lists", targ.mode))
    MRs <- list(Ms=star.Ms, Rs=star.Rs)
    save(MRs, file=paste0("save/MRs", targ.mode))
    
    ref.mod2 <- if (target.name=='CygA') get_model(freqs=freqs, 
        model.name='CygAwball', 
        target.name='CygA', 
        k.pair=k.pair, square.Ks=F) else ref.mod 
    make_plots(plot_inversion_lists_mean,
        filename=paste0("inv-lists-optim", targ.mode), 
        model=ref.mod2, inversion.lists=inv.lists, legend.spot=NULL,
        xlim=c(0, 0.45), 
        k.pair=k.pair, sampler=sampler)

}

#save(model.lists, file="model.lists")

#sampler <- c(F, T)#c(F, T, F, T, F, T, F, T, F)

ref.mod2 <- if (target.name=='CygA') get_model(freqs=freqs, 
        model.name='CygAwball', 
        target.name='CygA', 
        k.pair=k.pair, square.Ks=F) else ref.mod 

#sampler <- c(F, F, T, F, F, T, F, T, F, T, F, F) #
#sampler <- c(F, F, F, F, F, T, T, T, F, F, F, F)#c(F, T)
sampler <- T

#plot_inversion_lists_mean(model=ref.mod, inv.lists=inv.lists);dev.off;
make_plots(plot_inversion_lists_mean,
    filename=paste0("inv-lists-optim", targ.mode), 
    model=ref.mod2, inversion.lists=inv.lists, legend.spot=NULL,
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
    kernel.lists=cross.lists, cross=T, legend.spot="topleft",
    xlim=c(0, 0.45),
    ylim=c(-0.2, 1),
    model=ref.mod, k.pair=k.pair, sampler=sampler)#c(T))#

#for (ref.mod in model.list) plot_inversion_lists_mean(ref.mod=ref.mod, 
#    inv.lists=inv.lists)
#plot_inv_diffs_mean(model.list, inv.list)
#dev.off()

#make_plots(plot_inv_diffs_mean,
#    filename=paste0("inv-diffs-mean-optim", targ.mode), 
#    model.list=model.list, inv.list=inv.list)

#make_plots(plot_inversion_mean,
#    filename=paste0("inv-mean-optim", targ.mode), 
#    model=model.list[[5]], inv.list=inv.list)

if (F) {

    load(paste0("save/ref.mod", targ.mode))
    load(paste0("save/inv.lists", targ.mode))
    load(paste0("save/avg.kerns.lists", targ.mode))
    load(paste0("save/cross.lists", targ.mode))
    load(paste0("save/MRs", targ.mode))

}
