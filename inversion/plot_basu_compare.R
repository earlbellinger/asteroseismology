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


#star.names <- c(8006161, 6106415, 5774694, 5184732, 
#    7680114, 8379927, 8694723, 8938364, 
#    9098294, 9139163, 10454113, 10516096, ) 
star.names <- list.files('save')
star.names <- star.names[grepl('avg', star.names)]
star.names <- unname(sapply(star.names, function(x) 
    strsplit(strsplit(x, 'KIC_')[[1]][2], '-')[[1]][1]))
#initial.Ms <- c(1.014, 1.118, 1.000, 1.271) 
#sigma.Ms <- c(0.031, 0.052, 0.0357, 0.047) 
basu.index <- read.table(file.path('models', 'basu_legacy', 'index'),
    header=1)

redo <- c('8379927', '8694723', '9139163', '10454113')

for (ii in 1:length(star.names)) {#length(redo)) {#
    
    star.name <- star.names[ii]#redo[ii]#
    #initial.M <- initial.Ms[ii]#1.014
    #sigma.M <- sigma.Ms[ii]#0.031
    target.name <- paste0('KIC_', star.name)
    mode.set <- target.name
    error.set <- target.name
    ref.mod.name <- paste0(star.name, '_meanRmeanM')
    
    targ.mode <- paste0(
        '-p_', target.name, 
        '-m_', mode.set, 
        '-e_', error.set, 
        '-r_', ref.mod.name, 
        paste0("-", targ.kern.type))
    
    #model.names <- names.list[[star.name]]
    #model.names <- get(model.list.name) 
    #model.list <- parallelMap(function(model.name) 
    #        get_model(freqs=NULL, model.name=model.name, 
    #            target.name=target.name, k.pair=k.pair, square.Ks=F), 
    #    model.name=model.names)
    #names(model.list) <- model.names
    
    load(file.path('save', paste0('ref.mod',         targ.mode)))
    load(file.path('save', paste0('inv.lists',       targ.mode)))
    load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
    load(file.path('save', paste0('cross.lists',     targ.mode)))
    load(file.path('save', paste0('MRs',             targ.mode)))
    load(file.path('save', paste0('inv.params',      targ.mode)))
    
    index <- basu.index$KIC == star.name
    if (!any(index)) next 
    model.number <- basu.index$Model[index]
    
    nod.model <- file.path('models', 'basu_legacy', paste0(model.number, '.struc'))
    dif.model <- file.path('models', 'basu_legacy', paste0(model.number, '_dif.struc'))
    
    rm.first.kernel <- star.name %in% c('8379927', '8694723', '9139163', '10454113')
    rm.second.kernel <- star.name %in% c('10454113')
    rm.fifth.kernel <- star.name %in% c('10454113', '8694723')
    rm.last.kernel <- star.name %in% c('8379927', '8694723', '10454113')
    
    for (mod_i in 1:2) {
        
        ref.mod.name <- c(nod.model, dif.model)[mod_i]
        difno <- c('dif', 'nod')[mod_i]
        
        if (!file.exists(ref.mod.name)) next 
        
        #ref.mod.name <- paste0(star, 'Basu')
        model.prof <- read.table(ref.mod.name)[,c(2, 5, 6, 9)]
        names(model.prof) <- c('r', 'rho', 'P', 'eps')
        model <- list(k.pair=k.pair, 
            f1.spl=with(model.prof, splinefun(r, P/rho)))
        
        core.bound <- max(model.prof$r[model.prof$eps>0.01*max(model.prof$eps)])
        
        con <- file(ref.mod.name, "r")
        line <- readLines(con, n=2)
        close(con)
        mass <- as.numeric(strsplit(strsplit(line, " Msun")[[2]][1], "gm \\(")[[1]][2])
        
        inversion <- lists_to_inversion(model=model, rs=rs, 
            inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
            cross.lists=cross.lists, inv.params=inv.params, 
            kern.interp.xs=kern.interp.xs) 
        
        targ_kerns <- get_target_kernels(inversion$result$rs, 
            inversion$params[['widths']], 
            kern.interp.xs, r_f=0.2, 
            f.spl=ref.mod$cs.spl)
        
        penalties <- sapply(1:nrow(inversion$result), function(jj) {
            avg_kern <- inversion$avg_kerns[,jj]
            targ_kern <- targ_kerns[,jj]
            sintegral(kern.interp.xs, (avg_kern-targ_kern)**2)$value
        })
        
        penalties2 <- sapply(1:nrow(inversion$result), function(jj) {
            avg_kern <- inversion$avg_kerns[,jj]
            #targ_kern <- targ_kerns[,ii]
            mid <- inversion$result[jj,]$fwhm.mid
            FWHM <- with(inversion$result[jj,], fwhm.right - fwhm.left)
            a <- max(mid - FWHM, 0)
            b <- min(mid + FWHM, 1)
            first <- if (a>0.001) sintegral(kern.interp.xs[kern.interp.xs<a], 
                avg_kern[kern.interp.xs<a]**2)$value else 0
            second <- if (b<0.999) sintegral(kern.interp.xs[kern.interp.xs>b], 
                avg_kern[kern.interp.xs>b]**2)$value else 0
            first + second 
        })
        
        surf_term <- sapply(1:nrow(inversion$result), function(jj) {
            avg_kern <- inversion$avg_kerns[,jj]
            max(abs(avg_kern[kern.interp.xs > 0.9]))
        })
        
        sampler <- surf_term < 5
        if (rm.first.kernel) sampler[1] <- F
        if (rm.second.kernel) sampler[2] <- F
        if (rm.fifth.kernel) sampler[5] <- F
        if (rm.last.kernel) sampler[length(sampler)] <- F 
        
        make_plots_inversion_all(model, inversion, kern.interp.xs=kern.interp.xs,
            k.str=paste0("-Basu-", difno, '-', star.name), cross.inset="bottomright",
            inversion_ylim = c(-0.3, 0.3),
            caption=star.name, 
            caption2=as.expression(bquote('M' == .(signif(mass,3)))),
            sampler=sampler,
            core.bound=core.bound, 
            #sampler=penalties<4 & penalties2 > 1,
            cross_kern_ylim=c(-0.8, 0.3))
        
    }
}



