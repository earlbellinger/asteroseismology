#### Combine all modelled frequencies of the Sun 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

fgongs.dir <- 'high_ell'

combine <- function(filename, directories) {
    print(filename)
    directories <- directories[file.exists(file.path(directories, filename))]
    ker <- Reduce(function(...) merge(..., all=T, by='x'), Map(function(ell) {
        if (file.info(file.path(ell, filename))$size <= 0) 
            return(data.frame(x=0))
        DF <- read.table(file.path(ell, filename), header=1)
        if ("l..1_n..1" %in% names(DF)) {
            DF <- DF[,-which(names(DF) == "l..1_n..1")]
        }
        DF
    }, ell=directories))
    
    output.dir <- dirname(ells[1])
    print(paste('Writing', file.path(output.dir, filename)))
    write.table(ker, file.path(output.dir, filename), quote=F, row.names=F)
}

for (model in c('diffusion', 'no_diffusion')) {
    model.dir <- file.path(fgongs.dir, model)
    ells <- list.dirs(path = model.dir, full.names = TRUE, recursive = TRUE)
    ells <- ells[!grepl('maybe_sub_logs', ells)][-1]
    ells <- ells[!grepl('psi', ells)]#[-1]
    
    combine('E_K_u-Y.dat', ells)
    combine('E_K_Y-u.dat', ells)
    
    combine('E_K_c2-Gamma1.dat', ells)
    combine('E_K_Gamma1-c2.dat', ells)
    
    combine('E_K_c2-rho.dat', ells)
    combine('E_K_rho-c2.dat', ells)
    
    combine('E_K_rho-Gamma1.dat', ells)
    combine('E_K_Gamma1-rho.dat', ells)
    
    combine('E_K_u-Gamma1.dat', ells)
    combine('E_K_Gamma1-u.dat', ells)
    
    combine('E_K_rho-Y.dat', ells)
    combine('E_K_Y-rho.dat', ells)
    
    freqs <- do.call(rbind, Map(function(ell) {
            read.table(file.path(ell, paste0(model, '-freqs.dat')),
                col.names=c('l', 'n', 'nu', 'E'))
        }, ells))
    rownames(freqs) <- NULL
    
    write.table(unique(freqs[order(freqs$l, freqs$n),]), 
        file.path(model.dir, paste0(model, '-freqs.dat')),
        quote=F, row.names=F)
    
}

