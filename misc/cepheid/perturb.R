#### Generate Monte-Carlo perturbations of observed stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(parallel)
library(parallelMap)

out_dir <- file.path('perturb')
dir.create(out_dir, showWarnings=FALSE)

I <- read.table(file.path('data', 'galactic_I.dat'), header=1, 
    stringsAsFactors=F)
V <- read.table(file.path('data', 'galactic_V.dat'), header=1, 
    stringsAsFactors=F)
stars <- V$ID[V$ID %in% I$ID]

n_perturbations = 1000

### Obtain properties of real stars varied within their uncertainties 
perturb <- function(star, n_perturbations) {
    V_band <- V[V$ID == star,]
    I_band <- I[I$ID == star,]
    
    merged <- merge(V_band, I_band, by=c('ID', 'logP'), suffixes=c('.V', '.I'))
    err <- which(grepl('err', names(merged)))
    DF <- merged[,-err]
    DF.err <- merged[,err]
    
    perturb.names <- sub('_err', '', names(DF.err))
    
    num.qtys <- ncol(DF.err)
    qtys <- as.numeric(DF[perturb.names])
    errs <- as.numeric(DF.err)
    phis <- grepl('phi', perturb.names)
    
    res <- data.frame(do.call(rbind, Map(function(n) {
            perturbed <- rnorm(num.qtys, qtys, if (n==1) 0 else errs)
            perturbed[phis] <- perturbed[phis] %% (2*pi)
            perturbed
        }, 
        1:n_perturbations)))
    colnames(res) <- perturb.names
    
    merge(with(merged, data.frame(logP, Sk.V, Ac.V, Sk.I, Ac.I)), res)
}

process_star <- function(star, out_dir=file.path('perturb')) {
    print(paste("Processing", star))
    result <- perturb(star, n_perturbations)
    if (!is.null(result)) {
        write.table(result, 
            file.path(out_dir, paste0(star, ".dat")), 
            quote=FALSE, sep='\t', row.names=FALSE)
    }
}

## Perturb every star 10k times and save the results
parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
parallelMap(process_star, stars)
