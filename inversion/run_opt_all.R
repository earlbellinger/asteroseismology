#### Runs OLA_optimize_MR2 for all stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)
options(scipen=99999)

perturb.dir <- file.path('..', 'regression', 'perturb', 'inversions')
covs.dir <- file.path('..', 'regression', 'learn-inv', 'covs-inv', 'inversions')
cov.files <- list.files(covs.dir)

for (cov.file in cov.files) {
    star <- strtoi(strsplit(cov.file, '.dat')[[1]][1])
    if (is.na(star)) next 
    
    DF <- read.table(file.path(covs.dir, cov.file), header=1)
    
    if (!'radius' %in% names(DF)) {
        perturb.file <- file.path(perturb.dir, paste0(star, '_perturb.dat'))
        if (!file.exists(perturb.file)) next 
        perturb <- read.table(perturb.file, header=1)
        DF <- cbind(perturb, DF)
    }
    
    M.mean <- mean(DF$M) 
    R.mean <- mean(DF$radius) 
    cov.M <- cov(data.frame(M=DF$M, radius=DF$radius))
    
    command <- paste0('maybe_sub.sh -e -m 6000000 -p 9', 
        ' Rscript OLA_optimize_MR2.R',
        ' KIC_', star, 
        ' KIC_', star, 
        ' KIC_', star, 
        ' ', star, '_meanRmeanM', 
        ' ', star,
        ' mod_Gauss ',
        M.mean, ' ',
        R.mean, ' ',
        cov.M[[1]], ' ',
        cov.M[[4]], ' ',
        '128 F ',
        cov.M[[2]], ' ',
        cov.M[[3]])
    print(command) 
    #system(command) 
}

