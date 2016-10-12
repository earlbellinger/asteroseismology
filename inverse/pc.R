#### Total least squares via ``principal components'' regression
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(deming)

col.pal <- c('black', 'black', blue, '#900090', red)

data_dir <- file.path('learn', 'covs-simulations')
data_dir.1 <- file.path(data_dir, 'perturb')
data_dir.2 <- file.path(data_dir, 'kages')
data_dir.3 <- file.path(data_dir, 'legacy')

ml <- do.call(plyr:::rbind.fill, 
    Map(function(covs) {
        name <- sub('.dat', '', basename(covs))
        if (grepl("amp", name) || name=="Tagesstern") return(NULL)
        DF <- read.table(covs, header=1)
        results <- data.frame(Name=name)
        for (column in names(DF)) {
            vals <- DF[,column]
            result <- data.frame(mean(vals), sqrt(var(vals)), 
                                 median(vals) - quantile(vals, .16), 
                                 quantile(vals, .84) - median(vals))
            colnames(result) <- c(column, paste0(column, '_s'), 
                                  paste0('d', column, 'L'),
                                  paste0('d', column, 'H'))
            results <- cbind(results, result)
        }
        results
    }, c(file.path(data_dir.1, list.files(data_dir.1)),
         file.path(data_dir.2, list.files(data_dir.2)),
         file.path(data_dir.3, list.files(data_dir.3)))
    )
)
ml <- ml[!duplicated(ml$Name),]

