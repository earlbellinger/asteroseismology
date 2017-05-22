library(magicaxis)

data.dir <- file.path('perturb', 'kages')
all.files <- file.path(data.dir, list.files(data.dir))
all.files <- c(all.files, file.path('perturb', 
    c('Sun_perturb.dat', '16CygA_perturb.dat', '16CygB_perturb.dat')))

separations <- function(filename) {
    DF <- read.table(filename, header=1)
    data.frame(Teff=mean(DF$Teff), Teff_s=sqrt(var(DF$Teff)),
               Dnu=mean(DF$Dnu0), Dnu_s=sqrt(var(DF$Dnu0)),
               dnu=mean(DF$dnu02), dnu_s=sqrt(var(DF$dnu02)))
}

DF <- do.call(rbind, Map(separations, all.files))

