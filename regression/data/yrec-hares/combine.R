#### Combine YREC models into one big file 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Stellar Ages & Galactic Evolution Group
#### Max-Planck-Institut fur Sonnensystemforschung

source(file.path('/scratch', 'seismo', 'bellinger', 
    'asteroseismology', 'scripts', 'utils.R'))

fname <- 'yrec_hares.dat'

files <- list.files('.')
files <- files[grepl('full_', files)]
#if (fname %in% files) files <- files[-which(files == fname)]
dnus <- files[grepl('dnu', files)]
struc <- files[!grepl('dnu', files)]

DF <- do.call(rbind, Map(function(ii) {
    struc.i <- read.table(struc[ii], header=1)
    dnus.i <- read.table(dnus[ii], col.names=c("Dnu0"))
    cbind(struc.i, dnus.i)
}, ii=1:length(struc)))

nu_max_classic <- with(DF, scaling_nu_max(radius, M, Teff))
nu_max <- with(DF, scaling_nu_max_Viani(radius, M, Teff, mu))

DF <- cbind(DF, nu_max_classic=nu_max_classic, nu_max=nu_max)

FeH <- with(DF, data.frame( log10( Z_surf/X_surf/0.02293 ) ))
colnames(FeH) = "Fe/H"
DF <- cbind(DF, FeH)

write.table(DF, fname, col.names=T, row.names=F, quote=F, sep='\t')

