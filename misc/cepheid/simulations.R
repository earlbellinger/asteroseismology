models <- read.table(file.path('data', 'galactic_models.dat'), 
    header=1, stringsAsFactors=F)

I <- models[models$band == 'I',]
V <- models[models$band == 'V',]
DF <- merge(V, I, 
      by=c('Z', 'Y', 'Mass', 'Type', 'Temp', 'Luminosity', 'Radius', 'logP'), 
      suffixes=c('.V', '.I'))

write.table(DF, 'simulations.dat', quote=FALSE, sep='\t', row.names=FALSE)
