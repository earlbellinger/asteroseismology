#### Collect evolutionary models 
#### Author: Earl Patrick Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus

#library(parallel)
#library(parallelMap)

args <- commandArgs(TRUE)
sim_dir <- if (length(args)>0) args[1] else 'simulations'
output_fname <- paste0(sim_dir, '.dat')

simulations <- file.path(sim_dir, list.files(sim_dir))
simulations <- simulations[grep('.dat', simulations)]

load_data <- function(filename) read.table(filename, header=1, check.names=0)

# Load and combine data over all simulations 
#parallelStartMulticore(max(1, min(detectCores(), 62)))
seis.DF <- do.call(plyr:::rbind.fill, Map(load_data, simulations))
print(paste(nrow(seis.DF), "rows"))
print(sapply(seis.DF, fivenum))

# Save data
write.table(seis.DF, output_fname, quote=FALSE, sep='\t', row.names=FALSE)
