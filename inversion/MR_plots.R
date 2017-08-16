#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)

### LIBRARIES 
source('../scripts/utils.R') 

args <- commandArgs(TRUE)
MOLA            <- if (length(args)>0) as.logical(as.numeric(args[1])) else F
mode.set        <- if (length(args)>1)   args[2] else 'CygA'
error.set       <- if (length(args)>2)   args[3] else 'CygA'
target.name     <- if (length(args)>3)   args[4] else 'modmix'
targ.kern.type  <- if (length(args)>4)   args[5] else 'mod_Gauss'
n_trials        <- if (length(args)>5)   as.numeric(args[6]) else 128
model.list.name <- if (length(args)>6)   args[7] else 'perturbed.model.names'
initial.M       <- if (length(args)>7)   as.numeric(args[8]) else 1
initial.R       <- if (length(args)>8)   as.numeric(args[9]) else 1
sigma.M         <- if (length(args)>9)  as.numeric(args[10]) else 0.016
sigma.R         <- if (length(args)>10) as.numeric(args[11]) else 0.02
#'perturbed.model.names'

targ.mode <- paste0('-p_', target.name, '-m_', mode.set, '-e_', error.set,
    if (MOLA) "-MOLA" else paste0("-", targ.kern.type))

load(paste0("save/MRs", targ.mode))

M <- MRs$Ms
R <- MRs$Rs



