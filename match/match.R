
#args <- commandArgs(TRUE)
source(file.path('..', 'scripts', 'seismology.R'))


star <- read.table(file.path('..', 'regression', 'perturb', 'inversions', 
    '10068307_perturb.dat', header=1))

