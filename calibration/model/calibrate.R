#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

source('../scripts/seismology.R')

### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)
constraints <- list("log L" = c(log10(1.56), 0.05/1.56/2.3),
                    "log R" = c(log10(1.22), 0.02/1.22/2.3),
                    "Fe/H"  = c(0.096, 0.026),
                    "log T" = c(log10(5825), 50/5825/2.3),
                    "log g" = c(4.33, 0.07))

seismic.data <- read.table(file.path('..', 'inverse', 'perturb', 
    '16CygA_perturb.dat'), header=1)
seismic.means <- apply(seismic.data, 2, mean)
seismic.stds <- apply(seismic.data, 2, sd)
constraints[['r02']] <- c(log10(seismic.means[['r02']]), 
    seismic.stds[['r02']]/seismic.means[['r02']]/2.3)
constraints[['r10']] <- c(log10(seismic.means[['r10']]), 
    seismic.stds[['r10']]/seismic.means[['r10']]/2.3)

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
diffusion <- if (length(args)>0) 1 else 0
directory <- if (diffusion) "diffusion" else "no_diffusion"
print(directory)

### INITIALIZE PARAMETERS FOR SEARCH
param_names <- c("M", "age", "alpha", "Y", "Z")
param_init <- log10(c(1.08, 6.9, 1.86, 0.262, 0.022))

### DEFINE OBJECTIVE FUNCTION
chisq <- function() {
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    freqs <- read.table(file.path(directory, 'LOGS_MS', 'profile1-freqs.dat'), 
        col.names=c('l', 'n', 'nu', 'E'))
    seis.DF <- seismology(freqs, mdl$nu_max)
    if (!'r02' %in% colnames(seis.DF) || !'r10' %in% colnames(seis.DF)) 
        return(Inf)
    log_r02 <- log10(seis.DF[['r02']])
    log_r10 <- log10(seis.DF[['r10']])
    
    mdl.Fe.H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10_Z_div_X_solar
    mdl.vals <- c(
        (mdl$log_L    - constraints[['log L']][1]) / constraints[['log L']][2], 
        (mdl$log_R    - constraints[['log R']][1]) / constraints[['log R']][2], 
        (mdl.Fe.H     - constraints[['Fe/H']][1])  / constraints[['Fe/H']][2], 
        (mdl$log_Teff - constraints[['log T']][1]) / constraints[['log T']][2],
        (mdl$log_g    - constraints[['log g']][1]) / constraints[['log g']][2],
        (log_r02      - constraints[['r02']][1])   / constraints[['r02']][2],
        (log_r10      - constraints[['r10']][1])   / constraints[['r10']][2])
        #log10(mdl$star_age/10**9/4.572))
    
    chi2 <- mdl.vals**2
    
    cat("** Model values and chi^2\n")
    cat(paste(names(constraints), mdl.vals, chi2, "\n"))
    
    chi2
}

### SEARCH
iteration <<- 0
ran <<- 0
best_chisq <<- Inf
best_param <<- param_init

run <- function(params) {
    cat(paste(param_names, params, "\n"))
    
    params <- 10**params 
    M <- round(params[1], 8)
    age <- round(params[2], 8)
    alpha <- round(params[3], 8)
    Y <- round(params[4], 8)
    Z <- round(params[5], 8) 
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.2463  || Y > 0.33  || 
        Z < 0.015   || Z > 0.03 || 
        alpha < 1.5 || alpha > 2.5 || 
        M < 1       || M > 1.3 || 
        age < 4     || age > 9) {
        
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh",
        '-M', M, 
        '-Y', Y, 
        '-a', alpha, 
        '-Z', Z, 
        '-D', diffusion, 
        '-t', paste0(age, 'e9'), 
        '-d', directory)
    print(command)
    system(command)#, ignore.stdout=T)
    
    chi_sq <- sum(chisq())
    cat(paste('**** chi_sq =', chi_sq, "\n"))
    
    if (chi_sq < best_chisq) {
        best_chisq <<- chi_sq
        cat(paste(c("*****", param_names, params, "\n")))
        best_param <<- params
        cat("***** New record!\n")
    }
    
    chi_sq
}

result <- optim(par=param_init, fn=run, control=list(trace=999))

cat("Optimization terminated. Saving best result\n")
params <- 10**result$par 
M <- round(params[1], 8)
age <- round(params[2], 8)
alpha <- round(params[3], 8)
Y <- round(params[4], 8)
Z <- round(params[5], 8) 
command <- paste("./dispatch.sh",
    '-M', M, 
    '-Y', Y, 
    '-a', alpha, 
    '-Z', Z, 
    '-D', diffusion, 
    '-t', paste0(age, 'e9'), 
    '-d', directory) 
print(command) 
system(command) 
print(result) 
if (result$convergence == 0) cat("Optimization successful.\n") 


