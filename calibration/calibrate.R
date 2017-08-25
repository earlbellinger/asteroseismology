#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=0)

### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)
constraint.names = c("log L", "log R", "Fe/H")

### HELPERS
get_Z <- function(Y) {
    # log10(Z/X/0.02293) = 0
    # X+Y+Z = 1 hence...
    X <- -100000*(Y-1)/102293
    1-X-Y
}

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
diffusion <- if (length(args)>0) 1 else 0

### INITIALIZE PARAMETERS FOR SEARCH
if (diffusion) {
    directory <- file.path('sun', 'diffusion') 
    param_names <- c("Y", "alpha", "overshoot", "Z")
    param_init <- c(0.270303087498135, 1.8452097780554, 0.276191587521939, 
        0.0183545123346819)
} else {
    directory <- file.path('sun', 'no_diffusion')
    param_names <- c("Y", "alpha", "overshoot")
    param_init <- c(0.263608790582373, 1.66086916009315, 0.433637018173722)
}
#dir.create(directory, showWarnings = FALSE)
print(directory)

### DEFINE OBJECTIVE FUNCTION 
objective <- function() { 
    ## minimize sum(log(model values / solar values)**2) 
    # searches in LOGS_MS subdirectory of the global 'directory' variable 
    
    hstry_file <- file.path(directory, 'LOGS_MS', 'history.data')
    if (!file.exists(hstry_file)) return(Inf)
    hstry <- read.table(hstry_file, header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    # [Fe/H] = log10 ( Z / X / (Z/X)_Sun )
    mdl_Fe_H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10(Z_div_X_solar)
    mdl_vals <- c(mdl$log_L, mdl$log_R, mdl_Fe_H) 
    
    cat("*** Model values\n")
    cat(paste(constraint.names, mdl_vals, "\n"))
    
    sum(abs(mdl_vals))
}

### SEARCH
iteration <<- 0
ran <<- 0
best_val <<- Inf
best_param <<- param_init

run <- function(params) {
    cat(paste(param_names, params, "\n"))
    
    Y <- params[1]
    alpha <- params[2]
    overshoot <- params[3]
    Z <- if (diffusion) params[4] else get_Z(Y)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.2463  || Y > 0.33  || 
        Z < 0.016   || Z > 0.024 || 
        overshoot <= 0 || overshoot > 0.5 || 
        alpha < 1 || alpha > 3) {
        
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh", 
        '-Y', Y, 
        '-a', alpha, 
        '-o', overshoot, 
        '-Z', Z, 
        '-D', diffusion, 
        '-d', directory)
    print(command)
    system(command)#, ignore.stdout=T)
    
    obj_val <- objective() 
    cat(paste("**** objective value =", obj_val, "\n")) 
    
    if (obj_val < best_val) {
        best_val <<- obj_val
        cat(paste(c("*****", param_names, params, "\n")))
        best_param <<- params
        cat("***** New record!\n")
    }
    
    obj_val
}

result <- optim(par=param_init, fn=run, 
    control=list(trace=999, abstol=1e-8, maxit=10000,
                 reltol=sqrt(.Machine$double.eps)/2))

cat("Optimization terminated. Saving best result\n")
Y <- result$par[1]
alpha <- result$par[2]
overshoot <- result$par[3]
Z <- if (diffusion) result$par[4] else get_Z(Y)
command <- paste("./dispatch.sh", 
    '-Y', Y, 
    '-a', alpha, 
    '-o', overshoot, 
    '-Z', Z, 
    '-D', diffusion, 
    '-d', directory,
    '-s') 
print(command) 
system(command) 
if (result$convergence == 0) cat("Optimization successful.\n") 

