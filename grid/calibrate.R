#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=100000)
library(optimx)

### CONSTANTS 
Z_div_X_solar = 0.02293 #0.0245 # GN93 #0.02293 # GS98 
log10_Z_div_X_solar = log10(Z_div_X_solar) 
constraint.names = c("log L", "log R", "Fe/H")#, "bcz", "Dnu") 

args <- commandArgs(TRUE)
diffusion <- length(args) > 0

### INITIALIZE PARAMETERS FOR SEARCH 
directory <- file.path('calibrate') 
subdir <- ifelse(diffusion, 'diffusion', 'no_diffusion')
param_names <- c("Y", "alpha", "Z") 
#param_init <- c(0.273032382308896, 1.81835067664561, 0.0811842658081544, 
#    0.0196326242552884, 0.00585259361551818) 

param_init <- c(0.272804509755598, 1.83454562644668, 0.0185655106719405) 
if (!diffusion) param_init <- c(0.264382485315751, 1.68697171300493)

#print(directory)

get_Z <- function(Y) {
    # log10(Z/X/0.02293) = 0
    # X+Y+Z = 1 hence...
    X <- -100000/102293*(Y-1)
    1-X-Y
}

### DEFINE OBJECTIVE FUNCTION 
objective <- function() { 
    ## minimize sum(log(model values / solar values)**2) 
    # searches in LOGS_MS subdirectory of the global 'directory' variable 
    
    hstry_file <- file.path(directory, subdir, 'LOGS_3MS', 'history.data')
    if (!file.exists(hstry_file)) return(Inf)
    hstry <- read.table(hstry_file, header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    # [Fe/H] = log10 ( Z / X / (Z/X)_Sun )
    mdl_Fe_H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10(Z_div_X_solar) 
    #mdl_bcz <- log10( mdl$cz_bot_radius / 0.713 ) 
    #mdl_Dnu <- log10( mdl$Delta_nu / 134.863 )
    mdl_vals <- c(mdl$log_L, mdl$log_R, mdl_Fe_H) #, mdl_bcz, mdl_Dnu) 
    
    cat("*** Model values\n") 
    cat(paste(constraint.names, mdl_vals, "\n")) 
    cat(paste(c('L', 10**mdl$log_L, 'R', 10**mdl$log_R, "\n")))
    
    result <- sum(mdl_vals**2) 
    if (is.finite(result)) result else 10**10
}

### SEARCH
iteration <<- 0
ran <<- 0
best_val <<- Inf
best_param <<- param_init

run <- function(params) {
    Y         <- params[1]
    alpha     <- params[2]
    #overshoot <- params[3]
    Z         <- ifelse(diffusion, params[3], get_Z(Y))
    #f0        <- params[5]
    
    if (Y < 0.26 || Y > 0.31 || Z < 0.01 || Z > 0.03 || alpha < 1 || alpha > 3)
        return(10)
    
    #if (Y < 0.2463) Y <- 0.2463
    #if (Y > 0.33)   Y <- 0.33
    #if (Z < 0.01 && diffusion)   Z <- 0.01
    #if (Z > 0.03 && diffusion)   Z <- 0.03
    #if (overshoot < 0) overshoot <- 0
    #if (overshoot > 1) overshoot <- 1
    #if (alpha < 1) alpha <- 0.5
    #if (alpha > 3) alpha <- 3
    #if (f0 <= 0) f0 <- 10**-6
    #if (f0 > overshoot) f0 <- overshoot
    
    #cat(paste(param_names, c(Y, alpha, overshoot, Z, f0), '\n'))
    cat(paste(param_names, c(Y, alpha, Z), '\n'))
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh", 
        '-Y', Y, 
        '-a', alpha, 
        #'-o', overshoot, 
        #'-f', f0, 
        '-Z', Z, 
        '-D', ifelse(diffusion, 1, 0), 
        '-g', ifelse(diffusion, 1, 0),
        '-c', 4572000000,
        '-d', directory,
        '-n', subdir,
        '-p')
    print(command)
    system(command)
    
    obj_val <- log10(objective()) 
    cat(paste("**** objective value =", obj_val, "\n")) 
    
    if (obj_val < best_val) {
        best_val <<- obj_val
        cat(paste(c("*****", param_names, params, "\n")))
        best_param <<- params
        cat("***** New record!\n")
    }
    
    obj_val
}

#result <- optim(par=param_init, fn=run, method="BFGS", 
#    control=list(trace=999, abstol=1e-8, maxit=10000,
#                 reltol=sqrt(.Machine$double.eps)/2))
result <- optimx(par=param_init, fn=run, method="Nelder-Mead")

cat("Optimization terminated. Saving best result\n")
Y         <- result$p1#ar[1]
alpha     <- result$p2#ar[2]
#overshoot <- result$p3#ar[3]
Z         <- ifelse(diffusion, result$p3, get_Z(Y))
#f0        <- result$p5#ar[5]
command   <- paste("./dispatch.sh", 
    '-Y', Y, 
    '-a', alpha, 
    #'-o', overshoot, 
    #'-f', f0, 
    '-Z', Z, 
    '-D', ifelse(diffusion, 1, 0), 
    '-g', ifelse(diffusion, 1, 0),
    '-c', 4572000000, 
    '-d', directory,
    '-n', subdir,
    '-f') 
print(command) 
system(command) 
print(result)
#if (result$convergence == 0) cat("Optimization successful.\n") 

