#### Solar model calibration 
#### Find the Y_0, Z_0, and alpha that match log L, log R, [Fe/H] and log age 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

library(optimx)
options(scipen=99999)

### PARSE COMMAND LINE ARGS
args  <- commandArgs(TRUE)
star  <- if (length(args)>0)             args[1]  else 'Sun'
lowh  <- if (length(args)>1)             args[2]  else 'meanRmeanM'
mass  <- if (length(args)>2)  as.numeric(args[3]) else 1
log_R <- if (length(args)>3)  as.numeric(args[4]) else 0
age   <- if (length(args)>4)  as.numeric(args[5]) else 4.572*10**9
log_L <- if (length(args)>5)  as.numeric(args[6]) else 0
Fe_H  <- if (length(args)>6)  as.numeric(args[7]) else 0
Y_0   <- if (length(args)>7)  as.numeric(args[8]) else 0.27
Z_0   <- if (length(args)>8)  as.numeric(args[9]) else 0.02
alpha <- if (length(args)>9) as.numeric(args[10]) else 1.85
directory <- file.path('models', star, lowh)
dir.create(directory)#, showWarnings = FALSE)

### CONSTANTS 
Z_div_X_solar = 0.02293 # Grevesse & Sauval 1998

get_Z <- function(Y, Fe_H=0) {
    # [Fe/H] = log10(Z/X/0.02293)
    # X+Y+Z = 1 hence...
    (1-Y) / ( 1 + 10**( -Fe_H - log10(Z_div_X_solar) ) )
}

### INITIALIZE PARAMETERS FOR SEARCH 
constraint_names <- c("log L", "log R", "[Fe/H]", "log age") # solar units 
param_names <- c("Y_0", "alpha") # parameters to optimize 
param_init <- c(Y_0, alpha) # initial parameter values 

### DEFINE OBJECTIVE FUNCTION 
objective <- function() { 
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    # [Fe/H] = log10 ( Z / X / (Z/X)_Sun )
    mdl_Fe_H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10(Z_div_X_solar)
    
    mdl_vals <- c(mdl$log_L - log_L, # log (luminosity / solar luminosity) 
                  mdl$log_R - log_R, # same with radius 
                  mdl_Fe_H  - Fe_H,  # surface metallicity 
                  log10(mdl$star_age / age)) 
    
    cat("*** Model values\n") 
    cat(paste(constraint_names, mdl_vals, "\n")) 
    cat(paste(c('L', 10**log_L, 10**mdl$log_L, "\n")))
    cat(paste(c('R', 10**log_R, 10**mdl$log_R, "\n")))
    
    sum(mdl_vals**2)
}

### SEARCH 
precision <<- 8 # number of decimal places considered 
iteration <<- 0 # for bookkeeping: how many parameter sets attempted 
ran <<- 0 # for bookkeeping: how many MESA jobs actually have been run 
best_val <<- Inf 
best_param <<- param_init 

run <- function(params) {
    params <- round(params, precision) # truncate precision of search 
    cat("\n* Trying parameters:\n")
    cat(paste(param_names, params, "\n"))
    
    Y_0 <- round(params[1], precision)
    #Z_0 <- round(params[2], precision)
    alpha <- round(params[2], precision)
    Z_0 <- get_Z(Y_0, Fe_H)
    
    iteration <<- iteration + 1
    cat(paste("** iter:", iteration, "\n"))
    
    if (Y_0 < 0.15) Y_0 <- 0.15
    if (Y_0 > 0.45) Y_0 <- 0.45
    #if (Z_0 < 0.0001) Z_0 <- 0.0001
    #if (Z_0 > 0.1)    Z_0 <- 0.1
    if (alpha < 0.5) alpha <- 0.5
    if (alpha > 4) alpha <- 4
    
    # check if parameters are reasonable 
    #if (Y_0 < 0.2463 || Y_0 > 0.38 || 
    #    Z_0 < 0.001  || Z_0 > 0.1  || 
    #    alpha < 1    || alpha > 3) {
    #    
    #    cat("Input parameters out of bounds; skipping\n")
    #    return(Inf)
    #}
    
    ran <<- ran + 1 
    cat(paste("**  ran:", ran, "\n")) 
    
    # call 'dispatch.sh' shell script to run MESA with supplied parameters 
    command <- paste("./dispatch.sh", 
        '-M', mass, 
        '-Y', Y_0, 
        '-Z', Z_0, 
        '-a', alpha, 
        '-t', age, 
        '-D', 0, 
        '-d', directory) 
    print(command) 
    system(command) 
    
    # calculate objective function for this set of parameters 
    obj_val <- objective() 
    cat(paste("**** objective value =", obj_val, "\n")) 
    
    # save the best value (optional: optim takes care of this for us anyway) 
    if (obj_val < best_val) { 
        best_val <<- obj_val 
        cat(paste("*****", param_names, params, "\n")) 
        best_param <<- params 
        cat("***** New record! \n") 
    } 
    
    # return the objective value to the optimizer 
    log10(obj_val) 
}

cat("Initializing optimization\n") 
cat(paste("Outputting results in directory:", directory, "\n")) 

result <- optimx(par=param_init, fn=run, method="Nelder-Mead")#, 
    #lower=c(0.2463, 0.0001, 1), upper=c(0.38, 0.1, 3))

cat("Optimization terminated. Saving best result\n")
Y_0 <- result$p1
#Z_0 <- result$p2
alpha <- result$p2
Z_0 <- get_Z(Y_0, Fe_H)
command <- paste("./dispatch.sh", 
    '-M', mass, 
    '-Y', Y_0, 
    '-Z', Z_0, 
    '-a', alpha, 
    '-t', age, 
    '-D', 0, 
    '-s', 1, 
    '-d', directory) 
print(command) 
system(command) 
print(result) 
#if (result$convergence == 0) cat("Optimization successful.\n") 
objective()

