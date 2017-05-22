#### Solar model calibration 
#### Find the Y_0, Z_0, and alpha that match log L, log R, [Fe/H] and log age 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=99999)

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
mass <- if (length(args)>0) as.numeric(args[1]) else 1
log_R <- if (length(args)>1) as.numeric(args[2]) else 0
age <- if (length(args)>2) as.numeric(args[3]) else 4.572*10**9
log_L <- if (length(args)>3) as.numeric(args[4]) else 0
Fe_H <- if (length(args)>4) as.numeric(args[5]) else 0
directory <- file.path('calibrate2', 
    paste0('M=', mass, '_logR=', log_R, '_age=', age))
dir.create(directory, showWarnings = FALSE)

### CONSTANTS 
Z_div_X_solar = 0.02293 # Grevesse & Sauval 1998 

### INITIALIZE PARAMETERS FOR SEARCH 
#directory <- "result" # location for MESA to run 
constraint_names <- c("log L", "log R", "[Fe/H]", "log age") # solar units 
param_names <- c("Y_0", "Z_0", "alpha") # parameters to optimize 
param_init <- c(0.27095138, 0.01847179, 1.81945022) # initial parameter values 
#parscale <-   c(0.02,       0.005,      0.25) # for scaling the search 
#abstol <- 10**-7 # stop when the objective value is lower than this amount 

### DEFINE OBJECTIVE FUNCTION 
objective <- function() { 
    ## minimize sum(log(model values / solar values)**2) 
    # searches in LOGS_MS subdirectory of the global 'directory' variable 
    
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
    
    sum(mdl_vals**2)
}

### SEARCH 
precision <<- 12 # number of decimal places considered 
iteration <<- 0 # for bookkeeping: how many parameter sets attempted 
ran <<- 0 # for bookkeeping: how many MESA jobs actually have been run 
best_val <<- Inf 
best_param <<- param_init 

run <- function(params) {
    ## dispatch a stellar evolutionary track with the input parameters 
    # must be list with Y_0, Z_0, and alpha in that order 
    
    params <- round(params, precision) # truncate precision of search 
    cat("\n* Trying parameters:\n")
    cat(paste(param_names, params, "\n"))
    
    Y_0 <- params[1]
    Z_0 <- params[2]
    alpha <- params[3]
    
    iteration <<- iteration + 1
    cat(paste("** iter:", iteration, "\n"))
    
    # check if parameters are reasonable 
    if (Y_0 < 0.2    || Y_0 > 0.35 || 
        Z_0 < 0.001  || Z_0 > 0.04 || 
        alpha < 1    || alpha > 3) {
        
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**  ran:", ran, "\n")) 
    
    # call 'dispatch.sh' shell script to run MESA with supplied parameters 
    command <- paste("./dispatch.sh", 
        '-M', mass,
        '-Y', Y_0, 
        '-Z', Z_0, 
        '-a', alpha,
        '-t', age, 
        '-D', 1,
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
    obj_val 
}

cat("Initializing optimization for solar parameters \n") 
cat(paste("Outputting results in directory:", directory, "\n")) 

result <- optim(par=param_init, fn=run, 
    control=list(trace=999))#, abstol=abstol, parscale=parscale))

cat("Optimization terminated. Saving best result\n")
Y_0 <- result$par[1]
Z_0 <- result$par[2]
alpha <- result$par[3]
command <- paste("./dispatch.sh", 
    '-M', mass, 
    '-Y', Y_0, 
    '-Z', Z_0, 
    '-a', alpha, 
    '-t', age,
    '-D', 1,
    '-d', directory) 
print(command) 
system(command) 
print(result) 
if (result$convergence == 0) cat("Optimization successful.\n") 

