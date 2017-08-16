#### Solar model calibration 
#### Find the Y_0, Z_0, and alpha that match log L, log R, [Fe/H] and log age 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=99999)

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
spectral_type <- if (length(args)>0) args[1] else "G2V"
log_g <- if (length(args)>1) as.numeric(args[2]) else 4.438
Teff <- if (length(args)>2) as.numeric(args[3]) else 5800
M <- if (length(args)>3) as.numeric(args[4]) else 1
tau <- if (length(args)>4) as.numeric(args[5]) else 4.572*10**9
alpha_MLT <- if (length(args)>5) as.numeric(args[6]) else 1.67123747 
alpha_ov <- if (length(args)>6) as.numeric(args[7]) else 0.09104194 
directory <- file.path('models', spectral_type)
dir.create(directory, showWarnings = FALSE)

### CONSTANTS 
Z_div_X_solar = 0.02293 # Grevesse & Sauval 1998 

### INITIALIZE PARAMETERS FOR SEARCH 
#directory <- "result" # location for MESA to run 
constraint_names <- c("log g", "Teff") # solar units 
param_names <- c("mass", "log age", "alpha_MLT", "alpha_ov") 
param_init <- c(M, log10(tau), alpha_MLT, alpha_ov) # initial parameter values 

### DEFINE OBJECTIVE FUNCTION 
objective <- function() { 
    ## minimize sum(log(model values / solar values)**2) 
    # searches in LOGS_MS subdirectory of the global 'directory' variable 
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    # [Fe/H] = log10 ( Z / X / (Z/X)_Sun )
    #mdl_Fe_H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10(Z_div_X_solar)
    
    #mdl_vals <- c(mdl$log_L - log_L, # log (luminosity / solar luminosity) 
    #              mdl$log_R - log_R, # same with radius 
    #              mdl_Fe_H  - Fe_H,  # surface metallicity 
    #              log10(mdl$star_age / age)) 
    mdl_vals <- c(mdl$log_g - log_g,
                  mdl$log_Teff - log10(Teff))
    
    cat("*** Model values\n")
    cat(paste(constraint_names, mdl_vals, "\n"))
    
    sum(abs(mdl_vals))
}

### SEARCH 
precision <<- 6 # number of decimal places considered 
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
    
    #Y_0 <- params[1]
    #Z_0 <- params[2]
    #alpha <- params[3]
    M <- params[1]
    tau <- params[2]
    alpha_MLT <- params[3]
    alpha_ov <- params[4]
    
    iteration <<- iteration + 1
    cat(paste("** iter:", iteration, "\n"))
    
    # check if parameters are reasonable 
    if (M < 0.05 || M > 2 || tau > log10(13.8*10**9) || 
            alpha_MLT < 0.5 || alpha_MLT > 3 || 
            alpha_ov < 0 || alpha_ov > 1 ) {
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**  ran:", ran, "\n")) 
    
    # call 'dispatch.sh' shell script to run MESA with supplied parameters 
    command <- paste("./dispatch.sh", 
        '-M', M, 
        '-t', 10**tau, 
        '-a', alpha_MLT,
        '-o', alpha_ov,
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
    control=list(trace=999, abstol=0.000001))#, parscale=parscale))

cat("Optimization terminated. Saving best result\n")
M <- result$par[1]
tau <- result$par[2]
alpha_MLT <- result$par[3]
alpha_ov <- result$par[4]
command <- paste("./dispatch.sh", 
    '-M', M, 
    '-t', 10**tau, 
    '-a', alpha_MLT, 
    '-o', alpha_ov, 
    '-m', -1, # write out profile 
    '-f', 1, # write out FGONG file 
    '-d', directory) 
print(command) 
system(command) 
print(result) 
if (result$convergence == 0) cat("Optimization successful.\n") 

