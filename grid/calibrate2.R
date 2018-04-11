#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=0)
library(optimx)

### CONSTANTS
Z_div_X_solar = 0.02293 # GS98 # 0.0245 # GN93 #
log10_Z_div_X_solar = log10(Z_div_X_solar) 
constraint.names = c("log L", "log R", "Fe/H") 

### INITIALIZE PARAMETERS FOR SEARCH 
directory <- file.path('calibrate') 
param_names <- c("Y", "alpha", "Z") 
param_init <- c(0.273449170177157, 1.83413390909832, 0.0197444964340224) 
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
    mdl_Fe_H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10_Z_div_X_solar 
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
    Z         <- params[3]
    
    if (Y < 0.25)  Y <- 0.25
    if (Y > 0.32)  Y <- 0.32
    if (Z < 0.012) Z <- 0.012
    if (Z > 0.03)  Z <- 0.03
    if (alpha < 1) alpha <- 1
    if (alpha > 3) alpha <- 3
    
    cat(paste(param_names, c(Y, alpha, Z), '\n'))
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh", 
        '-Y', Y, 
        '-a', alpha, 
        '-o', 0, 
        #'-f', 0, 
        '-Z', Z, 
        '-D', 1, 
        '-g', 1,
        '-e', 0,
        '-c', "4572000000", 
        '-d', directory)
    print(command)
    system(command)
    
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

#result <- optim(par=param_init, fn=run, method="Nelder-Mead", #"BFGS", 
#    control=list(trace=999, abstol=1e-10, maxit=10000,
#                 reltol=sqrt(.Machine$double.eps)/2,
#                 parscale=param_init/50))

result <- optimx(par=param_init, fn=run, method="Nelder-Mead")#,#"L-BFGS-B", #
    #lower=c(0.26, 1.5, 0.012), upper=c(0.32, 2.5, 0.03))

cat("Optimization terminated. Saving best result\n")
Y         <- result$p1#par[1]
alpha     <- result$p2#par[2]
Z         <- result$p3#par[3]
command   <- paste("./dispatch.sh", 
    '-Y', Y, 
    '-a', alpha, 
    '-o', 0, 
    #'-f', 0, 
    '-Z', Z, 
    '-D', 1, 
    '-g', 1,
    '-e', 0,
    '-c', "4572000000", 
    '-d', directory) 
print(command) 
system(command) 
print(result)
#if (result$conv == 0) cat("Optimization successful.\n") 

