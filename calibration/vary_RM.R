#### Calibrate solar models with high and low mass and radius 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)
constraint.names = c("log L", "log R", "Fe/H", "Age")
solar_age = 4.572 * 10**9

### HELPERS
get_Z <- function(Y) {
    X <- -100000*(Y-1)/102293
    1-X-Y
}

### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
diffusion <- if (length(args)>0) as.numeric(args[1]) else 0 
mass <- if (length(args)>1) as.numeric(args[2]) else 1 
radius <- if (length(args)>2) as.numeric(args[3]) else 1 
directory <- paste0("models/M=",mass,"_R=",radius,"_D=",diffusion)
print(directory)

### INITIALIZE PARAMETERS FOR SEARCH
if (diffusion>0) {
    param_names <- c("Y", "alpha", "overshoot", "Z")
    param_init <- c(0.27202387, 1.84663590, 0.09104194, 0.01830403)
} else {
    param_names <- c("Y", "alpha", "overshoot")
    param_init <- c(0.2648378, 1.6974652, 0.2345875)
}

### DEFINE OBJECTIVE FUNCTION
chisq <- function() {
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    mdl.Fe.H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10_Z_div_X_solar
    mdl.vals <- c(mdl$log_L, 
        log10(10**mdl$log_R/radius), 
        mdl.Fe.H, 
        log10(mdl$star_age/solar_age)) 
    
    chi2 <- mdl.vals**2
    
    cat("** Model values and chi^2\n")
    cat(paste("Teff", 10**mdl$log_Teff, "\n"))
    cat(paste("log g", mdl$log_g, "\n"))
    cat(paste(constraint.names, mdl.vals, chi2, "\n"))
    
    chi2
}

### SEARCH
iteration <<- 0
ran <<- 0
best_chisq <<- Inf
best_param <<- param_init

run <- function(params) {
    cat(paste(param_names, params, "\n"))
    
    Y <- params[1]
    alpha <- params[2]
    overshoot <- params[3]
    Z <- if (diffusion) params[4] else get_Z(Y)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.2463     || Y > 0.33  || 
        Z < 0.01       || Z > 0.03 || 
        overshoot <= 0 || overshoot > 0.5 || 
        alpha < 1      || alpha > 3) {
        
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh",
        '-M', mass,
        '-Y', Y, 
        '-a', alpha, 
        '-Z', Z, 
        '-o', overshoot, 
        '-D', diffusion, '-d', directory)
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

print(result)

cat("Optimization terminated. Saving best result\n")
Y_0 <- result$par[1]
alpha <- result$par[2]
overshoot <- params[3]
Z_0 <- if (length(par)>=4) result$par[4] else get_Z(Y_0)
command <- paste("./dispatch.sh", 
    '-M', mass,
    '-Y', Y_0, 
    '-Z', Z_0, 
    '-o', overshoot,
    '-a', alpha, 
    '-d', directory, 
    '-m', mesh_delta_coeff, 
    '-v', varcontrol_target) 
print(command) 
system(command) 
print(result) 
if (result$convergence == 0) cat("Optimization successful.\n") 

