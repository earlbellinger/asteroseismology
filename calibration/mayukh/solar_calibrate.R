#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)
constraint.names = c("log L", "log R", "Fe/H", "Age")

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
directory <- if (diffusion) "diffusion" else "no_diffusion"
print(directory)

### INITIALIZE PARAMETERS FOR SEARCH
if (diffusion) {
    param_names <- c("Y", "alpha", "overshoot", "Z")
    param_init <- c(0.27202387, 1.84663590, 0.09104194, 0.01830403)
    #c(0.27203565118994, 1.84491361918352, 0.0831545426805924, 
    #    0.0183069639782441)
    #c(0.27238323, 1.83904669, 0.08285504, 0.01831911)
} else {
    param_names <- c("Y", "alpha", "overshoot")
    param_init <- c(0.2648378, 1.6974652, 0.2345875)
    #c(0.26510493, 1.69010270, 0.07192404)
}

### DEFINE OBJECTIVE FUNCTION
chisq <- function() {
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    mdl.Fe.H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10_Z_div_X_solar
    mdl.vals <- c(mdl$log_L, mdl$log_R, mdl.Fe.H, 
        log10(mdl$star_age/10**9/4.572))
    
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
    
    Y <- round(params[1], 8)
    alpha <- round(params[2], 8)
    overshoot <- round(params[3], 8)
    Z <- if (diffusion) round(params[4], 8) else round(get_Z(Y), 8)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.2463  || Y > 0.33  || 
        Z < 0.016   || Z > 0.021 || 
        overshoot <= 0 || overshoot > 0.3 || 
        alpha < 1.5 || alpha > 2.5) {
        
        cat("Input parameters out of bounds; skipping\n")
        return(Inf)
    }
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh",
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

