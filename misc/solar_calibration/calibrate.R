#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES
source('../../scripts/seismology.R')


### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)


### HELPERS
show10 <- function(x) format(round(x, 10), nsmall=10)

get_Z <- function(Y) {
    # log10(Z/X/0.02293) = 0
    # X+Y+Z = 1 hence...
    X <- -100000*(Y-1)/102293
    as.numeric(show10(1-X-Y))
}


### PARSE COMMAND LINE ARGS
args <- commandArgs(TRUE)
diffusion <- if (length(args)>0) 1 else 0
directory <- if (diffusion) "diffusion" else "no_diffusion"
print(directory)


### PARSE SOLAR DATA
sun_perturb <- read.table('../../inverse/perturb/Sun_perturb.dat', header=1)
sun_obs <- as.data.frame(Map(mean, sun_perturb))
sun_std <- sqrt(as.data.frame(Map(var, sun_perturb)))

constraint.names <- c(#"log Teff", "log g", 
    "log L", "log R", "Fe/H", "Age") #, "r02", "r13", "r10", "r01")
sun.log <- c(#1, 
    1, 1, 1, 0) #, 0, 0, 0, 0)

logger <- function(val, logged=sun.log) 
    ifelse(logged, val, log10(val))
std.logger <- function(val, std, logged=sun.log) 
    ifelse(logged, std, std/(val*log(10)))
log.mapper <- function(vals) 
    unlist(Map(logger, val=vals, logged=sun.log))
std.mapper <- function(vals, stds) 
    unlist(Map(std.logger, val=vals, std=stds, logged=sun.log))

# solar vals come from
# https://sites.google.com/site/mamajeksstarnotes/basic-astronomical-data-for-the-sun
sun.vals <- c(#3.76131, 4.43812, 
    0, 0, 0, 4.572)#, 
    #sun_obs$r02, sun_obs$r13, sun_obs$r10, sun_obs$r01)
sun.stds <- c(#0.0005, 0.00013, 
    1e-5, 1e-5, 0.002, 0.004)#, 
    #2*sun_std$r02, 2*sun_std$r13, 2*sun_std$r10, 2*sun_std$r01)

sun.vals.log <- log.mapper(sun.vals)
sun.stds.log <- std.mapper(sun.vals, sun.stds)

print("** Matching solar values:")
cat(paste(constraint.names, show10(sun.vals), "+/-", show10(sun.stds), "\n"))


### INITIALIZE PARAMETERS FOR SEARCH
overshoot <- 0.005

if (diffusion) {
    param_names <- c("Y", "alpha", "Z")
    param_init <- c(0.27404954, 2.03315750, 0.01867976) #c(0.263863, 1.994654, 0.016609)
} else {
    param_names <- c("Y", "alpha")
    param_init <- c(0.2651741, 1.8809009) # c(0.265414, 1.877870)
}

#param_init <- c(Y, Z, alpha)#, overshoot)
#param_names <- c("Y", "Z", "alpha")#, "ov")


### DEFINE OBJECTIVE FUNCTION
chisq <- function() {
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    #mdl.Fe.H <- log10(10**mdl$log_surf_cell_z / mdl$surface_h1 / Z_div_X_solar)
    mdl.Fe.H <- mdl$log_surf_cell_z-log10(mdl$surface_h1)-log10_Z_div_X_solar
    
    #mdl.freqs <- read.table(
    #    file.path(directory, 'LOGS_MS', 'profile1-freqs.dat'), 
    #    col.names=c('l', 'n', 'nu'))
    #seis <- seismology(mdl.freqs, nu_max=mdl$nu_max)
    
    prof <- read.table(file.path(directory, 'LOGS_MS', 'profile1.data'), 
        header=1, skip=5)
    mdl.R_cz <- prof$radius[which(diff(prof$mlt_mixing_type)==-1)[1]]
    
    mdl.vals <- c(#mdl$log_Teff, mdl$log_g, 
        mdl$log_L, mdl$log_R, mdl.Fe.H, 
        mdl$star_age/10**9)#, seis$r02, seis$r13, seis$r10, seis$r01)
    mdl.vals.log <- log.mapper(mdl.vals)
    
    chi2 <- (sun.vals.log - mdl.vals.log)**2 / sun.stds.log
    
    cat("** Model values and chi^2\n")
    #cat(paste("log R", show10(mdl$log_R), '\n'))
    cat(paste("Teff", show10(10**mdl$log_Teff), "\n"))
    cat(paste("log g", show10(mdl$log_g), "\n"))
    cat(paste("R_cz", show10(mdl.R_cz), "\n"))
    cat(paste(constraint.names, show10(mdl.vals), show10(chi2), "\n"))
    
    #c((log10(sun_obs$Teff) - mdl$log_Teff)**2   / (sun_std$Teff   / (log(10)*5)), #sun_obs$Teff)), 
    #   mdl$log_L**2                             / (sun_std$L      / (log(10)*1e-5)), #sun_obs$L)), 
    #  (sun_obs$log_g - mdl$log_g)**2            /  sun_std$log_g,
    #   mdl.Fe.H**2                              /  0.02, #sun_std$Fe.H,
    #   #mdl$log_R**2                             / (sun_std$radius / (log(10)*sun_obs$radius)),
    #  (4.57     - log10(mdl$star_age/10**9))**2 / (4.57           / (log(10)*5)), 
    #  (log10(sun_obs$r02) - log10(seis$r02))**2 / (sun_std$r02    / (log(10)*sun_obs$r02)), 
    #  (log10(sun_obs$r13) - log10(seis$r13))**2 / (sun_std$r13    / (log(10)*sun_obs$r13)), 
    #  (log10(sun_obs$r10) - log10(seis$r10))**2 / (sun_std$r10    / (log(10)*sun_obs$r10)) / 2, 
    #  (log10(sun_obs$r01) - log10(seis$r01))**2 / (sun_std$r01    / (log(10)*sun_obs$r01)) / 2)
    
    chi2
}


### SEARCH
iteration <<- 0
ran <<- 0
best_chisq <<- Inf
best_param <<- param_init

run <- function(params) {
    params <- as.numeric(show10(params))
    cat(paste(param_names, show10(params), "\n"))
    
    Y <- params[1]
    alpha <- params[2]
    Z <- if (diffusion) params[3] else get_Z(Y)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.24 || Y > 0.34 || Z < 0.005 || Z > 0.035
        || alpha < 1 || alpha > 3) return(Inf)
    
    ran <<- ran + 1 
    cat(paste("**** ran:", ran, "\n"))
    
    command <- paste("./dispatch.sh",
        '-Y', Y, 
        '-a', alpha,
        '-Z', Z,  
        '-o', overshoot, #params[4])
        '-D', diffusion, '-d', directory)
    print(command)
    system(command)#, ignore.stdout=T)
    
    chi_sq <- sum(chisq())
    cat(paste('**** chi_sq =', show10(chi_sq), "\n"))
    
    if (chi_sq < best_chisq) {
        best_chisq <<- chi_sq
        cat(paste(c("*****", param_names, show10(params), "\n"))) #, "ov"), params))
        best_param <<- params
        print("***** New record!")
        system(paste0("cp -r ", directory, " ", directory, "_best"))
    }
    
    chi_sq
}

result <- optim(par=param_init, fn=run, control=list(trace=999)) #,
#    bounds=list(c(0.2, 0.4), c(0.00001, 0.1), c(1, 3), c(0, 0.5)))

print(best_chisq)
print(best_param)
print(result)

