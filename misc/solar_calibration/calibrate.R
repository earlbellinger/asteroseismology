#### Calibrate a solar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES
#source('../../scripts/seismology.R')


### CONSTANTS
Z_div_X_solar = 0.02293
log10_Z_div_X_solar = log10(Z_div_X_solar)


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
    1e-6, 1e-6, 1e-6, 1e-10)#, 
    #2*sun_std$r02, 2*sun_std$r13, 2*sun_std$r10, 2*sun_std$r01)

sun.vals.log <- log.mapper(sun.vals)
sun.stds.log <- std.mapper(sun.vals, sun.stds)

print("** Matching solar values:")
cat(paste(constraint.names, sun.vals, "+/-", sun.stds, "\n"))


### INITIALIZE PARAMETERS FOR SEARCH
overshoot <- 0.005

if (diffusion) {
    param_names <- c("Y", "alpha", "overshoot", "Z")
    param_init <- c(0.271032787415861, 1.78385818489939, 0.269788190914154, 0.0184014107904859)
    #c(0.27082287, 1.82418236, 0.08263118, 0.01831911)
    #c(0.26929755, 1.80003929, 0.08212582, 0.01758726)
    #c(0.269284255505623, 1.8001409944105, 0.0798093575943592, 0.0175851807290451)
    #c(0.269293167805798, 1.80022524500671, 0.0797476449908921, 0.0175858314409672)
    #c(0.269303889796648, 1.80122252094736, 0.0755625620387972, 0.0175861063758966)
    #c(0.27012154215252, 1.78385615431702, 0.005, 0.0176650574205056)
    #c(0.274452820062591, 1.81511095685771, 0.0186941465808226) # Eddington
    #c(0.274263097979849, 2.17046834906094, 0.0186651396018428) # simple_photosphere
    #c(0.274407513903574, 2.02490269613316, 0.0186873947508156)
    #c(0.2742011377657840, 2.0288779958248542, 0.0187091666137872)
    #c(0.27404954, 2.03315750, 0.01867976) 
    #c(0.263863, 1.994654, 0.016609)
} else {
    param_names <- c("Y", "alpha", "overshoot")
    param_init <- c(0.262195489130086, 1.66478537654705, 0.234472799771336)
    #c(0.26510493, 1.69010270, 0.07192404)
    #c(0.265103082607226, 1.68954421641992, 0.0701074797612781)
    #c(0.265110964684205, 1.69041572979509, 0.0712533832973251)
    #c(0.265116471537859, 1.69115066597018, 0.0728631047775711)
    #c(0.265469923147422, 1.68717773375763, 0.005)
    #c(0.265461419951207, 1.87358013005573)
    #c(0.2651678, 1.8769572)
    #c(0.2651741, 1.8809009) 
    #c(0.265414, 1.877870)
}

#param_init <- c(Y, Z, alpha)#, overshoot)
#param_names <- c("Y", "Z", "alpha")#, "ov")


### DEFINE OBJECTIVE FUNCTION
chisq <- function() {
    hstry <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
        header=1, skip=5)
    mdl <- hstry[nrow(hstry),]
    
    #mdl.Fe.H <- log10(10**mdl$log_surf_cell_z / mdl$surface_h1 / Z_div_X_solar)
    mdl.Fe.H <- mdl$log_surf_cell_z-mdl$log_surface_h1-log10_Z_div_X_solar
    
    #mdl.freqs <- read.table(
    #    file.path(directory, 'LOGS_MS', 'profile1-freqs.dat'), 
    #    col.names=c('l', 'n', 'nu'))
    #seis <- seismology(mdl.freqs, nu_max=mdl$nu_max)
    
    #prof <- read.table(file.path(directory, 'LOGS_MS', 'profile1.data'), 
    #    header=1, skip=5)
    #mdl.R_cz <- prof$radius[which(diff(prof$mlt_mixing_type)==-1)[1]]
    
    mdl.vals <- c(#mdl$log_Teff, mdl$log_g, 
        round(mdl$log_L, 12), round(mdl$log_R, 12), round(mdl.Fe.H, 12), 
        mdl$star_age/10**9)#, seis$r02, seis$r13, seis$r10, seis$r01)
    mdl.vals.log <- log.mapper(mdl.vals)
    
    chi2 <- (sun.vals.log - mdl.vals.log)**2 #/ sun.stds.log
    
    cat("** Model values and chi^2\n")
    #cat(paste("log R", mdl$log_R, '\n'))
    cat(paste("Teff", 10**mdl$log_Teff, "\n"))
    cat(paste("log g", mdl$log_g, "\n"))
    #cat(paste("R_cz", mdl.R_cz, "\n"))
    cat(paste(constraint.names, mdl.vals, chi2, "\n"))
    
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
    #params <- as.numeric(params)
    cat(paste(param_names, params, "\n"))
    
    Y <- round(params[1], 12)
    alpha <- round(params[2], 12)
    overshoot <- round(params[3], 12)
    Z <- if (diffusion) round(params[4], 12) else round(get_Z(Y), 12)
    
    iteration <<- iteration + 1
    cat(paste("**** iter:", iteration, "\n"))
    
    if (Y < 0.26  || Y > 0.35  || 
        Z < 0.016 || Z > 0.025 || 
        overshoot < 0 || overshoot > 0.5 || 
        alpha < 1.5 || alpha > 2.5) return(Inf)
    
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
        cat(paste(c("*****", param_names, params, "\n"))) #, "ov"), params))
        best_param <<- params
        print("***** New record!")
        #system(paste0("cp -r ", directory, " ", directory, "_best"))
    }
    
    chi_sq
}

result <- optim(par=param_init, fn=run, control=list(trace=999)) #,
#    bounds=list(c(0.2, 0.4), c(0.00001, 0.1), c(1, 3), c(0, 0.5)))

print(result)

