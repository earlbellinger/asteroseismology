#### Generate Monte-Carlo perturbations of observed stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'seismology.R'))
library(parallel)
library(parallelMap)
library(pracma)

#options(warn=1)

dir.create('perturb', showWarnings=FALSE)

speed_of_light = 299792 # km/s
boltz_sigma = 5.670367e-5 # erg cm^-2 K^-4 s^-1
R_solar = 6.957e10 # cm
L_solar = 3.828e33 # erg s^-1
n_perturbations = 10000
classical = F 

### Obtain properties of real stars varied within their uncertainties 
perturb <- function(star, obs_data_file, freqs_data_file, n_perturbations) {
    obs_data <<- read.table(obs_data_file, header=TRUE) 
    obs.vals <<- data.frame(rbind(obs_data$value)) 
    colnames(obs.vals) <- obs_data$name 
    obs.unc <<- data.frame(rbind(obs_data$uncertainty)) 
    colnames(obs.unc) <- obs_data$name 
    #if ('Teff' %in% obs_data$name 
    #     & 'L' %in% obs_data$name 
    #& 'radius' %in% obs_data$name) {
    #    obs_data <<- obs_data[-which(obs_data$name == 'Teff'),]
    #}
    if (file.exists(freqs_data_file)) 
        freqs <<- read.table(freqs_data_file, header=TRUE)
    else freqs <<- data.frame()
    if (classical && !('V_R21' %in% obs_data$name && 'I_R21' %in% obs_data$name)) {
        print(paste0('Skipping ', obs_data_file))
        return(NULL)
    }
    if (nrow(freqs) > 0) {
        nu_max <<- obs_data[obs_data$name == 'nu_max',]$value
        seis.DF <- seismology(freqs, nu_max)
    } else seis.DF <- NULL
    cols <<- length(seis.DF)
    start.time <- proc.time()
    res <- do.call(plyr:::rbind.fill, parallelMap(rand_inst, 1:n_perturbations)) 
    total.time <- proc.time() - start.time
    time_per_perturb <- total.time[[3]]/n_perturbations
    print(paste("Total time:", total.time[[3]], 
                "; Time per perturbation:", time_per_perturb))
    return(list(res=res, time_per_perturb=time_per_perturb))
}

rand_inst <- function(n) {
    # Perturb observations by their uncertainties
    attach(obs_data)
    repeat {
        obs.DF <- data.frame(rbind(rnorm(nrow(obs_data), value, 
            if (n==1) 0 else uncertainty)))
        colnames(obs.DF) <- name
        
        if ('Dnu0' %in% names(obs.DF) && obs.DF$Dnu0 < 0) next 
        if ('nu_max' %in% names(obs.DF) && obs.DF$nu_max < 0) next 
        if ('V_R21' %in% names(obs.DF) && obs.DF$V_R21 < 0) next
        if ('I_R21' %in% names(obs.DF) && obs.DF$I_R21 < 0) next
        if ('V_P21' %in% names(obs.DF)) obs.DF$V_P21 <- obs.DF$V_P21 %% (2*pi) # && obs.DF$V_P21 < 0) next
        if ('I_P21' %in% names(obs.DF)) obs.DF$I_P21 <- obs.DF$I_P21 %% (2*pi) # && obs.DF$I_P21 < 0) next
        if ('V_P31' %in% names(obs.DF)) obs.DF$V_P31 <- obs.DF$V_P31 %% (2*pi) # && obs.DF$V_P31 < 0) next
        if ('V_P31' %in% names(obs.DF)) obs.DF$I_P31 <- obs.DF$I_P31 %% (2*pi) # && obs.DF$I_P31 < 0) next
        
        # Apply the Stefan-Boltzmann Law 
        if ('Teff' %in% names(obs.DF)#obs_data$name 
                & 'L' %in% names(obs.DF) 
                & 'radius' %in% names(obs.DF)) {
            obs.DF$Teff <- nthroot(obs.DF$L * L_solar / 
                (4*pi*(R_solar * obs.DF$radius)**2*boltz_sigma), 4)
        }
        
        if (any(obs.DF < obs.vals-4*obs.unc | obs.DF > obs.vals+4*obs.unc)) next
        
        # Correct frequencies for Doppler shift
        radial_velocity <- 
            if (any(grepl("radial_velocity", obs_data$name))) {
                rnorm(1, value[name=="radial_velocity"], 
                ifelse(n==1, 0, uncertainty[name=="radial_velocity"]))
            } else 0
        doppler_beta <- radial_velocity/speed_of_light
        doppler_shift <- sqrt((1+doppler_beta)/(1-doppler_beta))
        
        # Perturb frequencies
        if (cols > 0) {
            noisy_freqs <- freqs
            noisy_freqs$nu <- rnorm(nrow(freqs), 
                freqs$nu * doppler_shift, 
                if (n==1) 0 else freqs$dnu * doppler_shift)
            if (any(noisy_freqs$nu < 0)) next 
            if (any(noisy_freqs$nu < freqs$nu*doppler_shift-5*freqs$dnu*doppler_shift |
                    noisy_freqs$nu > freqs$nu*doppler_shift+5*freqs$dnu*doppler_shift)) next 
            # Calculate Dnu, dnus, and ratios
            seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
        } else seis.DF <- data.frame()
        if (all(!is.na(seis.DF)) && length(seis.DF) == cols) break
		#if (length(seis.DF) == cols) break
    }
    detach(obs_data)
    if (length(obs.DF) > 0 && length(seis.DF) > 0)
        merge(rbind(obs.DF), rbind(seis.DF))
    else if (length(obs.DF) > 0)
        rbind(obs.DF)
    else rbind(seis.DF) 
}

process_dir <- function(star_dir) {
    out_dir <- file.path('perturb', basename(star_dir))
    dir.create(out_dir, showWarnings=FALSE)
    fnames <- list.files(star_dir)
    times <- c()
    for (fname in fnames) {
        if (!grepl('-obs.dat', fname)) next
        star <- sub('-obs.dat', '', fname)
        times <- c(times, process_star(star, star_dir, out_dir=out_dir))
    }
    print(paste("Average time:", mean(times), "+/-", sqrt(var(times))))
    cat('\n\n')
}

process_star <- function(star, star_dir, out_dir="perturb") {
    print(paste("Processing", star))
    obs_data_file <- file.path(star_dir, paste0(star, "-obs.dat"))
    freqs_data_file <- file.path(star_dir, paste0(star, "-freqs.dat"))
    results <- perturb(star, obs_data_file, freqs_data_file, n_perturbations)
    result <- results$res
    if (!is.null(result)) {
        write.table(result, 
            file.path(out_dir, paste0(star, "_perturb.dat")), 
            quote=FALSE, sep='\t', row.names=FALSE)
    }
    results$time_per_perturb
}

## Perturb every star 10k times and save the results
#reg.finalizer(environment(), cleanup, onexit = FALSE)
parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
#process_dir(file.path('data', 'inversions'))
#process_dir(file.path('data', 'newsun'))
#process_dir(file.path('data', 'gangelou'))
#process_dir(file.path("data", "classical"))
#stop()
#process_dir(file.path("data", "sg-basu"))
#stop()
#process_dir(file.path("data", "inversions"))
#process_dir(file.path("data", "legacyRox2"))
#process_dir(file.path("data", "benard"))
#process_dir(file.path("data", "legacyRox"))
#process_dir(file.path("data", "sg-hares-mesa"))
#process_dir(file.path("data", "sg-hares"))
#process_dir(file.path("data", "procyon"))
#process_dir(file.path("data", "Dnu"))
#process_dir(file.path("data", "legacy"))
#process_dir(file.path("data", "kages"))
#process_dir(file.path("data", "basu"))
#process_dir(file.path("data", "hares"))

star_names <- c("16CygAamp", "16CygBamp", "16CygA", "16CygB", 
    "Tagesstern", "Sun")
star_dir <- file.path("data")
for (star in star_names) process_star(star, star_dir)

