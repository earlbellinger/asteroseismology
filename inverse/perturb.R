#### Generate Monte-Carlo perturbations of observed stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source('../scripts/seismology.R')
library(parallel)
library(parallelMap)

dir.create('perturb', showWarnings=FALSE)
dir.create('perturb/kages', showWarnings=FALSE)
dir.create('perturb/hares', showWarnings=FALSE)

speed_of_light = 299792 # km/s

### Obtain properties of real stars varied within their uncertainties 
monte_carlo_perturbations <- function(star, obs_data_file, freqs_data_file,
        n_perturbations=10000) {
    
    freqs <- read.table(freqs_data_file, header=TRUE)
    obs_data <- read.table(obs_data_file, header=TRUE)
    
    seismology(freqs, obs_data[obs_data$name == 'nu_max',]$value, outf=star,
        filepath=file.path('plots', 'perturb'))
    
    noisy_freqs <- freqs
    parallelStartMulticore(max(1, detectCores()))
    do.call(rbind, with(obs_data, {
        parallelMap(function(n) {
        #Map(function(n) {
            # Perturb observations by their uncertainties
            obs.DF <- data.frame(rbind(rnorm(nrow(obs_data), value, 
                if (n==1) 0 else uncertainty)))
            colnames(obs.DF) <- name
            
            # Correct frequencies for Doppler shift
            radial_velocity <- 
                if (any(grepl("radial_velocity", obs_data$name))) {
                    rnorm(1, value[name=="radial_velocity"], 
                    ifelse(n==1, 0, uncertainty[name=="radial_velocity"]))
                } else 0
            doppler_beta <- radial_velocity/speed_of_light
            doppler_shift <- sqrt((1+doppler_beta)/(1-doppler_beta))
            
            # Perturb frequencies
            noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu * doppler_shift, 
                ifelse(n==1, 0, freqs$dnu * doppler_shift))
            
            # Calculate Dnu, dnus, and ratios
            seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
            merge(rbind(obs.DF), rbind(seis.DF))
        }, 1:n_perturbations)
    }))
}

process <- function(star, star_dir, out_dir="perturb") {
    obs_data_file <- file.path(star_dir, paste0(star, "-obs.dat"))
    freqs_data_file <- file.path(star_dir, paste0(star, "-freqs.dat"))
    write.table(
        monte_carlo_perturbations(star, obs_data_file, freqs_data_file), 
        file.path(out_dir, paste0(star, "_perturb.dat")), 
        quote=FALSE, sep='\t', row.names=FALSE)
}

# Perturb every star 10k times and save the results
star_dir <- file.path("data", "hares")
for (fname in list.files(star_dir)) {
    if (!grepl('-obs.dat', fname)) next
    star <- sub('-obs.dat', '', fname)
    process(star, star_dir, out_dir="perturb/hares")
}

star_dir <- file.path("data", "kages")
for (fname in list.files(star_dir)) {
    if (!grepl('-obs.dat', fname)) next
    star <- sub('-obs.dat', '', fname)
    process(star, star_dir, out_dir="perturb/kages")
}

star_names <- c("Tagesstern", "16CygA", "16CygB", "Sun")
star_dir <- file.path("data")
for (star in star_names) process(star, star_dir)
