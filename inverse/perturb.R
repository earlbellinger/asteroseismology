#### Generate Monte-Carlo perturbations of observed stars 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'seismology.R'))
library(parallel)
library(parallelMap)

dir.create('perturb', showWarnings=FALSE)

speed_of_light = 299792 # km/s

### Obtain properties of real stars varied within their uncertainties 
perturb <- function(star, obs_data_file, freqs_data_file,
        n_perturbations=10000) {
    obs_data <<- read.table(obs_data_file, header=TRUE) 
    freqs <<- read.table(freqs_data_file, header=TRUE) 
    seis.DF <- seismology(freqs, obs_data[obs_data$name == 'nu_max',]$value, 
        outf=star, filepath=file.path('plots', 'perturb')) 
    cols <<- length(seis.DF)
    do.call(plyr:::rbind.fill, parallelMap(rand_inst, 1:n_perturbations)) 
}

rand_inst <- function(n) {
    # Perturb observations by their uncertainties
    attach(obs_data)
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
    noisy_freqs <- freqs
    repeat {
        noisy_freqs$nu <- rnorm(nrow(freqs), freqs$nu * doppler_shift, 
            ifelse(n==1, 0, freqs$dnu * doppler_shift))
        # Calculate Dnu, dnus, and ratios
        seis.DF <- seismology(noisy_freqs, obs.DF$nu_max)
        if (all(!is.na(seis.DF)) && length(seis.DF) == cols) break
    }
    detach(obs_data)
    merge(rbind(obs.DF), rbind(seis.DF))
}

process_dir <- function(star_dir) {
    out_dir <- file.path('perturb', basename(star_dir))
    dir.create(out_dir, showWarnings=FALSE)
    fnames <- list.files(star_dir)
    for (fname in fnames) {
        if (!grepl('-obs.dat', fname)) next
        star <- sub('-obs.dat', '', fname)
        process_star(star, star_dir, out_dir=out_dir)
    }
}

process_star <- function(star, star_dir, out_dir="perturb") {
    print(paste("Processing", star))
    obs_data_file <- file.path(star_dir, paste0(star, "-obs.dat"))
    freqs_data_file <- file.path(star_dir, paste0(star, "-freqs.dat"))
    result <- perturb(star, obs_data_file, freqs_data_file)
    if (!is.null(result))  
        write.table(result, 
            file.path(out_dir, paste0(star, "_perturb.dat")), 
            quote=FALSE, sep='\t', row.names=FALSE)
}

## Perturb every star 10k times and save the results
parallelStartMulticore(max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))
process_dir(file.path("data", "kages"))
process_dir(file.path("data", "basu"))
process_dir(file.path("data", "hares"))

star_names <- c("Tagesstern", "16CygA", "16CygB", "Sun")
star_dir <- file.path("data")
for (star in star_names) process_star(star, star_dir)

