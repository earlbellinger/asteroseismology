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
n_perturbations = 1000
classical = F 

### Obtain properties of real stars varied within their uncertainties 
perturb <- function(star, obs_data_file, freqs_data_file, 
        n_perturbations=n_perturbations, bias=0, imp=0, feh=T,
        feh.bias=0, teff.bias=0, feh.imp=0, teff.imp=0) {
    
    obs_data <<- read.table(obs_data_file, header=TRUE) 
    
    feh <- if (feh) obs_data$name == 'Fe/H' else obs_data$name == 'Teff'
    obs_data[feh,]$value <<- obs_data[feh,]$value + bias
    obs_data[feh,]$uncertainty <<- obs_data[feh,]$uncertainty + imp
    
    if (feh.bias != 0)
        obs_data[obs_data$name == 'Fe/H',]$value <<- 
            obs_data[obs_data$name == 'Fe/H',]$value + feh.bias
    
    if (teff.bias != 0)
        obs_data[obs_data$name == 'Teff',]$value <<- 
            obs_data[obs_data$name == 'Teff',]$value + teff.bias
    
    if (feh.imp != 0)
        obs_data[obs_data$name == 'Fe/H',]$uncertainty <<- 
            obs_data[obs_data$name == 'Fe/H',]$uncertainty + feh.imp
    
    if (teff.imp != 0)
        obs_data[obs_data$name == 'Teff',]$uncertainty <<- 
            obs_data[obs_data$name == 'Teff',]$uncertainty + teff.imp
    
    obs.vals <<- data.frame(rbind(obs_data$value)) 
    colnames(obs.vals) <- obs_data$name 
    
    obs.unc <<- data.frame(rbind(obs_data$uncertainty)) 
    colnames(obs.unc) <- obs_data$name 
    
    if (file.exists(freqs_data_file)) 
        freqs <<- read.table(freqs_data_file, header=TRUE)
    else freqs <<- data.frame()
    
    if (classical && 
            !('V_R21' %in% obs_data$name && 'I_R21' %in% obs_data$name)) {
        print(paste0('Skipping ', obs_data_file))
        return(NULL)
    }
    
    if (nrow(freqs) > 0) {
        nu_max <<- obs_data[obs_data$name == 'nu_max',]$value
        seis.DF <- seismology(freqs, nu_max, verbose=F)
    } else seis.DF <- NULL
    
    cols <<- length(seis.DF)
    start.time <- proc.time()
    set.seed(0)
    res <- do.call(plyr:::rbind.fill, 
        parallelMap(rand_inst, 1:n_perturbations)) 
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
        if ('V_P21' %in% names(obs.DF)) obs.DF$V_P21 <- obs.DF$V_P21 %% (2*pi) 
        # && obs.DF$V_P21 < 0) next
        if ('I_P21' %in% names(obs.DF)) obs.DF$I_P21 <- obs.DF$I_P21 %% (2*pi) 
        # && obs.DF$I_P21 < 0) next
        if ('V_P31' %in% names(obs.DF)) obs.DF$V_P31 <- obs.DF$V_P31 %% (2*pi) 
        # && obs.DF$V_P31 < 0) next
        if ('V_P31' %in% names(obs.DF)) obs.DF$I_P31 <- obs.DF$I_P31 %% (2*pi) 
        # && obs.DF$I_P31 < 0) next
        
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
            if (any(noisy_freqs$nu < freqs$nu*doppler_shift
                        - 5*freqs$dnu*doppler_shift |
                    noisy_freqs$nu > freqs$nu*doppler_shift
                        + 5*freqs$dnu*doppler_shift)) next 
            # Calculate Dnu, dnus, and ratios
            seis.DF <- seismology(noisy_freqs, obs.DF$nu_max, verbose=F)
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

process_dir <- function(star_dir, 
        feh.bias=0, teff.bias=0, 
        feh.imp=0, teff.imp=0,
        bias=0, 
        imp=0, feh=T) {
    out_dir <- file.path('perturb', basename(star_dir))
    if (bias != 0) out_dir <- paste0(out_dir, '_bias', bias)
    if (imp != 0) out_dir <- paste0(out_dir, '_imp', imp)
    if (feh.bias != 0) out_dir <- paste0(out_dir, '_fehbias', feh.bias)
    if (teff.bias != 0) out_dir <- paste0(out_dir, '_teffbias', teff.bias)
    if (feh.imp != 0) out_dir <- paste0(out_dir, '_fehimp', feh.imp)
    if (teff.imp != 0) out_dir <- paste0(out_dir, '_teffimp', teff.imp)
    if (file.exists(out_dir)) return()
    dir.create(out_dir, showWarnings=FALSE)
    fnames <- list.files(star_dir)
    times <- c()
    if (F) {
    for (fname in fnames) {
        if (!grepl('-obs.dat', fname)) next
        star <- sub('-obs.dat', '', fname)
        times <- c(times, 
            process_star(star, star_dir, out_dir=out_dir, 
                bias=bias, imp=imp, feh=feh, 
                feh.bias=feh.bias, teff.bias=teff.bias,
                feh.imp=feh.imp, teff.imp=teff.imp))
    }
    } else {
        times <- do.call(c, parallelMap(function(fname) {
            star <- sub('-obs.dat', '', fname)
            process_star(star, star_dir, out_dir=out_dir, 
                bias=bias, imp=imp, feh=feh,
                feh.bias=feh.bias, teff.bias=teff.bias,
                feh.imp=feh.imp, teff.imp=teff.imp)
        }, fname=fnames[grepl('obs.dat', fnames)]))
    }
    print(paste("Average time:", mean(times), "+/-", sqrt(var(times))))
    cat('\n\n')
}

process_star <- function(star, star_dir, out_dir="perturb", bias=0, imp=0,
        feh=T, feh.bias=0, teff.bias=0, feh.imp=0, teff.imp=0) {
    print(paste("Processing", star))
    obs_data_file <- file.path(star_dir, paste0(star, "-obs.dat"))
    freqs_data_file <- file.path(star_dir, paste0(star, "-freqs.dat"))
    results <- perturb(star, obs_data_file, freqs_data_file, n_perturbations,
        bias=bias, imp=imp, feh=feh, feh.bias=feh.bias, teff.bias=teff.bias,
        feh.imp=feh.imp, teff.imp=teff.imp)
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


for (feh.bias in c(-0.1, -0.075, -0.05, -0.025, -0.01,
                    0,
                    0.1,  0.075,  0.05,  0.025,  0.01)) {
    for (teff.bias in c(-10, 10)) {
        print(c(feh.bias, teff.bias))
        process_dir(file.path("data", "both"), 
            feh.bias=feh.bias, teff.bias=teff.bias)
    }
}
for (feh.bias in c(-0.01, 0.01)) {
    for (teff.bias in c(-100, -75, -50, -25, -10,
                         0,
                         100,  75,  50,  25,  10)) {
        print(c(feh.bias, teff.bias))
        process_dir(file.path("data", "both"), 
            feh.bias=feh.bias, teff.bias=teff.bias)
    }
}

stop()


#process_dir(file.path("data", "final"))
#process_dir(file.path("data", "feh"))

#stop()

# Teff experiment 

imps <- unique(c(#seq(0.1, 1, 0.1)#, 
                 #1:10, 
                 #seq(10, 100, 10),
                 seq(100, 1000, 100)
                 ))
#for (imp in imps) {
parallelMap(function(imp) {
    print(c('imp teff', imp))
    process_dir(file.path("data", "teff"), imp=imp, feh=F)
}, imp=imps)

#stop()

biases <- unique(c(seq(-1000, -100, 100)#,
                   #seq(-100, -10, 10),
                   #seq(-10, -1, 1),
                   #seq(-1, -0.1, 0.1)
                ))
biases <- sort(c(biases, abs(biases)))
#for (bias in biases) {
parallelMap(function(bias) {
    print(c('bias teff', bias))
    process_dir(file.path("data", "teff"), bias=bias, feh=F)
}, bias=biases)

# [Fe/H] experiment
biases <- unique(c(seq(-1, -0.1, 0.1)#,
                   #seq(-0.1, -0.01, 0.01), 
                   #seq(-0.01, 0.001, 0.001),
                   #seq(-0.001, 0.0001, 0.0001)
                   ))
biases <- c(biases, abs(biases))
#for (bias in biases) {
parallelMap(function(bias) {
    print(c('bias feh', bias))
    process_dir(file.path("data", "feh"), bias=bias)
}, bias=biases)
imps <- unique(c(#seq(0.0001, 0.001, 0.0001),
                 #seq(0.001, 0.01, 0.001),
                 #seq(0.01, 0.1, 0.01),
                 seq(0.1, 1, 0.1)
                ))
#for (imp in imps) {
parallelMap(function(imp) {
    print(c('imp feh', imp))
    process_dir(file.path("data", "feh"), imp=imp)
}, imp=imps)
#}



stop() 

for (teff.imp in c(0, 20, 40, 60, 80, 100)) {
    for (feh.imp in c(0, 0.02, 0.04, 0.06, 0.08, 0.1)) {
        print(c(feh.imp, teff.imp))
        process_dir(file.path("data", "both"), 
            feh.imp=feh.imp, teff.imp=teff.imp)
    }
}

for (feh.bias in c(-0.1, -0.075, -0.05, -0.025, 
                    0,
                    0.1,  0.075,  0.05,  0.025)) {
    for (teff.bias in c(-100, -75, -50, -25,
                         0,
                         100,  75,  50,  25)) {
        print(c(feh.bias, teff.bias))
        process_dir(file.path("data", "both"), 
            feh.bias=feh.bias, teff.bias=teff.bias)
    }
}

#stop()



#process_dir(file.path("data", "LMC_CEP"))
#process_dir(file.path('data', 'sun'))
stop()
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

