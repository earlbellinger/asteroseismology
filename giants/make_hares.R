#### Collect evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

Z_div_X_solar = 0.02293

indir <- file.path('hares',
'M=1.6_Y=0.274452820062591_Z=0.0186941465808226_alpha=1.81511095685771_overshoot=0.005_diffusion=0',
'LOGS_SG')

outdir <- file.path('..', 'inverse', 'data', 'sg-hares-mesa')
dir.create(outdir)

hstry <- read.table(file.path(indir, 'history.data'), skip=5, header=1)

truth <- with(hstry, data.frame(
    model_number, age=star_age/10**9, log_g, radius=10**log_R, 
    Y_surf=surface_he4+surface_he3,
    tau_MS=star_age/min(star_age)
))

write.table(truth, file=file.path(outdir, '..', 'sg-hares-mesa.dat'),
    quote=F, row.names=F)

obs <- with(hstry, data.frame(
    model_number, Teff=10**log_Teff, L=10**log_L, nu_max, 
    Fe_H=log10(10**log_surf_cell_z / surface_h1 / Z_div_X_solar)
))

files <- list.files(indir)
freq.files <- files[grepl('-freqs.dat$', files)]
#prof.files <- sub('-freqs.dat', '.data', freq.files)

indices <- read.table(file.path(indir, 'profiles.index'), skip=1,
    col.names=c('model_number', 'priority', 'profile_number'))
obs. <- merge(obs, indices, by="model_number")


for (freq.file in freq.files) {
    profile_number <- as.numeric(
        sub('-freqs.dat', '', sub('profile', '', freq.file)))
    vals <- obs.[obs.$profile_number == profile_number,]
    model_number <- vals$model_number
    Teff <- vals$Teff #round(rnorm(1, vals$Teff, 85), 1)
    L <- vals$L #round(rnorm(1, vals$L, .03*vals$L), 1)
    Fe_H <- vals$Fe_H #round(rnorm(1, vals$Fe_H, 0.09), 3)
    nu_max <- vals$nu_max #round(rnorm(1, vals$nu_max, .05*vals$nu_max), 2)
    
    write.table(data.frame(name=c("Teff", "L", "Fe_H", "nu_max"),
            value=c(Teff, L, Fe_H, nu_max),
            uncertainty=c(85, round(.03*L, 2), 0.09, round(.05*nu_max, 2))
        ), file=file.path(outdir, paste0(model_number, '-obs.dat')), 
        quote=F, row.names=F)
    
    freqs <- read.table(file.path(indir, freq.file), skip=5, header=1)
    dnu <- ifelse(freqs$l <= 1, 3.5e-05, .0007)
    nus <- rnorm(nrow(freqs), freqs$Re.freq., dnu)
    
    write.table(with(freqs, data.frame(l, n=n_pg, n_p, n_g, 
            nu=nus, dnu=nus*dnu)),
        file=file.path(outdir, paste0(model_number, '-freqs.dat')),
        quote=F, row.names=F)
}

