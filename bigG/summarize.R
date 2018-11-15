#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus

source(file.path('..', 'scripts', 'seismology.R'))
source(file.path('..', 'scripts', 'utils.R'))

Z_div_X_solar = 0.02293

directory <- commandArgs(1)[1]
print(directory)

trackfile <- file.path(directory, 'track')
params.DF <- read.table(trackfile, header=1)
params.DF['age'] <- as.numeric(sub('d', 'e', params.DF$age))

pro_header <- read.table(file.path(directory, 'LOGS', 'profile1.data'), 
    header=TRUE, nrows=1, skip=1)

if (pro_header['star_age'] != params.DF['age']) {
    print('Ages do not match')
    print('Requested age: ')
    print(params.DF$age)
    print('Actual age: ')
    print(pro_header$age)
    stop()
}

hst <- read.table(file.path(directory, 'LOGS', 'history.data'),
    header=TRUE, skip=5)
hst <- hst[nrow(hst),]

obs.DF <- NULL
obs.DF['Teff'] <- pro_header['Teff']
obs.DF['L'] <- pro_header['photosphere_L']
obs.DF['radius'] <- pro_header['photosphere_r']
obs.DF['nu_max'] <- hst['nu_max']
obs.DF['Fe_H'] <- log10(10**hst$log_surf_cell_z / 
            (hst$surface_h1 + hst$surface_h2) / Z_div_X_solar)

freqs <- parse_freqs(file.path(directory, 'LOGS', 'profile1-freqs.dat'), gyre=T)
seis.DF <- seismology(freqs, nu_max=obs.DF$nu_max, all_freqs=T, all_ratios=T)

DF <- cbind(params.DF, obs.DF, seis.DF)

write.table(DF, paste0(directory, '.dat'), quote=FALSE, sep='\t', 
    row.names=FALSE)
