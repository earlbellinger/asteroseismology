#### Collate Gaia data into a useful format for me 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 

gaia_filename <- 'gaia.dat'
con  <- file(gaia_filename, open = "r")
first_line <- readLines(con, n=1, warn=F)
close(con)
widths <- nchar(unlist(regmatches(first_line, gregexpr(".+?\\|", first_line))))
gaia.DF <- read.fwf(file=gaia_filename, widths=widths, skip=4)
colnames(gaia.DF) <- strsplit(
    gsub('\\s+', '', substring(first_line, 2)), '\\|')[[1]]


gaia.DF2 <- with(gaia.DF, data.frame(KIC=KIC_ID,
    L=LUM_VAL_GAIA, 
    e_L=ERROR_LUM_PERCENTILE_UPPER_GAIA-LUM_VAL_GAIA,
    radial_velocity=RADIAL_VELOCITY_GAIA,
    e_radial_velocity=ERROR_RADIAL_VELOCITY_GAIA,
    R=RADIUS_VAL_GAIA,
    e_R=(ERROR_RADIUS_PERCENTILE_UPPER_GAIA
           -ERROR_RADIUS_PERCENTILE_LOWER_GAIA)/2))
    #Teff=TEFF_VAL_GAIA,
    #e_Teff=(ERROR_TEFF_PERCENTILE_UPPER_GAIA
    #       -ERROR_TEFF_PERCENTILE_LOWER_GAIA)/2))
gaia.DF2 <- gaia.DF2[!(gaia.DF2$L < -9998),]

filenames <- list.files(file.path('data', 'gaia'), full.names=T)
obs <- filenames[grepl('obs', filenames)]
for (ob in obs) {
    KIC <- sub('-obs.dat', '', basename(ob))
    if (!KIC %in% gaia.DF2$KIC) next
    gaia <- gaia.DF2[gaia.DF2$KIC == KIC,]
    DF <- read.table(ob, header=1)
    if ('L' %in% DF$name) DF[-which(DF$name == 'L'),]
    DF <- rbind(DF, data.frame(name='L', value=gaia$L, uncertainty=gaia$e_L))
    
    if (gaia$radial_velocity > -9998) {
        if ('radial_velocity' %in% DF$name) 
            DF[-which(DF$name == 'radial_velocity'),]
        DF <- rbind(DF, data.frame(name='radial_velocity', 
            value=gaia$radial_velocity, uncertainty=gaia$e_radial_velocity))
    }
    
    if (gaia$R > -9998) {
        if ('radius' %in% DF$name) 
            DF[-which(DF$name == 'radius'),]
        DF <- rbind(DF, data.frame(name='radius', 
            value=gaia$R, uncertainty=gaia$e_R))
    }
    
    write.table(DF, ob, quote=F, row.names=F)
}

