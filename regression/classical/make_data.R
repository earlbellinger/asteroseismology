args <- commandArgs(TRUE)
filename <- ifelse(length(args) > 0, args[1], 
    file.path('data', 'cep_lmc_v.dat'))

Band <- toupper(substr(filename, nchar(filename)-4, nchar(filename)-4))

#filename <- 'Table_RRL.dat'
DF <- read.table(filename, header=1)

for (ii in 1:nrow(DF)) {#
    print(paste0(ii, ' of ', nrow(DF), ' (', signif(ii/nrow(DF)*100,4), '%)'))
    row <- DF[ii,] 
    
    star.fname <- file.path('..', 'data', 'classical-new', 
            paste0(row$ID, '-obs.dat'))
    
    attach(row)
    
    #mode <- if (Mode == 'FU') 0 else if (Mode == 'FO') 1
    star.info <- data.frame(
        name=c('Period'),
        value=c(Period), 
        uncertainty=c(0.001))
    
    star <- data.frame(
        name=paste0(Band, '_', 
             c('Amplitude', 'Skewness', 'acuteness', 'R21', 'R31', 'A1', 'A2', 'A3')),
        value=c(Amplitude,   Skewness,   acuteness,   R21,   R31,   A1,   A2,   A3),
        uncertainty=c(sqrt(2)*eA0, min(1e-3, rms), min(1e-3, rms),  eR21,  eR31,  eA1,  eA2,  eA3))
        #name=c(paste0(Band, '_A'), paste0(Band, '_Sk'), paste0(Band, '_Ac'), paste0(Band, '_mean_mag'), paste0(Band, '_R21'), paste0(Band, '_R31'), paste0(Band, '_P21'), paste0(Band, '_P31')),
        #value=c(              A,                  Sk,                  Ac,                  mean_mag,                  R21,                  R31,                  P21,                  P31),
        #uncertainty=c(        0,                   0,                   0,                  mean_magErr,               R21Err,               R31Err,               P21Err,               P31Err))
    
    if (!file.exists(star.fname)) {
        star <- rbind(star.info, star)
        #fileConn <- file(file.path('..', 'data', 'classical', 
        #        paste0(row$StarID, '-freqs.dat')))
        #writeLines(c("l n nu dnu"), fileConn)
        #close(fileConn)
    } else {
        star <- rbind(read.table(star.fname, header=1), star)
    }
    
    detach(row)
    
    write.table(star, 
        file=star.fname, 
        quote=F, row.names=F, col.names=T)
}


#        star <- data.frame(
#            name=c('Band', 'Mode', 'logP', 'A', 'Sk', 'Ac', 'mean_mag',  'R21',  'R31',  'P21',  'P31'),
#            value=c(band,   mode,   logP,   A,   Sk,   Ac,   mean_mag,    R21,    R31,    P21,    P31),
#            uncertainty=c(0,   0,      0,   0,    0,    0,   mean_magErr, R21Err, R31Err, P21Err, P31Err))
