library(BBmisc)
library(pracma)

boltz_sigma = 5.670367e-5 # erg cm^-2 K^-4 s^-1
R_solar = 6.957e10 # cm
L_solar = 3.828e33 # erg s^-1

args <- commandArgs(TRUE)
filename <- ifelse(length(args) > 0, args[1], 
    file.path('models', 'model_cep_fu-higherM.dat'))

DF <- read.table(filename, header=1)
DF2 <- NULL
#DF.infos <- NULL
#DF.models <- NULL


for (ii in 1:nrow(DF)) {#1:10) {#1:500) {#
    row <- DF[ii,]
    
    print(ii)
    
    if (!(row$Band == 'V' || row$Band == 'I')) next 
    
    attach(row)
    
    mode <- if (Mode == 'FU') 0 else if (Mode == 'FO') 1
    model.info <- if (grepl('cep', filename)) {
        data.frame(M=Mass, Z=Z, logR=Radius, 
            logL=Luminosity, 
            Teff=Temperature,#nthroot(10**Luminosity * L_solar / 
                #(4*pi*(R_solar * 10**Radius)**2*boltz_sigma), 4), 
            Mode=mode, Period=Period)
    } else if (grepl('rr', filename)) {
        data.frame(M=Mass, Y=Y, Z=Z, logR=Radius, 
            logL=Luminosity, Teff=Temperature, Mode=mode, Period=Period)
    } else {
        data.frame(M=M, Y=Y, Z=Z, logR=logR, 
            L=L, Teff=Teff, Mode=mode, logP=logP)
    }
    print(paste(Temperature, model.info$Teff, Temperature-model.info$Teff))
    
    already <- if (!is.null(DF2)) {
        matches <- sapply(1:nrow(DF2), 
            function(ii) all(model.info==DF2[ii,1:ncol(model.info)]))
        ifelse(any(matches), which(matches)[1], F)
    } else F
    
    model <- if (grepl('cep', filename) || grepl('rr', filename)) {
        data.frame(Amplitude, Skewness, acuteness, R21, R31, M0=A0, A1, A2, A3)
    } else {
        data.frame(A, Sk, Ac, mmag, R21, R31, P21, P31)
        #colnames(model) <- c(paste0(Band, '_A'), paste0(Band, '_Sk'), 
        #    paste0(Band, '_Ac'), paste0(Band, '_mean_mag'),
        #    paste0(Band, '_R21'), paste0(Band, '_R31'),
        #    paste0(Band, '_P21'), paste0(Band, '_P31'))
    }
    colnames(model) <- paste0(Band, '_', colnames(model))
    
    #already <- DF2$Z == Z & 
    #    DF2$Y == Y & 
    #    DF2$M == M & 
    #    DF2$Teff == Teff & 
    #    DF2$logP == logP & 
    #    DF2$logR == logR & 
    #    DF2$L == L
    
    if (is.null(DF2)) {
        DF2 <- cbind(model.info, model)
        #DF.models <- model
        #DF.infos <- model.info
    } else if (!any(already)) {
        DF2 <- plyr:::rbind.fill(DF2, cbind(model.info, model))
        #DF.models <- plyr:::rbind.fill(DF.models, model)
        #DF.infos <- plyr:::rbind.fill(DF.infos, model.info)
    } else {
        old <- DF2[already,]
        old <- old[,complete.cases(t(old))]
        new <- cbind(old, model)
        #if (ncol(DF2) == ncol(new)) {
        #    DF2[already,] <- new
        #} else {
            DF2 <- plyr:::rbind.fill(DF2[-already,], new)
        #}
    }
    
    detach(row)
}


#DF2$Color <- DF2$V_M0 - DF2$I_M0
DF2$W <- with(DF2, I_M0 - 1.55*(V_M0 - I_M0))
#write.table(sortByCol(DF2, names(DF2)), 'models_cep.dat', row.names=F, quote=F, sep='\t')
#write.table(sortByCol(DF2, names(DF2)), 'models_rrab.dat', row.names=F, quote=F, sep='\t')

