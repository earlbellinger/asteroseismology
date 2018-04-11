library(BBmisc)

DF <- read.table('Table_models.dat', header=1)
DF2 <- NULL


for (ii in 1:nrow(DF)) {#1:500) {#
    row <- DF[ii,]
    
    print(ii)
    
    if (!(row$Band == 'V' || row$Band == 'I')) next 
    
    attach(row)
    
    mode <- if (Mode == 'FU') 0 else if (Mode == 'FO') 1
    model.info <- data.frame(M=M, Y=Y, Z=Z, logR=logR, L=L, Teff=Teff, Mode=mode, logP=logP)
    model <- data.frame(A, Sk, Ac, mmag, R21, R31, P21, P31)
    colnames(model) <- c(paste0(Band, '_A'), paste0(Band, '_Sk'), 
        paste0(Band, '_Ac'), paste0(Band, '_mean_mag'),
        paste0(Band, '_R21'), paste0(Band, '_R31'),
        paste0(Band, '_P21'), paste0(Band, '_P31'))
    
    already <- DF2$Z == Z & 
        DF2$Y == Y & 
        DF2$M == M & 
        DF2$Teff == Teff & 
        DF2$logP == logP & 
        DF2$logR == logR & 
        DF2$L == L
    
    if (is.null(DF2)) {
        DF2 <- cbind(model.info, model)
    } else if (!any(already)) {
        DF2 <- plyr:::rbind.fill(DF2, cbind(model.info, model))
    } else {
        old <- DF2[already,]
        old <- old[,complete.cases(t(old))]
        new <- cbind(old, model)
        #if (ncol(DF2) == ncol(new)) {
        #    DF2[already,] <- new
        #} else {
            DF2 <- plyr:::rbind.fill(DF2[-which(already),], new)
        #}
    }
    
    detach(row)
}

write.table(sortByCol(DF2, names(DF2)), 'models.dat', row.names=F, quote=F, sep='\t')

