DF <- read.table('LOGS/history.data', skip=5, header=1)
decreasing_L <- which(diff(DF$log_L) < 0 & DF$star_age[-1] < 0.25*10**9)
if (any(decreasing_L)) {
    goes_back_up <- diff(decreasing_L) > 1
    pms <- ifelse(any(goes_back_up), 
               which(goes_back_up)[1] + 1, 
               max(decreasing_L))
    DF <- DF[-1:-pms,]
}
DF <- DF[floor(nrow(DF)*0.05):floor(nrow(DF)-nrow(DF)*0.05),]
#print(length(boxplot.stats(diff(DF$log_Teff), coef=3)$out) > 0)
cat(as.numeric(length(boxplot.stats(diff(DF$log_Teff), coef=3)$out) > 0))
