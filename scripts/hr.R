library(magicaxis)

options(scipen=5)
args <- commandArgs(TRUE)
filename <- if (length(args)>0) args[1] else file.path('LOGS', 'history.data')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    #print(paste("Clipping", pms, "points"))
    #DF <- DF[-1:-pms,]
}

png('HR.png', width=800, height=600, type='cairo')#, family='Palatino')
par(mar=c(5,5,3,1), family='Palatino', cex=1.25, cex.lab=1.25, cex.axis=1.25)
plot(DF$log_Teff, DF$log_L, 
    type='l', axes=F, tcl=0, 
    xlim=rev(range(DF$log_Teff)), 
    ylim=range(DF$log_L), 
    xlab=expression("Temperature"~log(T["eff"]/K)), 
    ylab=expression("Luminosity"~log(L/L["Sun"])))
magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
    logpretty=F)#, unlog='xy')
axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
    labels=c("O", "B", "A", "F", "G", "K", "M"))
points(log10(5777), 0)
points(log10(5777), 0, pch=20, cex=0.1)

points(DF$log_Teff[1], DF$log_L[1], pch=20, cex=0.5)
ZAMS <- which(10**DF$log_LH / 10**DF$log_L > 0.9)[1]
points(DF$log_Teff[ZAMS], DF$log_L[ZAMS], pch=1, cex=1)

dev.off()

