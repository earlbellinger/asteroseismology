library(magicaxis)

options(scipen=5)
options(bitmapType='cairo')

models <- Map(function(directory) {
        DF <- read.table(file.path(directory, 'LOGS_MS', 'history.data'), 
            header=1, skip=5)
        decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
        if (any(decreasing_L)) {
            pms <- max(decreasing_L)
            DF <- DF[-1:-pms,]
        }
        DF
    }, directory=list.dirs('models', recursive=F))

png('kiel.png', width=1260.248*3/3, height=753.973*3/2, 
    type='cairo', family='Times', res=400)
par(mar=c(3, 4, 1, 1), family='Times', mgp=c(2, 0.25, 0))

plot(NA, axes=F, tcl=0, 
    xlim=c(7000, 3500), #rev(10**range(sapply(models, function(model) model$log_Teff))), 
    ylim=c(5, 4.1), #range(sapply(models, function(model) model$log_g)), 
    xlab=expression("Temperature"~T["eff"]/K), 
    ylab=expression("Surface gravity"~log~g))

rect(7240, 4, 5920, 6, col='#ffffbf', border=NA) # F
rect(5920, 4, 5300, 6, col='#fff4e8', border=NA) # G
rect(5300, 4, 3850, 6, col='#ffddb4', border=NA) # K
rect(3850, 4, 0,    6, col='#ffbd6f', border=NA) # M

#points(5777, 4.43812, cex=0.5)
#points(5777, 4.43812, pch=20, cex=0.1)

for (model_i in 1:length(models)) {
    model <- models[[model_i]]
    lines(10**model$log_Teff, model$log_g, lwd=2, col="#0571b0")
    points(10**model$log_Teff[nrow(model)], model$log_g[nrow(model)], 
        pch=20, cex=0.5)
}

magaxis(side=1:4, family='Times', tcl=0.25, labels=c(1,1,0,0), las=1,
    logpretty=F)#, unlog='xy')
#axis(3, at=c(41000, 31000, 9500, 7240, 5920, 5300, 3850),
#    labels=c("O", "B", "A", "F", "G", "K", "M"))

dev.off()

