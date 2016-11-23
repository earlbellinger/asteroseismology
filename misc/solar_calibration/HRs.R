source('../../scripts/utils.R')

options(scipen=5)

filename <- file.path('diffusion', 'LOGS_MS', 'history.data')
DF1 <- read.table(filename, header=1, skip=5)

# clip PMS
decreasing_L <- which(diff(DF1$log_L) < 0 & DF1$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF1 <- DF1[-1:-pms,]
}

print("Temperature difference")
with(DF1, cat(paste(10**log_Teff[nrow(DF1)] - 10**log_Teff[1], "\n")))

filename <- file.path('no_diffusion', 'LOGS_MS', 'history.data')
DF2 <- read.table(filename, header=1, skip=5)

# clip PMS
decreasing_L <- which(diff(DF2$log_L) < 0 & DF2$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF2 <- DF2[-1:-pms,]
}

with(DF2, cat(paste(10**log_Teff[nrow(DF2)] - 10**log_Teff[1], "\n")))

plot_HRs <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font, 
        mar=c()) {
    par(family=font, cex=text.cex, cex.lab=text.cex, cex.axis=text.cex)
    plot(10**DF1$log_Teff, 10**DF1$log_L, 
        type='l', axes=F, tcl=0, 
        xlim=rev(10**range(DF1$log_Teff, DF2$log_Teff)), 
        ylim=10**range(DF1$log_L, DF2$log_L), 
        xlab=expression("Temperature"~T["eff"]/K), 
        ylab=expression("Luminosity"~L/L["☉"]))
        #xlab=expression("Temperature"~log(T["eff"]/K)), 
        #ylab=expression("Luminosity"~log(L/L["☉"])))
    lines(10**DF2$log_Teff, 10**DF2$log_L, lty=2)
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=1,
        logpretty=F, cex.axis=text.cex)
    axis(3, at=c(41000, 31000, 9500, 7240, 5920, 5300, 3850),
        labels=c("O", "B", "A", "F", "G", "K", "M"), cex=text.cex)
    points(5777.74, 1)
    points(5777.74, 1, pch=20, cex=0.1)
    #points(log10(5777), 0)
    #points(log10(5777), 0, pch=20, cex=0.1)
    
    legend("bottomleft", bty='n', lty=c(1,2), cex=text.cex,
        legend=c("Diffusion", "No diffusion"))
}

make_plots(plot_HRs, 'HRs')

#dev.off()

