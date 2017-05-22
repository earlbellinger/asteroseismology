source('../../scripts/utils.R')

options(scipen=5)

filename <- file.path('diffusion_best', 'LOGS_MS', 'history.data')
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

filename <- file.path('no_diffusion_best', 'LOGS_MS', 'history.data')
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
        xlab=expression("Effective Temperature"~T["eff"]/K), 
        ylab=expression("Luminosity"~L/L["Sun"]))
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


Mars.1 <- 10**DF1$log_Teff * sqrt( 10**DF1$log_R*696*10**6 / ( 2*1.524*149.598*10**9 ) )
Mars.2 <- 10**DF2$log_Teff * sqrt( 10**DF2$log_R*696*10**6 / ( 2*1.524*149.598*10**9 ) )
plot_BB <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font, 
        mar=c()) {
    par(family=font, cex=text.cex, cex.lab=text.cex, cex.axis=text.cex)
    plot(DF1$star_age/10**9, Mars.1, 
        type='l', axes=F, tcl=0, 
        xlim=c(0, 4.572),
        ylim=range(Mars.1, Mars.2),
        ylab=expression("Blackbody Temperature of Mars"~T/K), 
        xlab=expression("Mars Age"~tau/Gyr))
    lines(DF2$star_age/10**9, Mars.2, lty=2)
    abline(v=4.57-3.5, lty=3)
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=0,
        logpretty=F, cex.axis=text.cex)
    #axis(3, at=c(41000, 31000, 9500, 7240, 5920, 5300, 3850),
    #    labels=c("O", "B", "A", "F", "G", "K", "M"), cex=text.cex)
    #points(5777.74, 1)
    #points(5777.74, 1, pch=20, cex=0.1)
    #points(log10(5777), 0)
    #points(log10(5777), 0, pch=20, cex=0.1)
    
    legend("bottomright", bty='n', lty=c(1,2), cex=text.cex,
        legend=c("Solar model with diffusion", "No diffusion"))
}

#make_plots(plot_BB, 'Martian_T')


Earth.1 <- 10**DF1$log_Teff * sqrt( 10**DF1$log_R*696*10**6 / ( 2*149.598*10**9 ) )
Earth.2 <- 10**DF2$log_Teff * sqrt( 10**DF2$log_R*696*10**6 / ( 2*149.598*10**9 ) )
plot_BB <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font, 
        mar=c()) {
    par(family=font, cex=text.cex, cex.lab=text.cex, cex.axis=text.cex)
    plot(DF1$star_age/10**9, Earth.1, 
        type='l', axes=F, tcl=0, 
        xlim=c(0, 4.572),
        ylim=range(Earth.1, Earth.2, 288),
        ylab=expression("Blackbody Temperature of Earth"~T/K), 
        xlab=expression("Age"~tau/Gyr)) 
    lines(DF2$star_age/10**9, Earth.2, lty=2)
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), las=0,
        logpretty=F, cex.axis=text.cex)
    points(4.572, 288)
    points(4.572, 288, pch=3, cex=0.7)
    
    #points(4.572, tail(Earth.1, 1))
    points(4.572, tail(Earth.1, 1), pch=20, cex=0.7)
    #axis(3, at=c(41000, 31000, 9500, 7240, 5920, 5300, 3850),
    #    labels=c("O", "B", "A", "F", "G", "K", "M"), cex=text.cex)
    #points(5777.74, 1)
    #points(5777.74, 1, pch=20, cex=0.1)
    #points(log10(5777), 0)
    #points(log10(5777), 0, pch=20, cex=0.1)
    
    legend("bottomright", bty='n', lty=c(1,2), cex=text.cex,
        legend=c("Solar model with diffusion", "No diffusion"))
}

make_plots(plot_BB, 'Earth_T')

#dev.off()

