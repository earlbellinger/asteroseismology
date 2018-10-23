source('../../scripts/utils.R')

models <- read.table('models_cep.dat', header=1)

data_out <- read.table('../learn-classical/tables-models_cep/LMC_CEP.dat', header=1)
data_in  <- read.table('data/cep_lmc_i.dat', header=1)

names(data_in)[1] <- 'Name'
data <- merge(data_in[,c(1:2)], data_out, by='Name')

plot_both <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times", tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.2, -0.5, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.66)
    
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0.4, 1.6), 
         ylim=c(-6, -2),
         xlab="",
         ylab="")
    
    with(models,
        points(log10(Period), I_M0, pch=20, cex=3, lwd=0.5, col='darkgray'))
    with(models,
        points(log10(Period), V_M0, pch=20, cex=3, lwd=0.5, col='gray'))
    
    #if (F) {
    with(data,
        points(log10(Period), I_M0, pch=20, cex=0.33, col=blue))
    with(data,
        points(log10(Period), V_M0, pch=20, cex=0.33, col=orange)) 
    #}
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.3, 0), unlog='x',
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, majorn=3)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), 
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, las=1, majorn=3)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression("Period"~P/"days"))
    par(mgp=mgp+c(0.8, 0, 0))
    title(ylab=expression("Magnitude"))
    
}

plot_both_W <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times", tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.2, -0.5, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.66)
    
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=c(0.4, 1.6), 
         ylim=c(-8, -3.5),
         xlab="",
         ylab="")
    
    with(models,
        points(log10(Period), W, pch=20, cex=3, lwd=0.5, col='gray'))
    with(data,
        points(log10(Period), W, pch=20, cex=0.33, col=red))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.3, 0), unlog='x',
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, majorn=3)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), 
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, las=1, majorn=3)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression("Period"~P/"days"))
    par(mgp=mgp+c(0.8, 0, 0))
    title(ylab=expression("Wesenheit Index"))
    
}

plot_both_HR <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font="Times", tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.2, -0.5, -0.3), mgp=mgp+c(0, 0.4, 0), lwd=1.66)
    
    plot(NA, axes=F, 
         xaxs='i', yaxs='i',
         xlim=log10(c(6500, 4000)), 
         ylim=c(2.5, 4.5),
         xlab="",
         ylab="")
    
    with(models,
        points(log10(Teff), logL, pch=20, cex=3, lwd=0.5, col='gray'))
    with(data,
        points(log10(Teff), logL, pch=20, cex=0.33, col=red))
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.3, 0), #unlog='x',
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, majorn=3)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.43, 0), 
        cex.axis=text.cex, lwd.ticks=par()$lwd,
        family=font, las=1, majorn=3)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    title(xlab=expression("Temperature"~log[10](T["eff"]/K)))
    par(mgp=mgp+c(0.8, 0, 0))
    title(ylab=expression("Luminosity"~log[10](L/L["sun"])))
    
}


make_plots(plot_both, 'PL')
make_plots(plot_both_W, 'PW')
make_plots(plot_both_HR, 'HR')

