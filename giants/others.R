
    #Dnu_plot
    #Dnus <- get_separations('Dnu', freqs, 0)
    #plot(Dnus$nus, Dnus$separations, axes=F, tcl=0,
    #    las=1,
    #    xlab=expression('Frequency'~nu/mu*Hz),
    #    ylab=expression('Large separation'~Delta*nu/mu*Hz))
    #magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    
    
    
    
    #echelle
    #zerone <- freqs[freqs$l %in% 0:1,]
    #nus <- zerone$nu
    #ell <- zerone$l
    #plot(c(nus %% large_sep, (nus %% large_sep) + large_sep), 
    #     c(nus, nus), 
    #     tck=0, axes=FALSE, 
    #     pch=20, 
    #     col=c(red, 'black')[ell+1], 
    #     xlab="", 
    #     ylab=expression("Frequency"~nu/mu*Hz),
    #     xlim=c(0, large_sep*2),
    #     ylim=freq_range)
    #magaxis(side=1:4, tcl=-0.25, labels=c(0,1,0,0), las=1)
    #title(xlab=expression((nu ~ mod ~ Delta * nu[0]) / mu * Hz))
    #ticks <- axTicks(1)
    #ticks <- ticks[ticks < large_sep]
    #ticks <- c(ticks, ticks+large_sep)
    #axis(1, tcl=0.25, tick=F, at=ticks, labels=ticks%%large_sep)
    #abline(v=large_sep, lty=3)
    #legend("bottomleft", pch=20, col=c(red, 'black'), cex=0.5,
    #       inset=c(0.02, 0.02), legend=c(paste0("\u2113=", 0:1)))
