#### Plot relative frequency differences for an array of stellar models 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source('models.R')
source('frequencies.R')
source('kernels.R')
source(file.path('..', 'scripts', 'utils.R'))
#parallelStartMulticore(9)

models <- get_model_list()

k.pair <- u_Y#u_Y#rho_Gamma1#rho_c2#rho_c2#rho_c2# # 
target.name <- 'CygAwball'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set) 
model.list <- list()
for (perturbed.model.name in perturbed.CygA.names) {
    model <- get_model(freqs=freqs, 
        model.name=perturbed.model.name, 
        target.name=target.name, 
        k.pair=k.pair, square.Ks=F) 
    model.list[[perturbed.model.name]] <- model
}
proxy <- get_model(freqs=freqs, model.name=target.name,
    target.name=target.name, k.pair=k.pair, square.Ks=F)

yoffset <- c(.5, 0, 0)

plot_freq_diffs <- function(model.list, freqs, short=F,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], 
            range((nu.x-nu.y-dnu)/nu.y, (nu.x-nu.y+dnu)/nu.y))
    }))
    
    xlim <- range(freqs$nu)
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], nu.x-nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=c(2000, 3500), #xlim, 
         ylim=c(-0.05, 0.05),#ylim+c(-0.01, 0.01),
         xlab=expression("Frequency"~nu/mu*Hz), 
         ylab="")
    
    abline(h=0, lty=2)#, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0 & model$nus$nu.y >= 2000,], 
            points(nu.y, (nu.x-nu.y)/nu.y, col=color[ii], pch=pch[ii], 
                cex=0.66))
    }
    
    if (F) {
    legend('right', pch=pch[ordering], col=color[ordering], inset=c(0.125,0),
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp+yoffset)
    title(ylab=expression("Actual"~delta*nu/nu))
}

make_plots(plot_freq_diffs, paste0('freq_diffs'), 
    filepath=file.path('plots', 'echelle'), model.list=model.list, freqs=freqs)




plot_function_differences <- function(model.list, short=F,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        range(model$d.f1.true) 
    }))
    
    xlim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        range(model$r) 
    }))
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    
    lty <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 1
    })
    
    plot(NA, axes=F, xaxs='i',
         xlim=xlim, #c(0, 1.01),#c(0.1,0.3), 
         ylim=ylim, #c(-0.15, 0.15),#ylim, 
         xlab=expression("Radius"~r/R), 
         ylab="")
    
    abline(h=0, lty=2)
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        lines(model$r, (model$f1-model$m2.f1)/model$m2.f1, lty=lty[ii], 
            col=adjustcolor(color[ii], alpha.f=0.75), 
            lwd=2)
    }
    
    if (F) {
    legend('bottomright', lty=lty, inset=c(0.01,0.01), col=color, 
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )
    )
    }
    
    #legend('bottomright', pch=c(NA,NA,NA,20,20,20), lty=c(2,1,3,NA,NA,NA),
    #    ncol=2, cex=0.8*text.cex, col=c(1,1,1, blue, 'black', red), bty='n',
    #    inset=c(0.02, 0.02),
    #    legend=c(expression(M==0.984), expression(M==1), expression(M==1.016),
    #             expression(R==0.98), expression(R==1), expression(R==1.02)))
    legend(0.67, -0.115, pch=NULL, lty=c(2,1,3),
        cex=0.8*text.cex, col=c(1,1,1), bty='n',
        inset=c(0.02, 0.02),
        legend=c(expression(M==0.984), expression(M==1), expression(M==1.016)))
    legend(0.47, -0.115, pch=c(1,1,1), lty=NULL,
        #ncol=2, 
        cex=0.8*text.cex, col=c(blue, 'black', red), bty='n',
        inset=c(0.05, 0.02),
        legend=c(expression(R==0.98), expression(R==1), expression(R==1.02)))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    #legend('bottomright', inset=c(0.01, 0.01), 
    #    lty=c(2, 1, 3, 1,1,1), 
    #    col=c('black', 'black', 'black', blue, 'black', red),
    #    legend=c(
    #        expression(M==0.984),
    #        expression(M==1),
    #        expression(M==1.016),
    #        expression(R==0.98),
    #        expression(R==1),
    #        expression(R==1.02)))
    
    #par(mgp=mgp+c(1.5, 0, 0))
    #f1.name <- model$f1.name
    #f1.exp <- model$f1.exp
    #title(ylab=bquote(.(f1.name)~'difference'~
    #    (.(f1.exp)['ref']-.(f1.exp))/.(f1.exp)))
    par(mgp=mgp+yoffset)
    title(ylab=expression(delta*c^2/c^2))
}

make_plots(plot_function_differences, paste0('func_diffs'), 
    filepath=file.path('plots', 'echelle'), model.list=model.list)





plot_freq_diffs_dimless <- function(model.list, freqs, short=F,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0,], 
            range((nu.x/sqrt(model$mass/model$radius^3) - 
                   nu.y/sqrt(proxy$mass/proxy$radius^3)) / 
                   (nu.y/sqrt(proxy$mass/proxy$radius^3))))
    }))
    
    xlim <- range(freqs$nu)
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], nu.x-nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=c(2000, 3500), #xlim, 
         ylim=ylim,#c(-0.008, 0),#ylim, #c(-0.05, 0.05),#ylim+c(-0.01, 0.01),
         xlab=expression("Frequency"~nu/mu*Hz), 
         ylab="")
    
    #abline(h=0, lty=2)#, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        #with(model$nus[model$nus$l == 0 & model$nus$nu.y >= 2000,], 
        #    points(nu.y, (nu.x-nu.y)/nu.y, col=color[ii], pch=pch[ii], 
        #        cex=0.66))
        #model <- model.list[[ii]]
        #with(model$nus[model$nus$l == 0,],
        #    segments(nu.y, nu.x-nu.y-dnu, nu.y, nu.x-nu.y+dnu))
        with(model$nus[model$nus$l == 0 & model$nus$nu.y >= 2000,], 
            points(nu.y, (nu.x/sqrt(model$mass/model$radius^3)-
                          nu.y/sqrt(proxy$mass/proxy$radius^3))/
                          (nu.y/sqrt(proxy$mass/proxy$radius^3)), 
                   col=color[ii], pch=pch[ii], cex=0.66))
    }
    
    if (F) {
    legend('right', pch=pch[ordering], col=color[ordering], inset=c(0.125,0),
        legend=c(
            expression(R == 0.98*","~M == 0.984),
            expression(R == 0.98*","~M == 1),
            expression(R == 0.98*","~M == 1.016),
            expression(R == 1*","~   M == 0.984),
            expression(R == 1*","~   M == 1),
            expression(R == 1*","~   M == 1.016),
            expression(R == 1.02*","~M == 0.984),
            expression(R == 1.02*","~M == 1),
            expression(R == 1.02*","~M == 1.016)
        )[ordering]
    )
    }
    
    legend('topright', pch=c(1,1,1,2,20,3), lty=NULL,
        ncol=2, cex=0.8*text.cex, col=c(blue, 'black', red, 1,1,1), bty='n',
        inset=c(0.02, 0.02),
        legend=c(expression(R==0.98), expression(R==1), expression(R==1.02),
                 expression(M==0.984), expression(M==1), expression(M==1.016)))
    
    #legend('topright', lty=NULL, pch=c(2,20,3),
    #    cex=0.8*text.cex, col=c(1,1,1), bty='n',
    #    inset=c(0.02, 0.02),
    #    legend=c(expression(M==0.984), expression(M==1), expression(M==1.016)))
    #legend('top', pch=c(20,20,20), lty=NULL,
    #    #ncol=2, 
    #    cex=0.8*text.cex, col=c(blue, 'black', red), bty='n',
    #    inset=c(0.05, 0.02),
    #    legend=c(expression(R==0.98), expression(R==1), expression(R==1.02)))
    
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    #par(mgp=mgp+c(1.5, 0, 0))
    #title(ylab=)
    par(mgp=mgp+yoffset)
    title(ylab=expression("Dimensionless"~delta*nu/nu))
}

make_plots(plot_freq_diffs_dimless, paste0('freq_diffs-dimless'), 
    filepath=file.path('plots', 'echelle'), model.list=model.list, freqs=freqs)





plot_freq_diffs_offset <- function(model.list, freqs, short=F,
        ..., text.cex=1, mgp=utils.mgp, mar=utils.mar) {
    
    ylim <- range(sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0 & model$nus$nu.y >= 2000,], 
            range(r.diff))
    }))
    
    xlim <- range(freqs$nu)
    
    color <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$R < 0.99) blue else if (model$R > 1.01) red else 'black'
    })
    pch <- sapply(1:length(model.list), function(ii) {
        model <- model.list[[ii]]
        if (model$M < 1) 2 else if (model$M > 1) 3 else 20
    })
    
    ordering <- rev(order(sapply(1:length(model.list), function(ii) {
        with(model.list[[ii]]$nus[model.list[[ii]]$nus$l == 0,], nu.x-nu.y)[1]
    })))
    
    plot(NA, axes=F, 
         xlim=c(2000, 3500), #xlim, 
         ylim=ylim,#c(-0.005, 0),#ylim,#ylim+c(-0.01, 0.01),
         xlab=expression("Frequency"~nu/mu*Hz), 
         ylab="")
    
    #abline(h=0, lty=2)#, col='#00000099')
    
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        with(model$nus[model$nus$l == 0 & model$nus$nu.y >= 2000,], 
            points(nu.y, r.diff, 
                col=color[ii], pch=pch[ii], 
                cex=0.66))
    }
    
    legend('topright', pch=c(1,1,1,2,20,3), lty=NULL,
        ncol=2, cex=0.8*text.cex, col=c(blue, 'black', red, 1,1,1), bty='n',
        inset=c(0.02, 0.02),
        legend=c(expression(R==0.98), expression(R==1), expression(R==1.02),
                 expression(M==0.984), expression(M==1), expression(M==1.016)))
    
    #legend('topright', lty=NULL, pch=c(2,20,3),
    #    cex=0.8*text.cex, col=c(1,1,1), bty='n',
    #    inset=c(0.02, 0.02),
    #    legend=c(expression(M==0.984), expression(M==1), expression(M==1.016)))
    #legend('top', pch=c(20,20,20), lty=NULL,
    #    #ncol=2, 
    #    cex=0.8*text.cex, col=c(blue, 'black', red), bty='n',
    #    inset=c(0.05, 0.02),
    #    legend=c(expression(R==0.98), expression(R==1), expression(R==1.02)))
    
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    
    #par(mgp=mgp+c(1.5, 0, 0))
    #title(ylab=)
    par(mgp=mgp+yoffset)
    title(ylab=expression(delta*nu/nu-"offset"))
}

make_plots(plot_freq_diffs_offset, paste0('freq_diffs-offset'), 
    filepath=file.path('plots', 'echelle'), model.list=model.list, freqs=freqs)

