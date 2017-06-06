#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R')) 
source(file.path('..', 'scripts', 'seismology.R')) 
source('models.R') 
source('frequencies.R') 
models <- get_model_list()
perturb <- F
k.pair <- NULL

target.name <- 'BiSON'
ref.mod <- 'diffusion'
freqs <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 

target.name <- 'CygA'
ref.mod <- 'CygAwball'
freqs <- get_freqs(target.name=target.name) 
m2 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 

target.name <- 'CygB'
ref.mod <- 'CygBwball'
freqs <- get_freqs(target.name=target.name) 
m3 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 

model.list <- list(m2, m3) #m1, #

for (ii in 1:length(model.list)) {
    nus <- model.list[[ii]]$nus 
    
    diffs <- nus$nu.y - nus$nu.x 
    
    nu <- nus$nu.x / model.list[[ii]]$nu_ac 
    inertia <- nus$E #Q# 
    Xp <- matrix(c(nu**-1, nu**3) / (inertia * nus$dnu), ncol=2) 
    a. <- ginv( Xp ) %*% ( diffs / nus$dnu ) 
    F_surf <- ( a.[[1]]*nu**-1 + a.[[2]]*nu**3 ) / (inertia) 
    
    surf.diff <- (nus$nu.x+F_surf-nus$nu.y)
    sigma.dist <- (surf.diff/nus$dnu) 
    
    model.list[[ii]]$nus <- cbind(model.list[[ii]]$nus, 
        data.frame(surf.diff = surf.diff, 
                   surfless = sigma.dist, 
                   surf.unc = abs(sigma.dist) * nus$dnu/nus$nu.y,
                   r_ts = model.list[[ii]]$r_ts))
}




plot_surfless <- function(model.list, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", thin=F) {
    
    #par(mar=utils.mar+c(0, 0.05, 0, 0))
    par(mgp=mgp-c(0.5,0,0))
    
    for (ii in 1:length(model.list)) {
        nus <- model.list[[ii]]$nus
        
        diffs <- nus$nu.y - nus$nu.x
        
        #nu <- model$nus$nu.x / model$nu_ac
        #r.diff <- ( nu - m1$nus$nu.y/sqrt(m1$M/R**3) )/nu #m1$nus$r.diff #
        #inertia <- m1$nu_ac * m1$nus$Q_norm #* m1$nus$d.r.diff
        nu <- nus$nu.x / model.list[[ii]]$nu_ac
        inertia <- nus$E #Q# 
        Xp <- matrix(c(nu**-1, nu**3) / (inertia * nus$dnu), ncol=2)
        a. <- ginv( Xp ) %*% ( diffs / nus$dnu ) #m1$nus$d.r.diff )
        F_surf <- ( a.[[1]]*nu**-1 + a.[[2]]*nu**3 ) / (inertia) #* nus$dnu)
        
        sigma.dist <- ((nus$nu.x+F_surf-nus$nu.y)/nus$dnu)
        
        model.list[[ii]]$nus <- cbind(model.list[[ii]]$nus, 
            data.frame(surfless = sigma.dist,
                       surf.unc = abs(sigma.dist) * nus$dnu/nus$nu.y))
            #surfless = ((nus$nu.x+F_surf-nus$nu.y)/nus$dnu)**2))
            #nus$nu.x+F_surf-nus$nu.y,
                       #chi2 = 
    }
    
    #xlim <- range(sapply(model.list, function(model) with(model$nus, 
    #    range(surfless - surf.unc, surfless + surf.unc))))
    #xlim <- c(-12, 5)
    #ylim <- range(sapply(model.list, function(model) range(model$nus$nu.y)))
    ylim <- c(1400, 3450)
    
    
    plot(NA, axes=F, #log='y',
        ylim=ylim,#ylim, #c(0.01, 100),
        xlim=c(-13, 5),# ylim, ##ylim,#ylim,# c(-3, 3),#
        ylab="",#bquote("Frequency"~nu/mu*Hz),#bquote( "Goodness of fit" ~ chi^2 ),
        xlab=bquote("(Data - Model)/Uncertainty")) #~ 
            #(nu["data"]-nu["model"])/sigma))
        #bquote( (nu["model"]-nu["star"]-F["surf"])/mu*Hz ))
    #rect(-3, xlim[1]*0.5, -2, xlim[2]*1.5, col="gray", border=NA)
    #rect(2, xlim[1]*0.5, 3, xlim[2]*1.5, col="gray", border=NA)
    rect(-3, ylim[1]*0.5, 3, ylim[2]*1.5, 
        col=adjustcolor("gray", alpha.f=0.5), border=NA)
    rect(-2, ylim[1]*0.5, 2, ylim[2]*1.5, 
        col=adjustcolor("gray", alpha.f=0.5), border=NA)
    abline(v=0, lty=2)
    #abline(v=2, lty=3)
    #abline(v=-2, lty=3, col=adjustcolor("black", alpha.f=1))
    #abline(v=3, lty=4)
    #abline(v=-3, lty=4)
    
    col.pal <- c("#DB4D48", "#F29559", "#B8B08D")
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        
        with(model$nus, 
            arrows(surfless-surf.unc, nu.y, 
                   surfless+surf.unc, nu.y, 
                code=3, angle=90, length=0.01, lwd=0.5))
        with(model$nus, 
            arrows(surfless, nu.y-surf.unc, 
                   surfless, nu.y+surf.unc, 
                code=3, angle=90, length=0.001, lwd=0.5))
    }
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        
        points(model$nus$surfless, model$nus$nu.y, cex=1, #lwd=0.5, 
            col=1, lwd=0.6, #col.pal[model$nus$l+1], 
            bg=col.pal[ii], 
            pch=c(22,23,21,24)[model$nus$l+1])
            #pch=c(22,21,24,23)[model$nus$l+1])
    }
    
    legend(legend.spot, inset=c(0.04,-0.03), pch=c(20, 20, 22,23,21,24),
        col=c(col.pal[1:length(model.list)], 1, 1, 1, 1), cex=text.cex,
        legend=as.expression(c("CygA", "CygB", 
            bquote("\u2113"==0),
            bquote("\u2113"==1),
            bquote("\u2113"==2),
            bquote("\u2113"==3))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=thin, mgp=mgp-c(2,0.3,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=thin, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(0.5,0,0))
    title(ylab=bquote("Frequency"~nu/mu*Hz))

}

make_plots(plot_surfless, 
    paste0('surfless'), 
    filepath=file.path('plots', 'echelle'), 
    model.list=model.list, #xlim=c(40, 140), 
    legend.spot='left')




plot_surfless <- function(model.list, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", thin=F) {
    
    par(mgp=mgp-c(0.5,0,0))
    
    ylim <- c(1400, 3450)
    plot(NA, axes=F,
        ylim=c(0,1),
        xlim=c(-13, 5),
        ylab="",
        xlab=bquote("(Data - Model)/Uncertainty")) 
    rect(-3, ylim[1]*0.5, 3, ylim[2]*1.5, 
        col=adjustcolor("gray", alpha.f=0.5), border=NA)
    rect(-2, ylim[1]*0.5, 2, ylim[2]*1.5, 
        col=adjustcolor("gray", alpha.f=0.5), border=NA)
    abline(v=0, lty=2)
    
    col.pal <- c("#DB4D48", "#F29559", "#B8B08D")
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        
        with(model$nus, 
            arrows(surfless-surf.unc, r_ts, 
                   surfless+surf.unc, r_ts, 
                code=3, angle=90, length=0.01, lwd=0.5))
        #with(model$nus, 
        #    arrows(surfless, nu.y-surf.unc, 
        #           surfless, nu.y+surf.unc, 
        #        code=3, angle=90, length=0.001, lwd=0.5))
    }
    for (ii in 1:length(model.list)) {
        model <- model.list[[ii]]
        
        points(model$nus$surfless, model$nus$r_ts, cex=1, 
            col=1, lwd=0.6, 
            bg=col.pal[ii], 
            pch=c(22,23,21,24)[model$nus$l+1])
    }
    
    legend(legend.spot, inset=c(0.04,-0.03), pch=c(20, 20, 22,23,21,24), 
        col=c(col.pal[1:length(model.list)], 1, 1, 1, 1), cex=text.cex, 
        legend=as.expression(c("CygA", "CygB", 
            bquote("\u2113"==0), 
            bquote("\u2113"==1), 
            bquote("\u2113"==2), 
            bquote("\u2113"==3))))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0), 
            las=thin, mgp=mgp-c(2,0.3,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=thin, mgp=mgp+c(1,0,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(0.5,0,0))
    title(ylab=bquote("Frequency"~nu/mu*Hz))

}

make_plots(plot_surfless, 
    paste0('surfless2'), 
    filepath=file.path('plots', 'echelle'), 
    model.list=model.list, 
    legend.spot='left')


    
