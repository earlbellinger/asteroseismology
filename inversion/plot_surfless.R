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

target.name <- 'CygA'
ref.mod <- 'CygAwball'
freqs <- get_freqs(target.name=target.name) 
CygA <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
CygA$target.name <- '16 Cyg A      '

target.name <- 'CygB'
ref.mod <- 'CygBwball'
freqs <- get_freqs(target.name=target.name) 
CygB <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
CygB$target.name <- '16 Cyg B      '

plot_surfless <- function(model, xlim=c(1400, 3700), ylim=c(-13, 8),
        legend.spot="topright", ylabs=T, ..., mar=utils.mar, 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", thin=F) {
    
    par(mgp=mgp-c(0.5, 0, 0))
    #par(mar=mar+c(-1, -0.5, 0, 0))
    
    nus <- model$nus
    diffs <- nus$nu.y - nus$nu.x
    nu <- nus$nu.x / model$nu_ac
    inertia <- nus$E 
    Xp <- matrix(c(nu**-1, nu**3) / (inertia * nus$dnu), ncol=2)
    a. <- ginv( Xp ) %*% ( diffs / nus$dnu )
    F_surf <- ( a.[[1]]*nu**-1 + a.[[2]]*nu**3 ) / (inertia)
    sigma.dist <- ((nus$nu.x+F_surf-nus$nu.y)/nus$dnu)
    model$nus <- cbind(model$nus, 
        data.frame(surfless = sigma.dist,
                   surf.unc = abs(sigma.dist) * nus$dnu/nus$nu.y))
    
    plot(NA, axes=F, #log='y',
        ylim=ylim,
        xlim=xlim,
        ylab="",
        xlab="")
    rect(xlim[1]*0.5, -3, xlim[2]*1.5, 3, 
        col=adjustcolor("gray", alpha.f=0.5), border=NA)
    rect(xlim[1]*0.5, -2, xlim[2]*1.5, 2, 
        col=adjustcolor("gray", alpha.f=0.9), border=NA)
    abline(v=0, lty=2)
    
    #with(model$nus, arrows(nu.y-dnu, surfless, nu.y+dnu, surfless, 
    #    code=3, angle=90, length=0.01, lwd=0.5))
    #with(model$nus, arrows(nu.y, surfless-surf.unc, nu.y, surfless+surf.unc, 
    #    code=3, angle=90, length=0.01, lwd=0.5))
    
    #col.pal <- c("#323031", "#DB4D48", "#F29559", blue)
    col.pal <- c(blue, "#323031", "#F29559", "#DB4D48") #blue)
    ell.pch <- c(22,24,23,21) #c(22,23,21,24)
    points(model$nus$nu.y, model$nus$surfless, cex=1, 
        col=1, lwd=1, 
        bg=col.pal[model$nus$l+1], 
        pch=ell.pch[model$nus$l+1])
    
    if (!is.null(legend.spot)) {
        legend(legend.spot, inset=c(0.04, 0.14), 
            pch=ell.pch, pt.bg=col.pal, bg="white",
            col=1, cex=text.cex, #inset=c(0, 0.2),
            legend=c("\u2113 = 0",
                     "\u2113 = 1",
                     "\u2113 = 2",
                     "\u2113 = 3"))
    } else {
        text(1500, -0.05,  '\u2264 3\u03C3', cex=text.cex)
        text(1500, 4.70,   '\u003E 3\u03C3', cex=text.cex)
        text(1500, -4.75,  '\u003E 3\u03C3', cex=text.cex)
    }
    
    legend('bottom', cex=text.cex, bty='n', legend=model$target.name,
        inset=c(0, 0.02))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=T, mgp=mgp-c(2,0.3,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=ylabs, #usepar=1,
            las=T, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    title(xlab=bquote("Frequency"~nu/mu*Hz))
    if (ylabs) {
        par(mgp=mgp+c(0,0,0), xpd=NA)
        title(ylab=bquote("(Data" - "Model)" / "Uncertainty"))
        par(xpd=F)
    }
}


model.list <- list(CygA, CygB)
plot_surflesses <- function(model.list, ...) {
    par(mfrow=c(1,length(model.list)))
    for (ii in 1:length(model.list)) { 
        plot_surfless(model=model.list[[ii]], ylabs=ii==1,
            legend.spot=if(ii==1) "topright" else NULL, ...)
    }
}
make_plots(plot_surflesses, 'surfless', filepath=file.path('plots', 'echelle'),
    mar=c(2,0.5,1,0.5), oma=c(0,2.2,0,0), cex.paper=0.75,
    model.list=model.list, use.cairo=T)

#make_plots(plot_surfless, 
#    paste0('surfless-CygA'), 
#    filepath=file.path('plots', 'echelle'), 
#    model=CygA, xlim=c(1400, 3700), 
#    legend.spot='right', use.cairo=T)
#make_plots(plot_surfless, 
#    paste0('surfless-CygB'), 
#    filepath=file.path('plots', 'echelle'), 
#    model=CygB, xlim=c(1400, 3700), 
#    legend.spot=NULL, use.cairo=T)

