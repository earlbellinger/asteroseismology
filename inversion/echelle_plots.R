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
freqs.1 <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs.1, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

target.name <- 'CygA'
ref.mod <- 'CygAwball'
freqs.2 <- get_freqs(target.name=target.name) 
m2 <- get_model(freqs=freqs.2, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

target.name <- 'CygB'
ref.mod <- 'CygBwball'
freqs.3 <- get_freqs(target.name=target.name) 
m3 <- get_model(freqs=freqs.3, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

model.list <- list(m1, m2, m3)
freqs.list <- list(freqs.1, freqs.2, freqs.3)
nu_maxs <- c(3090, 2201, 2552)
xlims <- list( c(20, 155), c(30, 125), c(40, 140) )

plot_echelle <- function(model, freqs, xlim=NULL, nu_max=3090, 
        legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    
    text.cex <- text.cex*1.5
    par(mar=utils.mar+c(0, 0.05, 0, 0), cex.lab=text.cex, mgp=mgp+c(0.25,0,0))
    
    Delta.nu <- seismology(freqs, nu_max)$Dnu0
    #model$nu_max)$Dnu0
    
    if (is.null(xlim)) xlim <- c(0, Delta.nu)
    
    plot(NA, axes=F, xaxs='i', 
        xlim=xlim,
        ylim=c(1000, 4000),#range(model$nus$nu.x, model$nus$nu.y),
        xlab=bquote((nu~mod~.(round(Delta.nu,0)))/mu*Hz),
        ylab="")
    abline(h=nu_max, lwd=2, lty=3, col='gray')
    abline(v=Delta.nu, lwd=1, lty=2, col='black')
    
    col.pal <- c("#177E89", "#323031", "#FFBF3F", "#DB3A34")
    
    for (ii in 0:1) {
    points(model$nus$nu.x %% Delta.nu + Delta.nu*ii, model$nus$nu.x, 
        cex=0.8, lwd=1, 
        col=col.pal[model$nus$l+1],
        pch=c(0,5,1,2)[model$nus$l+1])
        #pch=c(0,1,2,5)[model$nus$l+1])
    
    keep <- with(model$nus, (nu.y+dnu)%%Delta.nu - (nu.y-dnu)%%Delta.nu) > 0
    with(model$nus[keep,], 
        arrows((nu.y-dnu)%%Delta.nu + Delta.nu*ii, nu.y, 
               (nu.y+dnu)%%Delta.nu + Delta.nu*ii, nu.y, 
            code=3, angle=90, length=0.01, lwd=0.5))
    #with(model$nus[keep,], 
    #    arrows(nu.y%%Delta.nu, nu.y-dnu, 
    #           nu.y%%Delta.nu, nu.y+dnu, 
    #        code=3, angle=90, length=0.001, lwd=0.5))
    
    points(model$nus$nu.y %% Delta.nu + Delta.nu*ii, model$nus$nu.y, cex=0.7,
        col=1, lwd=0.5, #col.pal[model$nus$l+1], 
        bg=col.pal[model$nus$l+1], 
        pch=c(22,23,21,24)[model$nus$l+1])
        #pch=c(22,21,24,23)[model$nus$l+1])
    }
    
    legend(legend.spot[1], legend.spot[2], bty='n', #lty=NA, pch=NA, 
        cex=text.cex, #inset=c(-0.5, 0),
        #title.adj=-0.2,
        legend=model$target.name)
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(0,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    axis(1, at=axTicks(1), tick=F, cex.axis=text.cex,
        labels=round(axTicks(1)%%round(Delta.nu,0),0))
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    par(mgp=mgp+c(1,0,0))
    title(ylab=bquote( "Frequency"~nu/mu*Hz ))

}

if (F) {
make_plots(plot_echelle, 
    paste0('echelle-', target.name), 
    filepath=file.path('plots', 'echelle'), 
    model=m1, freqs=freqs.1, xlim=c(20, 150),
    legend.spot='topleft')
make_plots(plot_echelle, 
    paste0('echelle-', target.name), 
    filepath=file.path('plots', 'echelle'), 
    model=m1, freqs=freqs, xlim=c(30, 120),
    legend.spot='topleft')
make_plots(plot_echelle, 
    paste0('echelle-', target.name), 
    filepath=file.path('plots', 'echelle'), 
    model=m1, freqs=freqs, xlim=c(40, 140),
    legend.spot='topleft')
}

plot_echelles <- function(model.list, freqs.list, 
        xlims=NULL, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    par(mfrow=c(1,length(model.list)))
    for (ii in 1:length(model.list)) {
        plot_echelle(model.list[[ii]], freqs.list[[ii]], 
            xlim=xlims[[ii]], nu_max=nu_maxs[ii], legend.spot=list(
                c(50, 1450), 
                c(58, 1450),
                c(70, 1450))[[ii]],
                #c("center", "bottom"), c("center", "bottom"))[[ii]],
            #legend.spot,
            text.cex=text.cex, mgp=mgp, font=font, short=short)
    }
}
make_plots(plot_echelles, 
    paste0('echelles'), 
    filepath=file.path('plots', 'echelle'), 
    model.list=model.list, freqs.list=freqs.list, xlims=xlims,
    legend.spot='bottom')

