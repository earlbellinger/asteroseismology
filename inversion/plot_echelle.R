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
ref.mod <- 'ModelS'
freqs.1 <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs.1, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=T) 
m1$target.name <- 'Sun      '

target.name <- 'CygA'
ref.mod <- 'CygAwball'
freqs.2 <- get_freqs(target.name=target.name) 
m2 <- get_model(freqs=freqs.2, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
m2$target.name <- '16 Cyg A      '

target.name <- 'CygB'
ref.mod <- 'CygBwball'
freqs.3 <- get_freqs(target.name=target.name) 
m3 <- get_model(freqs=freqs.3, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
m3$target.name <- '16 Cyg B      '

model.list <- list(m1, m2, m3)
freqs.list <- list(freqs.1, freqs.2, freqs.3)
nu_maxs <- c(3090, 2201, 2552)
xlims <- list( c(20, 160), c(30, 125), c(40, 145) )

plot_echelle <- function(model, freqs, xlim=NULL, nu_max=3090, ylabs=T,
        legend.spot="topleft", label.ells=F, ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    
    #par(mar=c(2.5, 1, 1, 1))#, oma=if (ylabs) c(1, 3, 1, 1) else c(1, 1, 1, 1))
    
    text.cex <- text.cex*1.5
    #par(mar=utils.mar+c(0, 0.05, 0, 0), cex.lab=text.cex, mgp=mgp+c(0.25,0,0))
    par(cex.lab=text.cex, mgp=mgp+c(0.5,0,0))
    
    Delta.nu <- seismology(freqs=freqs, nu_max=nu_max)$Dnu0
    Dnu <- round(Delta.nu, 0)
    #model$nu_max)$Dnu0
    
    if (is.null(xlim)) xlim <- c(0, Delta.nu)
    
    plot(NA, axes=F, xaxs='i', 
        xlim=xlim,
        ylim=c(1000, 4000),
        xlab=bquote((nu~mod~.(Dnu))/mu*Hz),
        ylab="")
    #abline(h=nu_max, lwd=2, lty=3, col='gray')
    #abline(v=Delta.nu, lwd=2, lty=3, col='gray')
    abline(v=Delta.nu, lwd=1, lty=2, col='black')
    
    #col.pal <- c("#177E89", "#323031", "#FFBF3F", "#DB3A34")
    col.pal <- c(blue, "#323031", "#F29559", "#DB4D48") #blue)
    opn.pch <- c( 0, 2, 5, 1) #c( 0, 5, 1, 2)
    ell.pch <- c(22,24,23,21) #c(22,23,21,24)
    
    for (ii in 0:1) {
        points(model$nus$nu.x %% Delta.nu + Delta.nu*ii, model$nus$nu.x, 
            cex=if (model$nus$l==2) 1 else 0.8, lwd=1, 
            col=col.pal[model$nus$l+1],
            pch=opn.pch[model$nus$l+1])
            #pch=c(0,1,2,5)[model$nus$l+1])
        
        keep <- with(model$nus, (nu.y+dnu)%%Delta.nu-(nu.y-dnu)%%Delta.nu) > 0
        with(model$nus[keep,], 
            arrows((nu.y-dnu)%%Delta.nu + Delta.nu*ii, nu.y, 
                   (nu.y+dnu)%%Delta.nu + Delta.nu*ii, nu.y, 
                code=3, angle=90, length=0.01, lwd=0.5))
        #with(model$nus[keep,], 
        #    arrows(nu.y%%Delta.nu, nu.y-dnu, 
        #           nu.y%%Delta.nu, nu.y+dnu, 
        #        code=3, angle=90, length=0.001, lwd=0.5))
        
        points(model$nus$nu.y %% Delta.nu + Delta.nu*ii, model$nus$nu.y, 
            cex=if (model$nus$l==2) 1 else 0.7,
            col=1, lwd=0.5, #col.pal[model$nus$l+1], 
            bg=col.pal[model$nus$l+1], 
            pch=ell.pch[model$nus$l+1])
            #pch=c(22,21,24,23)[model$nus$l+1])
    }
    
    #if (F) {
    legend('bottom', #legend.spot[1], legend.spot[2], 
        bty='n', #lty=NA, pch=NA, 
        cex=text.cex, #inset=c(-0.5, 0),
        #title.adj=-0.2,
        legend=model$target.name)
    #}
    
    if (ylabs) {
        legend('topleft', bty='n', cex=text.cex, inset=c(0.02, 0),
            legend=c('Model', 'Measurement'),
            pch=c(1, 19), col=c(1,1), bg=c(1,1))
    }
    
    if (label.ells) {
        text(x=43.5, y=3100, cex=text.cex, labels='\u2113 = 2')
        text(x=76.5, y=2650, cex=text.cex, labels='0')
        text(x=108, y=3000, cex=text.cex, labels='3')
        text(x=141, y=2250, cex=text.cex, labels='1')
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(0,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    axis(1, at=axTicks(1), tick=F, cex.axis=text.cex,
        labels=round(axTicks(1)%%round(Delta.nu,0),0))
    magaxis(side=2, tcl=0.25, labels=ylabs, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    if (ylabs) {
        par(mgp=mgp+c(1.5,0,0), xpd=NA)
        title(ylab=bquote( "Frequency"~nu/mu*Hz ))
        par(xpd=F)
    }
}

if (F) {
make_plots(plot_echelle, 
    paste0('echelle-BiSON'), 
    filepath=file.path('plots', 'echelle'), 
    #paper_pdf_width=6.97522,
    model=m1, freqs=freqs.1, xlim=c(20, 155), nu_max=3090,
    legend.spot='topleft', slides=F, make_png=F, use.cairo=T) 
make_plots(plot_echelle, 
    paste0('echelle-CygA'), 
    filepath=file.path('plots', 'echelle'), 
    paper_pdf_width=6.97522*2/3,
    model=m2, freqs=freqs.2, xlim=c(30, 125), nu_max=2201, 
    legend.spot='topleft', ylabs=F, wide=F, tall=F, slides=F, make_png=F,
    use.cairo=T) 
make_plots(plot_echelle, 
    paste0('echelle-CygB'), 
    filepath=file.path('plots', 'echelle'), 
    paper_pdf_width=6.97522*2/3,
    model=m3, freqs=freqs.3, xlim=c(40, 140), nu_max=2552, 
    legend.spot='topleft', ylabs=F, wide=F, tall=F, slides=F, make_png=F,
    use.cairo=T)
}

plot_echelles <- function(model.list, freqs.list, 
        xlims=NULL, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    par(mfrow=c(1,length(model.list)))
    for (ii in c(2,3,1)) { 
        plot_echelle(model.list[[ii]], freqs.list[[ii]], 
            ylabs=ii==2, 
            label.ells=ii==1,
            xlim=xlims[[ii]], nu_max=nu_maxs[ii], legend.spot=list(
                c(50, 1450), 
                c(58, 1450),
                c(70, 1450))[[ii]],
            text.cex=text.cex, mgp=mgp, font=font, short=short, ...)
    }
}
make_plots(plot_echelles, 
    paste0('echelles'), 
    filepath=file.path('plots', 'echelle'), thin=F, tall=F, slides=F, 
    model.list=model.list, freqs.list=freqs.list, xlims=xlims, 
    legend.spot='bottom', use.cairo=T, mar=c(3,1,1,1), oma=c(0,4,0,0))









target.name <- 'CygA'
ref.mod <- 'CygADiffusion'
freqs.1 <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs.1, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=T) 
m1$target.name <- '6.9 Gyr      '

target.name <- 'CygA'
ref.mod <- 'CygAyoung'
freqs.2 <- get_freqs(target.name=target.name) 
m2 <- get_model(freqs=freqs.2, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
m2$target.name <- '6 Gyr      '

target.name <- 'CygA'
ref.mod <- 'CygAyounger'
freqs.3 <- get_freqs(target.name=target.name) 
m3 <- get_model(freqs=freqs.3, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F, match.nl=F) 
m3$target.name <- '5 Gyr      '

model.list <- list(m3, m1, m2)
freqs.list <- list(freqs.3, freqs.1, freqs.2)
nu_maxs <- c(2201, 2201, 2201)
xlims <- list( c(35, 141), c(2, 105), c(35, 140) )
plot_echelles <- function(model.list, freqs.list, 
        xlims=NULL, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    par(mfrow=c(1,length(model.list)))
    for (ii in c(2,3,1)) { 
        plot_echelle(model.list[[ii]], freqs.list[[ii]], 
            ylabs=ii==2, 
            label.ells=F,
            xlim=xlims[[ii]], nu_max=nu_maxs[ii], legend.spot=list(
                c(50, 1450), 
                c(58, 1450),
                c(70, 1450))[[ii]],
            text.cex=text.cex, mgp=mgp, font=font, short=short, ...)
    }
}
make_plots(plot_echelles, 
    paste0('echelles-young'), 
    filepath=file.path('plots', 'echelle'), thin=F, tall=F, slides=F, 
    model.list=model.list, freqs.list=freqs.list, xlims=xlims, 
    legend.spot='bottom', use.cairo=T, mar=c(3,1,1,1), oma=c(0,4,0,0))









r02.W <- get_separations('r_sep', get_freqs('CygAwball'), 0, nu_max=2201)
r02.D <- get_separations('r_sep', get_freqs('CygADiffusion'), 0, nu_max=2201)
r02.6 <- get_separations('r_sep', get_freqs('CygAyoung'), 0, nu_max=2201)
r02.5 <- get_separations('r_sep', get_freqs('CygAyounger'), 0, nu_max=2201)
r02.A <- get_separations('r_sep', get_freqs('CygA'), 0, nu_max=2201)

plot_r02 <- function(..., 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", short=F) {
    xlim <- range(r02.A$nus)
    plot(NA, axes=F, yaxs='i', xaxs='i',
        xlim=c(1500, 3000),
        ylim=c(0.02, 0.08),
        xlab=expression("Frequency"~nu/mu*Hz),
        ylab="")
    
    par(xpd=NA)
    
    nus <- seq(min(r02.A$nus), max(r02.A$nus), 1)
    lines(nus, splinefun(r02.A$nus, r02.A$separations)(nus), 
        col="#4d4d4d", lwd=2, lty=3)
    
    points(r02.W$nus[r02.W$nus>1590 & r02.W$nus<3000], 
           r02.W$separations[r02.W$nus>1590 & r02.W$nus<3000], pch=1,
        cex=text.cex, lwd=1.5)
    points(r02.D$nus[r02.D$nus>1590 & r02.D$nus<3000], 
           r02.D$separations[r02.D$nus>1590 & r02.D$nus<3000], pch=4, 
        col=orange, cex=text.cex, lwd=1.5)
    points(r02.6$nus[r02.6$nus>1590 & r02.6$nus<3000], 
           r02.6$separations[r02.6$nus>1590 & r02.6$nus<3000], pch=3, 
        col='#551A8B', cex=text.cex, lwd=1.5)
    points(r02.5$nus[r02.5$nus>1590 & r02.5$nus<3000], 
           r02.5$separations[r02.5$nus>1590 & r02.5$nus<3000], pch=2, 
        col=blue, cex=text.cex, lwd=1.5)
    points(r02.A$nus, r02.A$separations, pch=20, cex=text.cex, lwd=1.5)
    
    par(xpd=F)
    
    legend('bottomleft', pch=rev(c(20, c(1,4,3,2))), inset=c(-0.005, -0.045),
        col=rev(c(1,1, orange, '#551A8B', blue)), #lwd=1.5, lty=NA, 
        legend=rev(c("16 Cyg A", 
            expression("Model"~italic('GOE')~"(7.15 Gyr)"), 
            "Model 6.9 Gyr", 
            "Model 6 Gyr", 
            "Model 5 Gyr")), bty='n', cex=text.cex)
    
    #magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
    #        las=short, mgp=mgp-c(0.5,0.15,0), las=1, 
    #        family=font, cex.axis=text.cex)
    
    #magaxis(side=2, tcl=0.25, labels=1,
    #        las=short, mgp=mgp, las=1, 
    #        family=font, cex.axis=text.cex)
    
    magaxis(side=1, tcl=-0.25, labels=1,
            las=short, mgp=mgp, las=1, 
            family=font, cex.axis=text.cex)
    
    magaxis(side=2, tcl=-0.25, labels=1,
            las=short, mgp=mgp+c(0,0.25,0), las=1, 
            family=font, cex.axis=text.cex)
    
    
    par(mgp=mgp+c(0.3, 0, 0))
    title(ylab=expression("Frequency ratio"~r[0*","~2]))
}
make_plots(plot_r02, filename='r02-CygA', 
    #paper_pdf_height=4.17309*6/11,
    filepath=file.path('plots', 'echelle'), use.cairo=T)
