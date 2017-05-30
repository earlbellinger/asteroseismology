#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R')) 
source('models.R') 
source('frequencies.R') 
source('OLA_invert.R')
models <- get_model_list()
perturb <- F
k.pair <- NULL

target.name <- 'BiSON.MDI'
ref.mod <- 'hl.diff'
freqs <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

target.name <- 'BiSON'
ref.mod <- 'diffusion'
freqs <- get_freqs(target.name=target.name) 
m2 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

target.name <- 'BiSON'
ref.mod <- 'diffusion'
mode.set <- 'CygA'
freqs <- get_freqs(target.name=target.name, mode.set=mode.set) 
m3 <- get_model(freqs=freqs, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

model.list <- list(m1, m2, m3)

plot_turning_points <- function(model.list, legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", thin=F) {
    
    fct <- 100
    
    plot(NA, axes=F, frame.plot=T, #xaxs='i', yaxs='i', #frame.plot=T, xaxs='n',
        xlim=range(sapply(model.list, function(model) range(model$nus$nu.x))),
        ylim=fct**c(0,1),
        xlab=bquote("Frequency"~nu/mu*Hz),
        ylab=bquote( "Lower turning point"~r[t]/R["*"] ))
    
    col.pal <- c(NA, "#283845", "#F29559")
    for (ii in 1:length(model.list)) {
        
        model <- model.list[[ii]]
        
        nus <- cbind(model$nus, data.frame(r_ts=model$r_ts))
        
        for (ell in rev(sort(unique(nus$l)))) {
            nus.ell <- nus[nus$l==ell,]
            if (nrow(nus.ell) <= 1) next 
            if (any(nus.ell$r_ts < 0.999)) {
                nu.seq <- with(nus.ell, seq(min(nu.x), max(nu.x), 1))
                r_ts.seq <- splinefun(nus.ell$nu.x, fct**nus.ell$r_ts)(nu.seq)
                lines(nu.seq, r_ts.seq, col="gray", lty=1, lwd=1)
            }
        }
        for (ell in rev(sort(unique(nus$l)))) {
            nus.ell <- nus[nus$l==ell,]
            cex <- ifelse(nus.ell$l < 4 | nus.ell$r_ts <= 0.95, 0.66,  
                    0.66*log(nus.ell$r_ts, base=fct)/log(0.95, base=fct))
            points(nus.ell$nu.x, fct**nus.ell$r_ts, 
                cex=cex,
                col=1, lwd=0.5, 
                bg=col.pal[ii], 
                pch=if (ell < 4) c(22,23,21,24)[ell+1] else 1)
        }
        
    }
    at <- axTicks(2)
    at[1] <- 1
    labs <- signif(log(at, base=fct), 2)
    at <- signif(fct**labs, 2)
    minors <- signif(fct**seq(0, 1, 0.01), 2)
    #labs[1] <- 0
    axis(2, tick=T, at=at, cex.axis=text.cex, las=thin,
        labels=labs, tcl=0.25/2)
    axis(4, tick=T, at=at, cex.axis=text.cex, las=1,
        labels=F, tcl=0.25/2)
    axis(2, tick=T, at=minors, tcl=0.125/2, labels=F)
    axis(4, tick=T, at=minors, tcl=0.125/2, labels=F)
    
    magaxis(side=c(1,3), tcl=0.25/2, labels=c(1,0),
            las=thin, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    
    #magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
    #        las=short, mgp=mgp-c(0.5,0.15,0), 
    #        family=font, cex.axis=text.cex)
    #magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
    #        las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
    #        family=font, cex.axis=text.cex)

}
make_plots(plot_turning_points, 
    paste0('turning_points'), 
    filepath=file.path('plots', 'echelle'), 
    model.list=model.list, legend.spot='topleft')
