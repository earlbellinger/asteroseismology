#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R')) 
source('models.R') 
source('frequencies.R')
source('kernels.R')  
source('OLA_invert.R') 
models <- get_model_list()
perturb <- F
k.pair <- u_Y

target.name <- 'BiSON'
ref.mod <- 'diffusion'
freqs.1 <- get_freqs(target.name=target.name) 
m1 <- get_model(freqs=freqs.1, model.name=ref.mod, target.name=target.name, 
                k.pair=k.pair, square.Ks=T) 

kernels <- Map(function(width.r_f) {
    Map(function(r_0) {
        get_target_kernel(r_0=r_0, r_f=0, width.r_f=width.r_f, 
            K=model$k1, f.spl=model$cs.spl, targ.kern.type='mod_Gauss')
    }, r_0=seq(0.05, 0.45, 0.1))
}, width.r_f=c(0.01, 0.02, 0.03))

plot_target_functions <- function(legend.spot="topleft", ..., 
        text.cex=1, mgp=utils.mgp, font="Times", short=F) {
    
    #text.cex <- text.cex*1.5
    #par(mar=utils.mar+c(0, 0.05, 0, 0), cex.lab=text.cex, mgp=mgp+c(0.25,0,0))
    
    plot(NA, axes=F, yaxs='i', #xaxs='i', 
        xlim=c(0, .5),
        ylim=c(0, 120),#round(range(kernels), -2),
        xlab=bquote("Radius"~r/R),
        ylab="Target Kernel Function")
    #abline(h=0, lty=2)
    
    col.pal <- c("#D00000", "#FFBA08", "#1C3144", "#DB3A34")
    
    for (kernel.width in 1:length(kernels)) {
        for (kernel.r_0 in 1:length(kernels[[kernel.width]])) {
            lines(model$k1$x, kernels[[kernel.width]][[kernel.r_0]],
                lty=c(2,1,3)[kernel.width], 
                col=col.pal[kernel.width], 
                lwd=c(2.5,2.5,3)[kernel.width])
        }
    }
    
    legend(legend.spot, #bty='n', #lty=NA, pch=NA, 
        cex=text.cex, inset=c(0.02, 0.04),
        lty=c(3,1,2), col=col.pal, lwd=c(2,2,2.5),
        legend=as.expression(c(
            bquote(Delta["f"]==0.01),
            bquote(Delta["f"]==0.02),
            bquote(Delta["f"]==0.03))
        ))
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=short, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, #usepar=1,
            las=short, mgp=mgp+c(1,0,0),#+c(1,0.2,0), 
            family=font, cex.axis=text.cex)
    #par(mgp=mgp+c(1,0,0))
    #title(ylab=bquote( "Frequency"~nu/mu*Hz ))

}

make_plots(plot_target_functions, 
    paste0('target_functions-', target.name), 
    filepath=file.path('plots', 'echelle'), 
    legend.spot='topleft')
