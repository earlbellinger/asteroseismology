#### Plots the diffusion factor as a function of M for Kepler stars
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(deming)

col.pal <- c('black', 'black', blue, '#900090', red)

data_dir <- file.path('learn', 'covs-simulations')
data_dir.1 <- file.path(data_dir, 'perturb')
data_dir.2 <- file.path(data_dir, 'kages')
data_dir.3 <- file.path(data_dir, 'legacy')

ml <- do.call(plyr:::rbind.fill, 
    Map(function(covs) {
        name <- sub('.dat', '', basename(covs))
        if (grepl("amp", name) || name=="Tagesstern") return(NULL)
        DF <- read.table(covs, header=1)
        results <- data.frame(Name=name)
        for (column in names(DF)) {
            vals <- DF[,column]
            result <- data.frame(mean(vals), sqrt(var(vals)), 
                                 median(vals) - quantile(vals, .16), 
                                 quantile(vals, .84) - median(vals))
            colnames(result) <- c(column, paste0(column, '_s'), 
                                  paste0('d', column, 'L'),
                                  paste0('d', column, 'H'))
            results <- cbind(results, result)
        }
        results
    }, c(file.path(data_dir.1, list.files(data_dir.1)),
         file.path(data_dir.2, list.files(data_dir.2)),
         file.path(data_dir.3, list.files(data_dir.3)))
    )
)
ml <- ml[!duplicated(ml$Name),]


plot_diffusion <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {
    distance <- abs(ml$diffusion - 1) / ifelse(
        ml$diffusion > 1, ml$ddiffusionL, ml$ddiffusionH)
    distance <- ifelse(distance > 5, 5, 1+floor(distance))
    
    plot(NA, axes=F, ylab="", xlab="", log='y', #yaxs='i', 
        ylim=range(ml$diffusion - ml$ddiffusionL, 
                   ml$diffusion + ml$ddiffusionH),
        xlim=range(ml$M - ml$dML, ml$M + ml$dMH))
    
    abline(h=1, lty=2)
    arrows(ml$M, ml$diffusion-ml$ddiffusionL, 
           ml$M, ml$diffusion+ml$ddiffusionH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#00000028")
    arrows(ml$M-ml$dML, ml$diffusion,
           ml$M+ml$dMH, ml$diffusion,
        length=0.01, lwd=1.5, angle=90, code=3, col="#00000028")
    points(ml$diffusion ~ ml$M, lwd=1.5)
    with(ml[ml$Name == 'Sun',], points(diffusion ~ M, cex=0.1, pch=20))
    
    title(xlab=get_label("M"))
    title(ylab=get_label("diffusion"))
    
    mdl <- deming(diffusion~M, data=ml, xstd=M_s, ystd=diffusion_s, conf=.50) 
    
    #abline(mdl, untf=T)
    
    slope <- mdl[[1]][[2]]
    slope.se <- sqrt(mdl$variance[[4]])#/sqrt(nrow(ml))
    print(paste("P-value:", pt(slope / slope.se, mdl$n-2)))
    
    intercept <- mdl[[1]][[1]]
    M.new <- seq(0.7, 1.6, 0.001)
    D.new <- slope * M.new + intercept
    D.new[D.new < 0] <- 10^-10
    lines(D.new ~ M.new)
    
    newMs <- seq(0.7, 1.6, 0.01)
    lower <- mdl$ci[[1]] + mdl$ci[[2]] * newMs
    upper <- mdl$ci[[3]] + mdl$ci[[4]] * newMs
    lower[lower<0] <- 10^-10
    lines(newMs, lower, lty=3)
    lines(newMs, upper, lty=3)
    
    sig.figs <- gumr(slope, slope.se)
    slope <- sig.figs$value
    slope.se <- sig.figs$uncert
    
    intercept.se <- sqrt(mdl$variance[[1]])
    sig.figs <- gumr(intercept, intercept.se)
    intercept <- sig.figs$value
    intercept.se <- sig.figs$uncert
    
    cat(paste("\\text{D} = (", 
        signif(intercept, 3), "\\pm", signif(intercept.se, 3), 
        ")", ifelse(slope>=0, "+", "-"), "(", 
        abs(signif(slope, 3)), "\\pm", signif(slope.se, 3), 
        ") \\cdot \\text{M}/\\text{M}_\\odot\n"))
    
    legend.slope <- gumr(slope, slope.se)$value
    legend.inter <- gumr(intercept, intercept.se)$value
    legend("bottomleft", cex=text.cex, pch=c(3, NA, NA, NA), bty='n', 
           inset=c(0.02, 0.01), 
           lty=c(NA,2,1,3), col=c("darkgray", "black", "black", 'black'), 
           legend=c(
        "KOIs",
        expression(D==1), 
        bquote(D==
            .(signif(intercept,3)) -
            .(abs(signif(slope,3)))%*%M),
        "50% Confidence interval"))
    legend("bottomleft", inset=c(0.02, 0.01), 
        cex=text.cex, pch=c(1,NA,NA,NA), legend=c("","","",""), 
        lty=c(NA,NA,1,NA), bty='n')
    
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
}

make_plots(plot_diffusion, paste0('diffusion-big'), mar.paper=c(2.5, 3, 1, 1),
     filepath=file.path('plots', 'comparison'))

plot_Y_Z <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {
    plot(NA, axes=F, ylab="", xlab="", #log='y', #yaxs='i', 
        xlim=range(ml$Z - ml$dZL, ml$Z + ml$dZH),
        ylim=range(ml$Y - ml$dYL, ml$Y + ml$dYH))
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    title(ylab=get_label("Y"))
    title(xlab=get_label("Z"))
    arrows(ml$Z, ml$Y-ml$dYL, 
           ml$Z, ml$Y+ml$dYH, 
        length=0.01, lwd=1.5, angle=90, code=3, col="#00000028")
    arrows(ml$Z-ml$dZL, ml$Y,
           ml$Z+ml$dZH, ml$Y,
        length=0.01, lwd=1.5, angle=90, code=3, col="#00000028")
    points(ml$Y ~ ml$Z, lwd=1.5)
    with(ml[ml$Name == 'Sun',], points(Y ~ Z, cex=0.1, pch=20))
    
    
    mdl <- deming(Y~Z, data=ml, xstd=Z_s, ystd=Y_s, conf=.50) 
    
    slope <- mdl[[1]][[2]]
    slope.se <- sqrt(mdl$variance[[4]])#/sqrt(nrow(ml))
    print(paste("P-value:", pt(slope / slope.se, mdl$n-2)))
    
    intercept <- mdl[[1]][[1]]
    M.new <- seq(0.7, 1.6, 0.001)
    D.new <- slope * M.new + intercept
    D.new[D.new < 0] <- 10^-10
    lines(D.new ~ M.new)
    
    newMs <- seq(0.7, 1.6, 0.01)
    lower <- mdl$ci[[1]] + mdl$ci[[2]] * newMs
    upper <- mdl$ci[[3]] + mdl$ci[[4]] * newMs
    lower[lower<0] <- 10^-10
    lines(newMs, lower, lty=3)
    lines(newMs, upper, lty=3)
    
    sig.figs <- gumr(slope, slope.se)
    slope <- sig.figs$value
    slope.se <- sig.figs$uncert
    
    intercept.se <- sqrt(mdl$variance[[1]])
    sig.figs <- gumr(intercept, intercept.se)
    intercept <- sig.figs$value
    intercept.se <- sig.figs$uncert
    
}

make_plots(plot_Y_Z, paste0('y_z-big'), mar.paper=c(2.5, 3, 1, 1),
     filepath=file.path('plots', 'comparison'))

