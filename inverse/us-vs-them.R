#### Plot KAGES masses and ages against what we get with machine learning 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

col.pal <- c('black', 'black', blue, '#900090', red)

#############################################################
### Comparison with KAGES ###################################
#############################################################
kages <- read.table(file.path('data', 'kages.dat'), header=1)

data_dir <- file.path('learn_covs-simulations', 'kages')
ml <- do.call(plyr:::rbind.fill, Map(function(cov) {
        name <- sub('.dat', '', basename(cov))
        if (!name %in% kages$KIC) 
            return(NULL)
        DF <- read.table(cov, header=1)
        data.frame(Name = as.integer(name),
                   Age = median(DF$age),
                   dAgeL = median(DF$age) - quantile(DF$age, .16),
                   dAgeH = quantile(DF$age, .84) - median(DF$age),
                   Age_s = sqrt(var(DF$age)),
                   Mass = median(DF$M),
                   dMassL = median(DF$M) - quantile(DF$M, .16),
                   dMassH = quantile(DF$M, .84) - median(DF$M),
                   Mass_s = sqrt(var(DF$M)),
                   Radius = median(DF$radius),
                   dRadiusL = median(DF$radius) - quantile(DF$radius, .16),
                   dRadiusH = quantile(DF$radius, .84) - median(DF$radius),
                   Radius_s = sqrt(var(DF$radius)),
                   Luminosity = median(DF$L),
                   dLuminosityL = median(DF$L) - quantile(DF$L, .16),
                   dLuminosityH = quantile(DF$L, .84) - median(DF$L),
                   Luminosity_s = sqrt(var(DF$L)),
                   logg = median(DF$log_g),
                   dloggL = median(DF$log_g) - quantile(DF$log_g, .16),
                   dloggH = quantile(DF$log_g, .84) - median(DF$log_g),
                   logg_s = sqrt(var(DF$log_g)),
                   D = median(DF$diffusion),
                   dDL = median(DF$diffusion) - quantile(DF$diffusion, .16),
                   dDH = quantile(DF$diffusion, .84) - median(DF$diffusion),
                   D_s = sqrt(var(DF$diffusion)),
                   Z = median(DF$Z),
                   dZL = median(DF$Z) - quantile(DF$Z, .16),
                   dZH = quantile(DF$Z, .84) - median(DF$Z),
                   Z_s = sqrt(var(DF$Z)),
                   Y = median(DF$Y),
                   dYL = median(DF$Y) - quantile(DF$Y, .16),
                   dYH = quantile(DF$Y, .84) - median(DF$Y),
                   Y_s = sqrt(var(DF$Y)),
                   ov = median(DF$overshoot),
                   dovL = median(DF$overshoot) - quantile(DF$overshoot, .16),
                   dovH = quantile(DF$overshoot, .84) - median(DF$overshoot),
                   ov_s = sqrt(var(DF$overshoot))
                   #dAge = sqrt(var(DF$age))
                  )
    }, file.path(data_dir, list.files(data_dir)) ))

kages <- kages[kages$KIC %in% ml$Name,]
kages <- kages[order(kages$KIC),]
ml <- ml[order(ml$Name),]

plot_comparison <- function(qty, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    lims <- range(kages[[qty]]-kages[[low]], 
                  kages[[qty]]+kages[[high]],
                     ml[[qty]]-ml[[low]], 
                     ml[[qty]]+ml[[high]])
    
    #ml.sigma <- (ml[[qty]] * (ml[[high]]-ml[[low]])) / 2
    #kg.sigma <- (kages[[qty]] * (kages[[high]]-kages[[low]])) / 2
    sigma <- ml[[paste0(qty, '_s')]] + 
        ((kages[[qty]]+kages[[high]]) - (kages[[qty]]-kages[[low]])) / 2
    distance <- abs(ml[[qty]] - kages[[qty]]) / sigma
    #distance <- if (distance<2*ml.sigma) 0 else distance
    #dist2 <- abs(ml[[qty]] - kages[[qty]]) / kg.sigma
    #distance <- (dist1+dist2)/2
    #col.pal <- colorRampPalette(c(red, 'black'))(100)
    distance <- ifelse(distance > 5, 5, 1+floor(distance))
    #col.pal <- c('black', 'black', blue, 'darkred', red)
    
    name <- if (qty == 'Age') { 'age'
    } else if (qty == 'Mass') { 'M'
    } else if (qty == 'Luminosity') { 'L'
    } else if (qty == 'Radius') { 'radius'
    } else if (qty == 'logg') { 'log_g' }
    
    plot(NA, axes=F, ylab="", xlab="",
        xlim=lims, ylim=lims)
    abline(coef=c(0,1), lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(ml[[qty]]-ml[[low]], kages[[qty]], 
           ml[[qty]]+ml[[high]], kages[[qty]], 
        length=0, angle=90, code=3, col="darkgray")
    arrows(ml[[qty]], kages[[qty]]-kages[[low]], 
           ml[[qty]], kages[[qty]]+kages[[high]], 
        length=0, angle=90, code=3, col="darkgray")
    points(kages[[qty]] ~ ml[[qty]], pch=1, cex=0.5, 
        col=col.pal[distance])
        #[1+floor(99*dnorm(distance)/dnorm(0))])
    title(xlab=bquote("ML"~.(get_label(name))))
    title(ylab=bquote("KAGES"~.(get_label_nameless(name))))
}

for (qty in c("Age", "Mass", "Luminosity", "Radius", "logg")) {
    make_plots(plot_comparison, paste0('kages-', qty),
         filepath=file.path('plots', 'comparison'),
         qty=qty)
}

plot_diffusion <- function(..., text.cex=1, mgp=utils.mgp, font='Palatino') {
    distance <- abs(ml$D - 1) / ifelse(ml$D > 1, ml$dDL, ml$dDH)
    distance <- ifelse(distance > 5, 5, 1+floor(distance))
    
    plot(NA, axes=F, ylab="", xlab="", log='y', #yaxs='i', 
        ylim=range(ml$D - ml$dDL, ml$D + ml$dDH),
        xlim=range(ml$Mass - ml$dMassL, ml$Mass + ml$dMassH))
    
    abline(h=1, lty=2)
    arrows(ml$Mass, ml$D-ml$dDL, 
           ml$Mass, ml$D+ml$dDH, 
        length=0, angle=90, code=3, col="darkgray")
    arrows(ml$Mass-ml$dMassL, ml$D,
           ml$Mass+ml$dMassH, ml$D,
        length=0, angle=90, code=3, col="darkgray")
    points(ml$D ~ ml$Mass)
    #    col=colorRampPalette(c(blue, 'black', red))(1000)[
    #        floor(999*normalize(ml$Y)+1)])
    #brewer.pal(10, "Spectral")[floor(9*normalize(ml$ov)+1)])
    #col=col.pal[distance])
    
    title(xlab=get_label("M"))
    title(ylab=get_label("diffusion"))
    
    mdl <- deming(D~Mass, data=ml, xstd=Mass_s, ystd=D_s, conf=.68) 
    #lm(D ~ Mass, data=ml)
    #coefs <- summary(mdl)$coefficients
    abline(mdl, untf=T)
    
    #newMs <- seq(0.7, 1.6, 0.01)
    #lower <- mdl$ci[[1]] + mdl$ci[[2]] * newMs
    #upper <- mdl$ci[[3]] + mdl$ci[[4]] * newMs
    #lines(newMs, lower, lty=3)
    #lines(newMs, upper, lty=3)
    #prd <- predict(mdl, newdata=data.frame(Mass=newMs), 
    #    interval=c('confidence'), level=0.99, type='response')
    #prd[prd<0] <- 10**-10
    #lines(newMs, prd[,2], col='black', lty=3)
    #lines(newMs, prd[,3], col='black', lty=3)
    
    #legend("bottomleft", col=c(col.pal[2:5], 1, 1), pch=c(1,1,1,1,NA,NA), 
    #       cex=text.cex,
    #       lty=c(NA,NA,NA,NA,2,1), bty='n',
    #       legend=expression(""<2*sigma ~ "(consistent with D=1)",
    #           ""<3*sigma, ""<4*sigma, "">=4*sigma,
    #           D==1, D==(4.15%+-%0.55)-(2.7%+-%0.46)%*%M))
    legend("bottomleft", cex=text.cex, pch=c(3, NA, NA), #bty='n', 
           lty=c(NA,2,1), col=c("darkgray", "black", "black"),
           legend=c("KOIs", expression(D==1), 
        bquote(D==
            .(signif(mdl$coefficients[1],3)) -
            .(abs(signif(mdl$coefficients[2],3)))%*%M)))
            #(.(signif(coefs[[1]],3)) %+-% .(signif(coefs[[3]],3))) -
            #(.(abs(signif(coefs[[2]],3))) %+-% .(signif(coefs[[4]],3)))%*%M)))
    legend("bottomleft", cex=text.cex, pch=c(1,NA,NA), legend=c("","",""), 
        bty='n', lty=c(NA,NA,1))
    
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
}

make_plots(plot_diffusion, paste0('diffusion'),
     filepath=file.path('plots', 'comparison'))

#############################################################
### Comparison with Hare-and-Hound ##########################
#############################################################
basu <- read.table(file.path('data', 'basu.dat'), header=1)
data_dir <- file.path('learn_covs-simulations', 'basu')
ml <- do.call(plyr:::rbind.fill, Map(function(cov) {
        name <- sub('.dat', '', basename(cov))
        if (!name %in% basu$Model) 
            return(NULL)
        DF <- read.table(cov, header=1)
        data.frame(
            Name         = as.integer(name),
            Age          = median(DF$age),
            dAgeL        = median(DF$age) - quantile(DF$age, .16),
            dAgeH        = quantile(DF$age, .84) - median(DF$age),
            Age_s        = sqrt(var(DF$age)),
            Mass         = median(DF$M),
            dMassL       = median(DF$M) - quantile(DF$M, .16),
            dMassH       = quantile(DF$M, .84) - median(DF$M),
            Mass_s       = sqrt(var(DF$M)),
            Radius       = median(DF$radius),
            dRadiusL     = median(DF$radius) - quantile(DF$radius, .16),
            dRadiusH     = quantile(DF$radius, .84) - median(DF$radius),
            Radius_s     = sqrt(var(DF$radius)),
            Overshoot    = median(DF$overshoot),
            dOvershootL  = median(DF$overshoot) - quantile(DF$overshoot, .16),
            dOvershootH  = quantile(DF$overshoot, .84) - median(DF$overshoot),
            Overshoot_s  = sqrt(var(DF$overshoot)),
            Diffusion    = median(DF$diffusion),
            dDiffusionL  = median(DF$diffusion) - quantile(DF$diffusion, .16),
            dDiffusionH  = quantile(DF$diffusion, .84) - median(DF$diffusion),
            Diffusion_s  = sqrt(var(DF$diffusion))
        )
    }, file.path(data_dir, list.files(data_dir)) )
)

basu <- basu[basu$Model %in% ml$Name,]
basu <- basu[order(basu$Model),]
ml <- ml[order(ml$Name),]

plot_comparison <- function(qty, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    lims <- range(basu[[qty]], ml[[qty]]-ml[[low]], ml[[qty]]+ml[[high]])
    
    #ml.sigma <- (ml[[qty]] * (ml[[high]] - -ml[[low]])) / 2
    ml.sigma <- ml[[paste0(qty, '_s')]]
    distance <- abs(ml[[qty]] - basu[[qty]]) / ml.sigma
    distance <- ifelse(distance > 5, 5, 1+floor(distance))
    #distance <- if (distance<2*ml.sigma) 0 else distance
    #col.pal <- colorRampPalette(c('black', 'black', blue, 'darkred', red))(100)
    #col.pal <- c('black', 'black', blue, 'darkred', red)
    
    name <- if (qty == 'Age') { 'age'
     } else if (qty == 'Mass') { 'M'
     } else if (qty == 'Overshoot') { 'overshoot'
     } else if (qty == 'Radius') { 'radius'
     } else if (qty == 'Diffusion') { 'diffusion' }
    
    plot(NA, axes=F, ylab="", xlab="", xlim=lims, ylim=lims)
    abline(coef=c(0,1), lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(basu[[qty]], ml[[qty]]-ml[[low]], 
           basu[[qty]], ml[[qty]]+ml[[high]], 
        length=0, angle=90, code=3, col="darkgray")
    points(ml[[qty]] ~ basu[[qty]], pch=1, cex=0.75, 
        col=col.pal[distance])#[1+floor(99*dnorm(distance)/dnorm(0))])
    title(ylab=bquote("ML"~.(get_label_nameless(name))))
    title(xlab=bquote("Hare-and-Hound"~.(get_label(name))))
}

for (qty in c("Age", "Mass", "Overshoot", "Radius", "Diffusion")) {
    make_plots(plot_comparison, paste0('basu-', qty),
         filepath=file.path('plots', 'comparison'),
         qty=qty)
}

