#### Plot KAGES masses and ages against what we get with machine learning 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(deming)

col.pal <- c('black', 'black', blue, '#900090', red)

#############################################################
### Comparison with KAGES ###################################
#############################################################
kages <- read.table(file.path('data', 'kages.dat'), header=1)

data_dir <- file.path('learn', 'covs-simulations', 'kages')
ml <- do.call(plyr:::rbind.fill, Map(function(covs) {
        name <- sub('.dat', '', basename(covs))
        if (!name %in% kages$KIC) 
            return(NULL)
        DF <- read.table(covs, header=1)
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
    arrows(kages[[qty]]-kages[[low]],  ml[[qty]], 
           kages[[qty]]+kages[[high]], ml[[qty]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    arrows(kages[[qty]], ml[[qty]]-ml[[low]], 
           kages[[qty]], ml[[qty]]+ml[[high]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    #arrows(ml[[qty]]-ml[[low]], kages[[qty]], 
    #       ml[[qty]]+ml[[high]], kages[[qty]], 
    #    length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    #arrows(ml[[qty]], kages[[qty]]-kages[[low]], 
    #       ml[[qty]], kages[[qty]]+kages[[high]], 
    #    length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(ml[[qty]] ~ kages[[qty]], pch=1, cex=0.5)
    title(xlab=bquote("KAGES"~.(get_label(name))))
    title(ylab=bquote("ML"~.(get_label_nameless(name))))
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
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    arrows(ml$Mass-ml$dMassL, ml$D,
           ml$Mass+ml$dMassH, ml$D,
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(ml$D ~ ml$Mass)
    
    title(xlab=get_label("M"))
    title(ylab=get_label("diffusion"))
    
    mdl <- deming(D~Mass, data=ml, xstd=Mass_s, ystd=D_s) 
    
    #abline(mdl, untf=T)
    
    slope <- mdl[[1]][[2]]
    slope.se <- sqrt(mdl$variance[[4]])#/sqrt(nrow(ml))
    print(paste("P-value:", pt(slope / slope.se, mdl$n-2)))
    
    intercept <- mdl[[1]][[1]]
    M.new <- seq(0.7, 1.6, 0.001)
    D.new <- slope * M.new + intercept
    D.new[D.new < 0] <- 10^-10
    lines(D.new ~ M.new)
    
    intercept.se <- sqrt(mdl$variance[[1]])
    
    cat(paste("\\text{D} = (", 
        signif(intercept, 3),      "\\pm", signif(intercept.se, 3), 
        ")", ifelse(slope>=0, "+", "-"), "(", 
        abs(signif(slope,     3)), "\\pm", signif(slope.se,     3), 
        ") \\cdot \\text{M}/\\text{M}_\\odot\n"))
    
    legend("bottomleft", cex=text.cex, pch=c(3, NA, NA), bty='n', 
           inset=c(0.02, 0.01), 
           lty=c(NA,2,1), col=c("darkgray", "black", "black"), 
           legend=c("KOIs", expression(D==1), 
        bquote(D==
            .(signif(mdl$coefficients[1],3)) -
            .(abs(signif(mdl$coefficients[2],3)))%*%M)))
    legend("bottomleft", inset=c(0.02, 0.01), 
        cex=text.cex, pch=c(1,NA,NA), legend=c("","",""), 
        lty=c(NA,NA,1), bty='n')
    
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
}

make_plots(plot_diffusion, paste0('diffusion'),
     filepath=file.path('plots', 'comparison'))

#######################################################################
### Comparison with Sarbani's Hare-and-Hound ##########################
#######################################################################
age_mass_radius <- function(covs) {
    name <- sub('.dat', '', basename(covs))
    if (!name %in% basu$Model) 
        return(NULL)
    DF <- read.table(covs, header=1)
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
        Radius_s     = sqrt(var(DF$radius))
    )
}

basu <- read.table(file.path('data', 'basu.dat'), header=1)
data_dir <- file.path('learn', 'covs-simulations', 'basu')
ml <- do.call(plyr:::rbind.fill, Map(age_mass_radius, 
    file.path(data_dir, list.files(data_dir))))

basu <- basu[basu$Model %in% ml$Name,]
basu <- basu[order(basu$Model),]
ml <- ml[order(ml$Name),]

plot_comparison <- function(qty, other=basu, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    lims <- range(other[[qty]], 
                  ml[[qty]]-ml[[low]], 
                  ml[[qty]]+ml[[high]])
    
    name <- if (qty == 'Age') { 'age'
     } else if (qty == 'Mass') { 'M'
     } else if (qty == 'Overshoot') { 'overshoot'
     } else if (qty == 'Radius') { 'radius'
     } else if (qty == 'Diffusion') { 'diffusion' }
     else { qty }
    
    plot(NA, axes=F, ylab="", xlab="", xlim=lims, ylim=lims)
    abline(coef=c(0,1), lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(other[[qty]], ml[[qty]]-ml[[low]], 
           other[[qty]], ml[[qty]]+ml[[high]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(ml[[qty]] ~ other[[qty]], pch=1, cex=0.5)
    title(ylab=bquote("Predicted"~.(get_label_nameless(name))))
    title(xlab=bquote("True"~.(get_label(name))))
}

for (qty in c("Age", "Mass", "Radius")) {
    make_plots(plot_comparison, paste0('basu-', qty),
         filepath=file.path('plots', 'comparison'),
         qty=qty)
}

######################################################################
### Comparison with SpaceInn Hare-and-Hound ##########################
######################################################################
hares <- read.table(file.path('data', 'hares.dat'), header=1)
data_dir <- file.path('learn', 'covs-simulations', 'hares')
ml <- do.call(plyr:::rbind.fill, Map(function(covs) {
        name <- sub('.dat', '', basename(covs))
        if (!name %in% hares$Model) 
            return(NULL)
        DF <- read.table(covs, header=1)
        results <- data.frame(Name=name)
        for (column in names(hares)[-1]) {
            vals <- DF[,column]
            result <- data.frame(median(vals), 
                                 median(vals) - quantile(vals, .16), 
                                 quantile(vals, .84) - median(vals))
            colnames(result) <- c(column, 
                                  paste0('d', column, 'L'),
                                  paste0('d', column, 'H'))
            results <- cbind(results, result)
        }
        results
    }, 
    file.path(data_dir, list.files(data_dir))))
hares <- hares[hares$Model %in% ml$Name,]
hares <- hares[order(hares$Model),]
ml <- ml[order(ml$Name),]
for (qty in names(hares)[-1]) {
    make_plots(plot_comparison, paste0('hares-', qty),
         filepath=file.path('plots', 'comparison'),
         qty=qty, other=hares)
}

