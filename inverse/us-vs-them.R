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
    magaxis(side=2:4, family=font, tcl=0.25, labels=c(1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    magaxis(side=1, family=font, tcl=0.25, labels=1,
        cex.axis=text.cex, las=1, mgp=mgp-c(0, 0.1, 0))
    arrows(kages[[qty]]-kages[[low]],  ml[[qty]], 
           kages[[qty]]+kages[[high]], ml[[qty]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    arrows(kages[[qty]], ml[[qty]]-ml[[low]], 
           kages[[qty]], ml[[qty]]+ml[[high]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(ml[[qty]] ~ kages[[qty]], pch=1, cex=0.5)
    title(ylab=bquote("ML"~.(get_label_nameless(name))))
    par(mgp=mgp-c(0.2, 0, 0))
    label <- if (qty == 'logg') "Surface gravity log g" else get_label(name)
    title(xlab=bquote("KAGES"~.(label)))
}

plot_rel_diff <- function(qty, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    std <- paste0(qty, '_s')
    xlims <- range(kages[[qty]]-kages[[low]], 
                   kages[[qty]]+kages[[high]])
    
    kages.std <- (kages[[high]]+kages[[qty]]-(kages[[qty]]-kages[[low]]))/2
    
    dif <- kages[[qty]] - ml[[qty]]
    dif.unc <- sqrt( kages.std^2 + ml[[std]]^2 )
    
    rel.dif <- dif / kages[[qty]]
    rel.unc <- sqrt( (dif.unc/dif)^2 + kages.std^2 ) * abs(kages[[qty]])
    
    ylims <- range(rel.dif + rel.unc, rel.dif - rel.unc)
    #ylims[1] <- max(ylims[1], -1)
    #ylims[2] <- min(ylims[2], 1)
    
    name <- if (qty == 'Age') { 'age'
    } else if (qty == 'Mass') { 'M'
    } else if (qty == 'Luminosity') { 'L'
    } else if (qty == 'Radius') { 'radius'
    } else if (qty == 'logg') { 'log_g' }
    
    plot(NA, axes=F, ylab="", xlab="", xlim=xlims, ylim=ylims)
    #abline(coef=c(0,1), lty=2)
    abline(h=0, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(kages[[qty]]-kages[[low]],  rel.dif, 
           kages[[qty]]+kages[[high]], rel.dif, 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    arrows(kages[[qty]], rel.dif - rel.unc, 
           kages[[qty]], rel.dif + rel.unc, 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(rel.dif ~ kages[[qty]], pch=1, cex=0.5)
    title(xlab=bquote("KAGES"~.(get_label(name))))
    title(ylab=bquote(
        (.(seis.labs[[name]])["KAGES"] - .(seis.labs[[name]])["ML"]) / 
         .(seis.labs[[name]])["KAGES"]
    ))
}

plot_abs_diff <- function(qty, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    std <- paste0(qty, '_s')
    xlims <- range(kages[[qty]]-kages[[low]], 
                   kages[[qty]]+kages[[high]])
    
    kages.std <- (kages[[high]]+kages[[qty]]-(kages[[qty]]-kages[[low]]))/2
    dif <- kages[[qty]] - ml[[qty]]
    dif.unc <- sqrt( kages.std^2 + ml[[std]]^2 )
    ylims <- range(dif + dif.unc, dif - dif.unc)
    
    name <- if (qty == 'Age') { 'age'
    } else if (qty == 'Mass') { 'M'
    } else if (qty == 'Luminosity') { 'L'
    } else if (qty == 'Radius') { 'radius'
    } else if (qty == 'logg') { 'log_g' }
    
    plot(NA, axes=F, ylab="", xlab="", xlim=xlims, ylim=ylims)
    #abline(coef=c(0,1), lty=2)
    abline(h=0, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(kages[[qty]]-kages[[low]],  dif, 
           kages[[qty]]+kages[[high]], dif, 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    arrows(kages[[qty]], dif - dif.unc, 
           kages[[qty]], dif + dif.unc, 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(dif ~ kages[[qty]], pch=1, cex=0.5)
    title(xlab=bquote("KAGES"~.(get_label(name))))
    title(ylab=bquote(
        .(seis.labs[[name]])["KAGES"] - .(seis.labs[[name]])["ML"]
    ))
}

aspect <- 4.65014666667
for (qty in c("Age", "Mass", "Luminosity", "Radius", "logg")) {
    make_plots(plot_comparison, paste0('kages-', qty),
         filepath=file.path('plots', 'comparison', 'kages'), qty=qty, 
         paper_pdf_width=aspect,  slides_pdf_width=aspect,
         paper_pdf_height=aspect, slides_pdf_height=aspect)
    make_plots(plot_rel_diff, paste0('kages-', qty, '-rel'),
         filepath=file.path('plots', 'comparison', 'kages'), qty=qty, 
         paper_pdf_width=aspect,  slides_pdf_width=aspect,
         paper_pdf_height=aspect, slides_pdf_height=aspect)
    make_plots(plot_abs_diff, paste0('kages-', qty, '-diff'),
         filepath=file.path('plots', 'comparison', 'kages'), qty=qty, 
         paper_pdf_width=aspect,  slides_pdf_width=aspect,
         paper_pdf_height=aspect, slides_pdf_height=aspect)
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

make_plots(plot_diffusion, paste0('diffusion'), mar.paper=c(2.5, 3, 1, 2),
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
    magaxis(side=2:4, family=font, tcl=0.25, labels=c(1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    magaxis(side=1, family=font, tcl=0.25, labels=1,
        cex.axis=text.cex, las=1, mgp=mgp-c(0, 0.1, 0))
    arrows(other[[qty]], ml[[qty]]-ml[[low]], 
           other[[qty]], ml[[qty]]+ml[[high]], 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(ml[[qty]] ~ other[[qty]], pch=1, cex=0.5)
    title(ylab=bquote("Predicted"~.(get_label_nameless(name))))
    par(mgp=mgp-c(0.2, 0, 0))
    title(xlab=bquote("True"~.(get_label(name))))
}

plot_rel_diff <- function(qty, other=basu, ..., 
        text.cex=1, mgp=utils.mgp, font='Palatino') {
    low <- paste0('d', qty, 'L')
    high <- paste0('d', qty, 'H')
    xlims <- range(other[[qty]])
    
    name <- if (qty == 'Age') { 'age'
    } else if (qty == 'Mass') { 'M'
    } else if (qty == 'Luminosity') { 'L'
    } else if (qty == 'Radius') { 'radius'
    } else if (qty == 'logg') { 'log_g' }
    
    std <- paste0(qty, '_s')
    dif <- other[[qty]] - ml[[qty]]
    rel.dif <- dif / other[[qty]]
    rel.unc <- ((ml[[qty]]+ml[[high]])-(ml[[qty]]-ml[[low]]))/2
    ylims <- range(rel.dif + rel.unc, rel.dif - rel.unc)
    
    plot(NA, axes=F, ylab="", xlab="", xlim=xlims, ylim=ylims)
    #abline(coef=c(0,1), lty=2)
    abline(h=0, lty=2)
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
        cex.axis=text.cex, las=1, mgp=mgp)
    arrows(other[[qty]], rel.dif - rel.unc, 
           other[[qty]], rel.dif + rel.unc, 
        length=0.01, lwd=1.5, angle=90, code=3, col="darkgray")
    points(rel.dif ~ other[[qty]], pch=1, cex=0.5)
    title(xlab=bquote("True"~.(get_label(name))))
    title(ylab=bquote(
        (.(seis.labs[[name]])["True"] - .(seis.labs[[name]])["ML"]) / 
         .(seis.labs[[name]])["True"]
    ))
}

for (qty in c("Age", "Mass", "Radius")) {
    make_plots(plot_comparison, paste0('basu-', qty),
         filepath=file.path('plots', 'comparison', 'basu'), 
         mar.paper=c(2.5, 3, 1, 2),
         qty=qty)
    make_plots(plot_rel_diff, paste0('basu-', qty, '-rel'),
         filepath=file.path('plots', 'comparison', 'basu'), 
         mar.paper=c(2.5, 3, 1, 2),
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
         filepath=file.path('plots', 'comparison', 'hares'), 
         mar.paper=c(2.5, 3, 1, 2),
         qty=qty, other=hares)
}

