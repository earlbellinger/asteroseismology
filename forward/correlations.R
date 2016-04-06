#### Mesh and scatterplot analysis of evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source(file.path('..', 'scripts', 'utils.R'))
#source('grid-meta.R')

library(magicaxis)
library(RColorBrewer)
library(parallel)
library(parallelMap)
library(data.table)
library(lattice)
library(corrplot)
library(ggplot2)
library(GGally)
library(scales)

col.pal <- adjustcolor(colorRampPalette(brewer.pal(11, "Spectral"))(21), 
    alpha.f=0.75)

## Load data
seis.DF <- data.table(read.table('simulations.dat', header=1))
setkey(seis.DF, M, Y, Z, alpha, overshoot, diffusion)
keys <- key(seis.DF)

solar_vals <- read.table(
        file.path('..', 'inverse', 'perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(Map(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]





metals <- unlist(Map(function(i) min(merge(seis.DF, combos[i,])$Fe.H), 
    1:nrow(combos)))

init.FeH <- unlist(Map(function(i) merge(seis.DF, combos[i,])[1,]$Fe.H, 
    1:nrow(combos)))

metals2 <- metals
metals2[metals2 < -8] <- -8

plot_diffusion <- function(Z=metals2, ..., 
        mar=c(3,3,1,6), text.cex=1, mgp=utils.mgp, thin=F) {
    
    thin.mar <- mar-c(0,0,0,2)
    if (thin) par(mar=thin.mar)
    
    col.pal <- brewer.pal(9, 'Spectral')
    #colorRampPalette(brewer.pal(10, 'Spectral'))(10)
    
    plot(10**-4+combos$diffusion, combos$M,
         axes=FALSE,
         log='x', xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[floor(9*normalize(cbind(metals2, init.FeH)))+1], 
         xlim=c(1e-04, 1e+2), 
         ylim=c(0.7, 1.6),
         xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
        mgp=mgp, cex.axis=text.cex, las=1)
    
    Z_max <- max(metals2, init.FeH)
    Z_min <- min(metals2, init.FeH)
    X_range <- 4000#10**diff(par()$usr)[1]
    Zlabs <- signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs[1] <- paste("<", Zlabs[1])
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext(expression("Final Surface Metallicity [Fe/H]"["TAMS"]), 
        4, line=ifelse(thin, 2.75, 4.5), cex=text.cex)
    
    points(1, 1, pch=1, cex=1)
    points(1, 1, pch=20, cex=0.1)
    
    title(xlab=get_label('diffusion'))
    title(ylab=get_label('M'))
}
make_plots(plot_diffusion, "Fe_H_M_D", thin.hack=1, mar=c(3, 3, 1, 6))

plot_diffusion2 <- function(Z=init.FeH, ..., 
        mar=c(3,3,1,6), text.cex=1, mgp=utils.mgp, thin=F) {

    thin.mar <- mar-c(0,0,0,2)
    if (thin) par(mar=thin.mar)

    col.pal <- brewer.pal(9, 'Spectral')
    #colorRampPalette(brewer.pal(10, 'Spectral'))(10)
    
    Z2 <- cbind(metals2, init.FeH)
    
    plot(10**-4+combos$diffusion, combos$M,
         axes=FALSE,
         log='x', xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[
             floor( 9 * ( (Z-min(Z2)) / (max(Z2)-min(Z2)) ) ) + 1
         ],
         xlim=c(1e-04, 1e+2), 
         ylim=c(0.7, 1.6),
         xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
        mgp=mgp, cex.axis=text.cex, las=1)
    
    Z_max <- max(metals2, init.FeH)
    Z_min <- min(metals2, init.FeH)
    X_range <- 4000#10**diff(par()$usr)[1]
    Zlabs <- signif(quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs[1] <- paste("<", Zlabs[1])
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext(expression("Initial Surface Metallicity [Fe/H]"["ZAMS"]), 
        4, line=ifelse(thin, 2.75, 4.5), cex=text.cex)
    
    points(1, 1, pch=1, cex=1)
    points(1, 1, pch=20, cex=0.1)
    
    title(xlab=get_label('diffusion'))
    title(ylab=get_label('M'))
}
make_plots(plot_diffusion2, "FeH0_M_D", thin.hack=1, mar=c(3, 3, 1, 6))

plot_diffusion3 <- function(..., text.cex=1, mgp=utils.mgp) {
    col.pal <- brewer.pal(9, 'Spectral')
    Z <- log10(10**-4 + combos$diffusion)
    
    plot(10**init.FeH, 10**metals, 
         xlab=expression("Initial Metallicity"~"[Fe/H]"["ZAMS"]),
         ylab=expression("Final Metallicity"~"[Fe/H]"["TAMS"]),
         axes=FALSE,
         #xaxs='i', yaxs='i', 
         pch=20, cex=0.5, 
         col=col.pal[
             floor( 9 * ( normalize(Z) ) ) + 1 
         ])#,
         #xlim=c(1e-04, 1e+2), 
         #ylim=c(0.7, 1.6),
         #xlab="", ylab="")
    
    magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(0,0,0,0),
        mgp=mgp, cex.axis=text.cex, las=1, unlog='xy')
    at <- seq(0, 4, 0.5)
    axis(1, at=at, labels=c(signif(log10(at), 2)), tick=F, 
        cex.axis=text.cex)
    axis(2, at=at, labels=c(signif(log10(at), 2)), tick=F, 
        cex.axis=text.cex, las=1)
    
    Z_max <- max(Z)
    Z_min <- min(Z)
    X_range <- diff(par()$usr)[1]
    Zlabs <- signif(10**quantile(seq(Z_min, Z_max, length=1000), 
                                 c(0, 0.25, 0.5, 0.75, 1)), 1)
    Zlabs <- paste0(" ", Zlabs)
    color.legend(par()$usr[2]+0.05*X_range, par()$usr[3], 
                 par()$usr[2]+0.10*X_range, par()$usr[4], 
                 Zlabs, 
                 col.pal[1:length(col.pal)], 
                 cex=text.cex, gradient='y', align='rb')
    mtext("Diffusion factor D", 4, line=4.5, cex=text.cex)
}
make_plots(plot_diffusion3, "D_FeH0_FeHf", thin.hack=1, mar=c(3, 3, 1, 6))




## Make scatter and contour plots of X-R and C-D diagrams for all model vars
scatter_mesh('Teff', 'L', 'M', mesh=F, thin=F, short=F)
scatter_mesh('Teff', 'L', 'radius', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0_median', 'dnu02_median', 'M', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0_median', 'dnu02_median', 'age', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0_median', 'dnu02_median', 'X_c', mesh=F, thin=F, short=F)
scatter_mesh('r_sep02_median', 'r_sep13_median', 'X_c', mesh=F, thin=F, short=F)
scatter_mesh('r_sep02_median', 'r_avg01_median', 'age', mesh=F, thin=F, short=F)
scatter_mesh('r_sep02_median', 'r_avg01_median', 'M', mesh=F, thin=F, short=F)

for (Z in names(seis.DF)[1:9]) {
    scatter_mesh('Teff', 'L', Z)
    scatter_mesh('Dnu0_median', 'dnu02_median', Z)
}
scatter_mesh('Dnu0_median', 'L', 'M')
scatter_mesh('Dnu0_median', 'log_g', 'age')

## Make correlation plot 
method <- 'spearman'
sig_level <- 10^-5

M <- cor(seis.DF, method=method)
cor.mtest <- function(mat, test) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], method=test)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
    return(p.mat)
}
p.values <- cor.mtest(seis.DF, method)

cairo_pdf('plots/corr-spearman.pdf', width=7, height=7, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1,0,0,0), cex=1, cex.lab=1)
a <- corrplot(M, diag=1, type='lower', order="FPC", 
    p.mat=p.values, pch.cex=0.95, cl.cex=1, 
    tl.cex=0.35, tl.col=rgb(1,1,1,0), tl.srt=90, 
    sig.level=sig_level/length(M))
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
cols <- ifelse(grepl('M|Y|Z|alpha|diffusion|overshoot', names(seis.DF)[pos]), 
    '#800080', 'black')
text(1:ncol(seis.DF)-0.3, (ncol(seis.DF)+1):2-0.4, col=cols, pos=4, srt=90,
    as.expression(seis.labs[colnames(a)]))
text(0.8, ncol(seis.DF):1-0.1, col=cols, pos=2, 
    as.expression(seis.labs[colnames(a)]))
mtext(expression("  Spearman correlation coefficient"~rho), outer=1,
    side=1)
dev.off()

cairo_pdf('plots/corr-spearman-slides.pdf', width=2*6.22665, height=2*4.1511, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0.5, 0, 0, 0), 
    cex=1.5, cex.lab=1.5)
a <- corrplot(M, diag=1, type='lower', order="FPC", 
    p.mat=p.values, pch.cex=0.95, 
    cl.cex=1, cl.pos='r', cl.ratio=0.25, cl.offset=0, 
    tl.col=rgb(1,1,1,0), tl.cex=0.35, tl.srt=90, 
    sig.level=sig_level/length(M))
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
cols <- ifelse(grepl('M|^Y$|Z|alpha|overshoot|diffusion', names(seis.DF)[pos]), 
    '#800080', 'black')
text(1:ncol(seis.DF)-0.3, (ncol(seis.DF)+1):2-0.4, col=cols, pos=4, srt=90,
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[colnames(a)])))
text(0.8, ncol(seis.DF):1-0.1, col=cols, pos=2, 
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[colnames(a)])))
mtext(expression("Spearman correlation coefficient"~rho~"       "), outer=1,
    side=4, line=-5, cex=1.75)
age_pos <- length(seis.DF)-which(colnames(a)=="age")+1
#segments(1, age_pos, x1=age_pos+2, lty=2, col='darkgray', lwd=2)
#segments(age_pos+2, age_pos, y1=1, lty=2, col='darkgray', lwd=2)
dev.off()



## Make PCA plot
pca <- prcomp(seis.DF[,-1:-15, with=0], center=TRUE, scale.=TRUE)
pcs <- pca$x[,cumsum(pca$sdev)/sum(pca$sdev)<0.99]
vars <- unlist(Map(function(x) signif(x, 3), (pca$sdev / sum(pca$sdev) * 100)))
pclabs <- as.expression(lapply(1:ncol(pcs), 
    function(ii) bquote("("*.(vars[ii])*"%)"~PC[.(ii)])))


M2 <- cor(pcs, seis.DF[,1:15, with=0], method=method)
cor.mtest2 <- function(mat1, mat2, test) {
    p.mat <- matrix(NA, ncol(mat1), ncol(mat2))
    for (i in 1:nrow(p.mat)) {
        for (j in 1:ncol(p.mat)) {
            p.mat[i, j] <- cor.test(mat1[, i], mat2[, j], method=test)$p.value
        }
    }
    return(p.mat)
}
p.values2 <- cor.mtest2(pcs, as.data.frame(seis.DF[,1:15, with=0]), method)

cairo_pdf('plots/corr-pca.pdf', width=8, height=4.17309, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1, 0, 0, 0), 
    cex=1, cex.lab=1)
a <- corrplot(M2, cl.cex=1, cl.pos='b', 
    tl.col=rgb(1,1,1,0), tl.cex = 0.3, tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|^Y$|Z|alpha|overshoot|diffusion', 
    names(seis.DF)[1:15]), '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2)
text(1:15, ncol(pcs)+.4, 
    as.expression(Map(function(x) bquote(.(x)[""]), 
    seis.labs[names(seis.DF)[1:15]])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3)
mtext(expression("    Spearman correlation coefficient"~rho), outer=1,
    side=1)
dev.off()

cairo_pdf('plots/corr-pca-slides.pdf', width=6.22665, height=4.1511, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0.5, 0, 0, 0), 
    cex=1, cex.lab=1)
a <- corrplot(M2, pch.cex=0.95,
    cl.cex=1, cl.pos='r', cl.ratio=0.25, cl.offset=0,
    tl.col=rgb(1,1,1,0), tl.cex=0.3, tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[1:9]), 
    '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2)
text(1:9, ncol(pcs)+.3, 
    as.expression(Map(function(x) bquote(.(x)[""]), 
    seis.labs[names(seis.DF)[1:15]])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3)
mtext(expression("Spearman correlation coefficient"~rho~"    "), line=-5.5,
    side=4, cex=1.1)
dev.off()



M2 <- cor(pcs, seis.DF[,-1:-15, with=0], method=method)
cor.mtest2 <- function(mat1, mat2, test) {
    p.mat <- matrix(NA, ncol(mat1), ncol(mat2))
    for (i in 1:nrow(p.mat)) {
        for (j in 1:ncol(p.mat)) {
            p.mat[i, j] <- cor.test(mat1[, i], mat2[, j], method=test)$p.value
        }
    }
    return(p.mat)
}
p.values2 <- cor.mtest2(pcs, as.data.frame(seis.DF[,-1:-15, with=0]), method)

cairo_pdf('plots/corr-pca-observable.pdf', width=6, height=4.17309, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1, 0, 0, 0), 
    cex=1, cex.lab=1)
a <- corrplot(M2, cl.cex=1, cl.pos='b', 
    tl.col=rgb(1,1,1,0), tl.cex = 0.3, tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|^Y$|Z|alpha|overshoot|diffusion', 
    names(seis.DF)[-1:-15]), '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2) 
text(1:ncol(M2)+0.5, ncol(pcs)+.75, srt=45, 
    as.expression(Map(function(x) bquote(.(x)[""]), 
    seis.labs[names(seis.DF)[-1:-15]])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3)
mtext(expression("    Spearman correlation coefficient"~rho), outer=1,
    side=1)
dev.off()

cairo_pdf('plots/corr-pca-slides.pdf', width=6.22665, height=4.1511, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0.5, 0, 0, 0), 
    cex=1, cex.lab=1)
a <- corrplot(M2, pch.cex=0.95,
    cl.cex=1, cl.pos='r', cl.ratio=0.25, cl.offset=0,
    tl.col=rgb(1,1,1,0), tl.cex=0.3, tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[1:9]), 
    '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2)
text(1:9, ncol(pcs)+.3, 
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[1:9])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3)
mtext(expression("Spearman correlation coefficient"~rho~"    "), line=-5.5,
    side=4, cex=1.1)
dev.off()


