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
seis.DF <- data.table(read.table(file.path('..', 'forward', 'simulations.dat'), 
    header=1))
setkey(seis.DF, M, Y, Z, alpha, overshoot, diffusion)
keys <- key(seis.DF)

solar_vals <- read.table(
        file.path('..', 'inverse', 'perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

# select only stars having more than half their track with Teff < 7000
combos <- unique(seis.DF[,keys, with=0])
keeps <- c()
for (i in 1:nrow(combos)) {
    slice <- merge(seis.DF, combos[i,])
    keep <- sum(slice$Teff<7000) > nrow(slice)/2
    keeps <- append(keeps, rep(keep, nrow(slice)))
}
seis.DF <- seis.DF[keeps,]

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(Map(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]

# remove l=3 modes
seis.DF <- seis.DF[,-grep("3", names(seis.DF)), with=0]

## Make scatter and contour plots of X-R and C-D diagrams for all model vars
scatter_mesh('Teff', 'L', 'M', mesh=F, thin=F, short=F)#, log='xy')
scatter_mesh('Teff', 'L', 'radius', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0', 'dnu02', 'M', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0', 'dnu02', 'age', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0', 'dnu02', 'X_c', mesh=F, thin=F, short=F)
scatter_mesh('Dnu0', 'radius', 'M', mesh=F, thin=F, short=F)#, log='xy')
scatter_mesh('M', 'radius', 'Dnu0', mesh=F, thin=F, short=F)#, log='xy')

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

cairo_pdf('plots/corr-spearman.pdf', width=6.97522, height=6.97522, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1,0,0,0), cex=1, cex.lab=0.8)
a <- corrplot(M, diag=1, type='lower', order="FPC", 
    p.mat=p.values, pch.cex=0.8, cl.cex=0.8, 
    tl.cex=0.8, tl.col=rgb(1,1,1,0), tl.srt=90, 
    sig.level=sig_level/length(M))
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
cols <- ifelse(grepl('M|Y|Z|alpha|diffusion|overshoot', names(seis.DF)[pos]), 
    '#800080', 'black')
text(1:ncol(seis.DF)-0.3, (ncol(seis.DF)+1):2-0.4, col=cols, pos=4, srt=90,
    as.expression(seis.labs[colnames(a)]), cex=0.8)
text(0.8, ncol(seis.DF):1-0.1, col=cols, pos=2, 
    as.expression(seis.labs[colnames(a)]), cex=0.8)
mtext(expression("  Spearman correlation coefficient"~rho), outer=1,
    side=1, cex=0.8)
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

cairo_pdf('plots/corr-pca.pdf', width=6.97522, height=4.17309, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1, 0, 0, 0), 
    cex=1, cex.lab=0.8)
a <- corrplot(M2, pch.cex=0.8, cl.cex=0.8, tl.cex=2, cl.pos='b', 
    tl.col=rgb(1,1,1,0), tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|^Y$|Z|alpha|overshoot|diffusion', 
    names(seis.DF)[1:15]), '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2, cex=0.8)
text(1:15, ncol(pcs)+.4, 
    as.expression(Map(function(x) bquote(.(x)[""]), 
    seis.labs[names(seis.DF)[1:15]])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3, cex=0.8)
mtext(expression("    Spearman correlation coefficient"~rho), outer=1,
    side=1, cex=0.8)
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

cairo_pdf('plots/corr-pca-observable.pdf', width=6.97522, height=4.17309, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1, 0, 0, 0), 
    cex=0.8, cex.lab=0.8)
a <- corrplot(M2, cl.pos='b', pch.cex=0.8, cl.cex=1, tl.cex=2,
    tl.col=rgb(1,1,1,0), tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|^Y$|Z|alpha|overshoot|diffusion', 
    names(seis.DF)[-1:-15]), '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2, cex=1) 
text(1:ncol(M2)+0.25, ncol(pcs)+.75, srt=45, 
    as.expression(Map(function(x) bquote(.(x)[""]), 
    seis.labs[names(seis.DF)[-1:-15]])),
    #as.expression(seis.labs[1:9]), 
    col=colors, pos=3, cex=1)
mtext(expression("    Spearman correlation coefficient"~rho), outer=1,
    side=1, cex=0.8)
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




pca <- prcomp(seis.DF[,-1:-15, with=0], center=TRUE, scale.=TRUE)
seis.DF <- cbind(seis.DF, pca$x)
scatter_mesh('PC1', 'PC2', 'age', mesh=F, thin=F, short=F, make_pdf=F)

scatter_mesh('Z', 'age', 'PC1', mesh=F, thin=F, short=F, make_pdf=F, log='x',
    quartiles=T)
scatter_mesh('M', 'age', 'PC2', mesh=F, thin=F, short=F, make_pdf=F,
    quartiles=T)
scatter_mesh('M', 'Z', 'PC2', log='y', mesh=F, thin=F, short=F, make_pdf=F,
    quartiles=T)


#scatter_mesh('r_sep02_median', 'r_sep13_median', 'X_c', mesh=F, thin=F, 
#    short=F)
#scatter_mesh('r_sep02_median', 'r_avg01_median', 'age', mesh=F, thin=F, 
#    short=F)
#scatter_mesh('r_sep02_median', 'r_avg01_median', 'M', mesh=F, thin=F, short=F)

#for (Z in names(seis.DF)[1:9]) {
#    scatter_mesh('Teff', 'L', Z)
#    scatter_mesh('Dnu0_median', 'dnu02_median', Z)
#}
#scatter_mesh('Dnu0_median', 'L', 'M')
#scatter_mesh('Dnu0_median', 'log_g', 'age')

## see the trends within individual tracks
combos <- combos[combos$diffusion > 0]
corrs <- Map( function(i) {
    print(i)
    track <- merge(seis.DF, combos[i,])[, c(-1:-6, -9), with=0]
    cor(track, method='s')
}, 1:nrow(combos) )



cairo_pdf('plots/corr-tracks.pdf', width=6.97522, height=6.97522, 
    family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1,0,0,0), cex=1, cex.lab=0.8)
a <- corrplot(apply(simplify2array(corrs), 1:2, mean), 
    diag=1, type='lower', order='FPC',
    pch.cex=0.8, cl.cex=0.8, tl.cex=0.8, tl.col=rgb(1,1,1,0), tl.srt=90)
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
text(1:ncol(a)-0.3, (ncol(a)+1):2-0.4, pos=4, srt=90,
    as.expression(seis.labs[colnames(a)]), cex=0.8)
text(0.7, ncol(a):1-0.1, pos=2, 
    as.expression(seis.labs[colnames(a)]), cex=0.8)
mtext(expression("  Spearman correlation coefficient"~rho), outer=1,
    side=1, cex=0.8)
dev.off()

