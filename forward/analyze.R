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

col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(21)

## Load data
seis.DF <- data.table(read.table('simulations.dat', header=1))
setkey(seis.DF, M, Y, Z, alpha)
keys <- key(seis.DF)

solar_vals <- read.table(
        file.path('..', 'inverse', 'perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

cygA_stds <- sqrt(diag(var(read.table(
        file.path('..', 'inverse', 'perturb', '16CygA_perturb.dat'), 
    header=1))))

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(Map(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]

# Make inputs diagram
inputs_plot <- function(text.cex, ...) {
    H <- 1-combos$Y-combos$Z
    varmax <- max(H)
    varmin <- min(H)
    
    #print(ggpairs(data=combos, #font=font, #axisLabels="show", 
    #        columnLabels=sapply(names(seis.DF)[1:4], 
    #                            function (name) get_label(name))))
    
    print(splom(combos,
        cex=0.5, pch=20,
        col=col.pal[floor((H-varmin)/(varmax-varmin)*(length(col.pal)-1))+1],
        xlab=NULL, ylab=NULL, 
        axis.text.cex=3*text.cex/4, 
        axis.text.lineheight=0.0001,
        axis.line.tck=0.25,
        axis.half = 0,
        axis.check.overlap = 1,
        xaxs='n', yaxs='n',
        varname.cex=3*text.cex/4, 
        varnames=as.expression(seis.labs[1:4])))
        #sapply(names(seis.DF)[1:4], function (name) get_label(name))))
}
make_plots(inputs_plot, "inputs", filepath=file.path("plots", "mesh"))

## Make scatter and contour plots of H-R and C-D diagrams for all model vars
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

plot_corr <- function() {
    par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0,0,0,0), cex=1, cex.lab=1)
    a <- corrplot(M, diag=1, type='lower', order="FPC", 
        p.mat=p.values, pch.cex=0.95, cl.cex=1, 
        tl.cex=0.3, tl.col='white', tl.srt=90, sig.level=sig_level/length(M))
    pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
    cols <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[pos]), 
        '#800080', 'black')
    text(1:ncol(seis.DF)-0.3, (ncol(seis.DF)+1):2-0.4, col=cols, pos=4, srt=90,
        as.expression(seis.labs[colnames(a)]))
    text(0.8, ncol(seis.DF):1-0.1, col=cols, pos=2, 
        as.expression(seis.labs[colnames(a)]))
}

cairo_pdf('plots/corr-spearman.pdf', width=8, height=7, family='Palatino')
plot_corr()
dev.off()

png('plots/corr-spearman.png', width=2890.8, height=2529.45, family='Palatino', 
    res=400)
plot_corr()
dev.off()

## Make PCA plot
pca <- prcomp(seis.DF[,-1:-9, with=0], center=TRUE, scale.=TRUE)
pcs <- pca$x[,cumsum(pca$sdev)/sum(pca$sdev)<0.99]
vars <- formatC(pca$sdev / sum(pca$sdev) * 100, digits=4)
pclabs <- as.expression(lapply(1:ncol(pcs), 
    function(ii) bquote("("*.(vars[ii])*"%)"~PC[.(ii)])))
M2 <- cor(pcs, seis.DF[,1:9, with=0], method=method)
cor.mtest2 <- function(mat1, mat2, test) {
    p.mat <- matrix(NA, ncol(mat1), ncol(mat2))
    for (i in 1:nrow(p.mat)) {
        for (j in 1:ncol(p.mat)) {
            p.mat[i, j] <- cor.test(mat1[, i], mat2[, j], method=test)$p.value
        }
    }
    return(p.mat)
}
p.values2 <- cor.mtest2(pcs, as.data.frame(seis.DF[,1:9, with=0]), method)

plot_PCA <- function() {
    par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0,0,0,0), cex=1, cex.lab=1)
    a <- corrplot(M2, tl.col='white', tl.cex = 0.3, cl.cex=1, tl.srt=90,
        cl.pos='b', p.mat = p.values2, sig.level = sig_level/length(M2))
    colors <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[1:9]), 
        '#800080', 'black')
    text(0.6, ncol(pcs):1, pclabs, pos=2)
    text(1:9, ncol(pcs)+.4, as.expression(seis.labs[1:9]), col=colors, pos=3)
}

cairo_pdf('plots/corr-pca.pdf', width=8, height=7, family='Palatino')
plot_PCA()
dev.off()

png('plots/corr-pca.png', width=2890.8, height=2529.45, family='Palatino', 
    res=400)
plot_PCA()
dev.off()










