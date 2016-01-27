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
X <- 1-combos$Y-combos$Z
p = ggpairs(data=combos, axisLabels="show", #font=font, #
        columnLabels=sapply(names(seis.DF)[1:4], 
                            function (name) get_label_nameless(name)))
for (subplot in grep('points', p$plots)) {
    p$plots[[subplot]] <- sub("))", ", colour=X))", p$plots[[subplot]])
}
for (ii in 1:4) {
    for (jj in (ii+1):4) {
        if (jj > 4) break
        p <- putPlot(p, getPlot(p, ii, jj) + 
            theme(axis.ticks=element_blank(), axis.line=element_blank(), 
                  axis.text=element_blank(), panel.grid.major= element_blank()),
            ii, jj)
    }
}
for (col_j in 1:2) {
    zplot <- getPlot(p, 3, col_j)
    log_zplot <- zplot + scale_y_log10(limits=c(0.0001, 0.04))
    log_zplot$subtype <- 'logpoints'
    log_zplot$type <- 'logcontinuous'
    p <- putPlot(p, log_zplot, 3, col_j)
}
for (row_i in 3:4) {
    zplot <- getPlot(p, row_i, 3)
    log_zplot <- zplot + scale_x_log10(limits=c(0.0001, 0.04))
    log_zplot$subtype <- 'logpoints'
    log_zplot$type <- 'logcontinuous'
    p <- putPlot(p, log_zplot, row_i, 3)
}
inputs_plot <- function(text.cex, ...) {
    print(p, leftWidthProportion=0.3, bottomHeightProportion=0.3)
}
make_plots(inputs_plot, "inputs", filepath=file.path("plots", "inputs"))

## Make scatter and contour plots of X-R and C-D diagrams for all model vars
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

cairo_pdf('plots/corr-spearman.pdf', width=8, height=7, family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1,0,0,0), cex=1, cex.lab=1)
a <- corrplot(M, diag=1, type='lower', order="FPC", 
    p.mat=p.values, pch.cex=0.95, cl.cex=1, 
    tl.cex=0.3, tl.col=rgb(1,1,1,0), tl.srt=90, 
    sig.level=sig_level/length(M))
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(seis.DF))))
cols <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[pos]), 
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
cols <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[pos]), 
    '#800080', 'black')
text(1:ncol(seis.DF)-0.3, (ncol(seis.DF)+1):2-0.4, col=cols, pos=4, srt=90,
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[colnames(a)])))
text(0.8, ncol(seis.DF):1-0.1, col=cols, pos=2, 
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[colnames(a)])))
mtext(expression("Spearman correlation coefficient"~rho~"       "), outer=1,
    side=4, line=-5, cex=1.75)
age_pos <- length(seis.DF)-which(colnames(a)=="age")+1
segments(1, age_pos, x1=age_pos+2, lty=2, col='darkgray', lwd=2)
segments(age_pos+2, age_pos, y1=1, lty=2, col='darkgray', lwd=2)
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

cairo_pdf('plots/corr-pca.pdf', width=8, height=7, family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(1, 0, 0, 0), 
    cex=1, cex.lab=1)
a <- corrplot(M2, cl.cex=1, cl.pos='b', 
    tl.col=rgb(1,1,1,0), tl.cex = 0.3, tl.srt=90, 
    p.mat=p.values2, sig.level=sig_level/length(M2))
colors <- ifelse(grepl('M|Y|Z|alpha', names(seis.DF)[1:9]), 
    '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2)
text(1:9, ncol(pcs)+.4, 
    as.expression(Map(function(x) bquote(.(x)[""]), seis.labs[1:9])),
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


