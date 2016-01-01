library(corrplot)
library(stargazer)
library(magicaxis)
library(RColorBrewer)

method <- 'spearman'
sig.level <- 10^-5
clevel <- 0.99

font <- 'Palatino'

plot_width = 6.97522
plot_height = 4.17309

DF <- read.table('grids/deleter.dat', header=1)
exclude <- which(grepl('Dnu_|mass|r_avg10|intercept', names(DF)))
X <- DF[-exclude]
M <- cor(X, method=method)

labs <- expression(M, Y[0], Z[0], alpha["MLT"], tau, "mass", R, 
    H[c], "X(He)", log~g, L, T["eff"], "Fe"/"H", 
    "<"*Delta*nu*">", "<"*d*Delta*nu/d*nu*">", #"<"*Delta*nu^b*">", 
    "<"*Delta*nu[0]*">", "<"*d*Delta*nu[0]/d*nu*">", #"<"*Delta*nu[0]^b*">", 
    "<"*delta*nu[0*","*2]*">", "<"*d*delta*nu[0*","*2]/d*nu*">", 
        #"<"*delta*nu[0*","*2]^b*">", 
    "<"*r[0*","*2]*">", "<"*d*r[0*","*2]/d*nu*">", #"<"*r[0*","*2]^b*">", 
    "<"*r[0*","*1]*">", "<"*d*r[0*","*1]/d*nu*">", #"<"*r[0*","*1]^b*">", 
    "<"*delta*nu[1*","*3]*">", "<"*d*delta*nu[1*","*3]/d*nu*">", 
        #"<"*delta*nu[1*","*3]^b*">", 
    "<"*r[1*","*3]*">", "<"*d*r[1*","*3]/d*nu*">",
    "<"*r[1*","*0]*">", "<"*d*r[1*","*0]/d*nu*">"#, "<"*r[0*","*1]^b*">", 
)

latex_labs <- c("M", "$Y_0$", "$Z_0$", "$\\alpha_{\\text{\"MLT\"}}$", 
    "$\\tau$", "mass", "R", "$H_c$", "X(He)", "$\\log g$", "L", 
    "$T_{\text{\"eff\"}}$", "Fe/H", 
    
    "$\\langle\\Delta\\nu\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu}{d\nu}\\rangle$", 
    
    "$\\langle\\Delta\\nu_0\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu_0}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{02}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{02}\\rangle$", 
    "$\\langle\\frac{dr_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{01}\\rangle$", 
    "$\\langle\\frac{dr_{01}}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{13}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{13}\\rangle$", 
    "$\\langle\\frac{dr_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{10}\\rangle$", 
    "$\\langle\\frac{dr_{10}}{d\nu}\\rangle$"
)

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

res1 <- cor.mtest(X, method)

cairo_pdf('plots/corr-spearman.pdf', width=8, height=7, family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0,0,0,0), cex=1, cex.lab=1)
a <- corrplot(M, diag=1, tl.col='white', type='lower', order="FPC", 
    tl.cex = 0.3, cl.cex=1, tl.srt=90, 
    p.mat = res1, sig.level = sig.level/length(M))
    #hclust.method="median"
pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(X))))
colors <- ifelse(grepl('M|Y|Z|alpha', names(X)[pos]), '#800080', 'black')
text(1:ncol(X)-0.3, (ncol(X)+1):2-0.4, col=colors, 
    labs[-exclude][pos], pos=4, srt=90)
text(0.8, ncol(X):1-0.1, labs[-exclude][pos], col=colors, pos=2)
dev.off()

#rownames(M) <- latex_labs[-exclude]
#colnames(M) <- latex_labs[-exclude]
#stargazer(M)



#exclude2 <- which(grepl('Dnu0|mass', names(DF)))
#X2 <- DF[-exclude2]
pca <- prcomp(X[,-1:-8], center=TRUE, scale.=TRUE)
pcs <- pca$x[,cumsum(pca$sdev)/sum(pca$sdev)<0.99]
vars <- formatC(pca$sdev / sum(pca$sdev) * 100, digits=4)
pclabs <- as.expression(lapply(1:ncol(pcs), 
    function(ii) bquote("("*.(vars[ii])*"%)"~PC[.(ii)])))


M2 <- cor(pcs, X[,1:8], method=method)
#colnames(M2) <- labs[-exclude2][1:8]
#rownames(M2) <- pclabs
#stargazer(M2)

cor.mtest2 <- function(mat1, mat2, test) {
    #mat1 <- as.matrix(mat1)
    #mat2 <- as.matrix(mat2)
    p.mat <- matrix(NA, ncol(mat1), ncol(mat2))
    for (i in 1:nrow(p.mat)) {
        for (j in 1:ncol(p.mat)) {
            p.mat[i, j] <- cor.test(mat1[, i], mat2[, j], method=test)$p.value
        }
    }
    return(p.mat)
}

res2 <- cor.mtest2(pcs, X[,1:8], method)

cairo_pdf('plots/corr-pca.pdf', width=8, height=7, family='Palatino')
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0,0,0,0), cex=1, cex.lab=1)
a <- corrplot(M2, tl.col='white', tl.cex = 0.3, cl.cex=1, tl.srt=90,
    cl.pos='b',
    p.mat = res2, sig.level = sig.level/length(M2))
colors <- ifelse(grepl('M|Y|Z|alpha', names(X)[1:8]), '#800080', 'black')
text(0.6, ncol(pcs):1, pclabs, pos=2)
text(1:8, ncol(pcs)+.4, labs[-exclude][1:8], col=colors, pos=3)
dev.off()

X3 <- cbind(X2[,1:8], pcs)
x3_labs <- c(labs[-exclude2][1:8], pclabs)
M3 <- cor(X3, method=method)
res2 <- cor.mtest(X, method)

col.pal <- brewer.pal(10, "Spectral")

png('plots/L-dnus-R.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- min(DF$radius)
varmax <- max(DF$radius)
cols <- round(10*((DF$radius-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(DF, 
    plot(log10(L), Dnu_slope, pch=3, cex=0.01, tcl=0,
         xlim=log10(range(L)),
         ylim=range(Dnu_slope),
         xlab=expression(L/L['\u0298']),
         ylab=expression("<"*d*Delta*nu/d*nu*">"),
         xaxs='i', xaxt='n', yaxt='n', #yaxs='i', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), 
    mgp=c(2, 0.25, 0), unlog='x')
dev.off()

png('plots/dnu-dnus-logg.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- min(DF$log_g)
varmax <- max(DF$log_g)
cols <- round(10*((DF$log_g-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(DF, 
    plot(Dnu_median, Dnu_slope, pch=3, cex=0.01, tcl=0,
         xlim=range(Dnu_median),
         ylim=range(Dnu_slope),
         xlab=expression(Delta*nu),
         ylab=expression("<"*d*Delta*nu/d*nu*">"),
         xaxs='i', xaxt='n', yaxt='n', #yaxs='i', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), 
    mgp=c(2, 0.25, 0))
dev.off()

png('plots/M-R-PC1.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC1, c(0.01))
varmax <- quantile(X3$PC1, c(0.99))
cols <- round(10*((X3$PC1-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(M, radius, pch=3, cex=0.01, tcl=0,
         xlim=range(M),
         ylim=range(radius),
         xlab=expression(M/M['\u0298']),
         ylab=expression(R/R['\u0298']),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/M-pc1-pc2.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- 0.7
varmax <- 1.3
cols <- round(10*((X3$M-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(PC1, PC5, pch=3, cex=0.01, tcl=0,
         xlim=quantile(PC1, c(0.01, 0.99)),
         ylim=quantile(PC5, c(0.01, 0.99)),
         xlab=expression(PC[1]),
         ylab=expression(PC[5]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/Hc-pc1-pc2.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- 0
varmax <- 0.78
cols <- round(10*((X3$Hc-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(PC1, PC2, pch=3, cex=0.01, tcl=0,
         xlim=quantile(PC1, c(0.01, 0.99)),
         ylim=quantile(PC2, c(0.01, 0.99)),
         xlab=expression(PC[1]),
         ylab=expression(PC[2]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/Hc-resid.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- 0.22#quantile(He, 0.01)
varmax <- 0.34#quantile(He, 0.99)
cols <- round(10*((X3$Y-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(predict(lm(Hc~pcs)), Hc, pch=3, cex=0.01, tcl=0,
         xlim=c(0, 0.8),
         ylim=c(0, 0.8),
         ylab=expression(H[c]),
         xlab=expression(Sigma[i]*c[i]*PC[i]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
#abline(lm(Hc~pcs))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/hc-age-ratio.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
#varmin <- 0.7
#varmax <- 1.3
#cols <- round(10*((X3$M-varmin)/(varmax-varmin)))
#cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(age ~ PC2, pch=3, cex=0.01, tcl=0,
         xlim=quantile(PC2, c(0.01, 0.99)),
         ylim=c(0, 14),#round(quantile(Hc/age, c(0.01, 0.99))),
         xlab=expression(PC[2]),
         ylab=expression(H[c]/tau),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n'))
         #col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/age-hc-PC2.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC2, c(0.01))
varmax <- quantile(X3$PC2, c(0.99))
cols <- round(10*((X3$PC2-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(age, Hc, pch=3, cex=0.01, tcl=0,
         xlim=range(age),
         ylim=range(Hc),
         xlab=expression(tau),
         ylab=expression(H[c]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=rev(col.pal)[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/M-Hc-PC1.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC1, c(0.01))
varmax <- quantile(X3$PC1, c(0.99))
cols <- round(10*((X3$PC1-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(M, Hc, pch=3, cex=0.01, tcl=0,
         xlim=round(range(M), 2),
         ylim=round(range(Hc), 1),
         xlab=expression(M/M['\u0298']),
         ylab=expression(H[c]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
color.legend(1.35, 0, 1.4, 0.77968, 
    signif(seq(varmin, varmax, length=10), 2), 
    col.pal[1:10], gradient='y', align='rb')
mtext(expression(PC[1]), 4, line=4.5)
dev.off()

png('plots/y0-alpha-PC8.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC8, c(0.01))
varmax <- quantile(X3$PC8, c(0.99))
cols <- round(10*((X3$PC8-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(Y, alpha, pch=20, cex=.1, tcl=0,
         xlim=range(Y),
         ylim=range(alpha),
         xlab=expression(Y[0]),
         ylab=expression(alpha["MLT"]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/hc-xhe-PC1.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC1, c(0.01))
varmax <- quantile(X3$PC1, c(0.99))
cols <- round(10*((X3$PC1-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(Hc, He, pch=3, cex=0.01, tcl=0,
         xlim=range(Hc),
         ylim=range(He),
         xlab=expression(H[c]),
         ylab=expression(X(He)),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/hc-r-PC1.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC1, c(0.01))
varmax <- quantile(X3$PC1, c(0.99))
cols <- round(10*((X3$PC1-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(Hc, radius, pch=3, cex=0.01, tcl=0,
         xlim=range(Hc),
         ylim=range(radius),
         xlab=expression(H[c]),
         ylab=expression(R/R['\u0298']),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()

png('plots/pc1-hc-pc2.png', width=250*plot_width, height=250*plot_height, 
    res=400, family=font)
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
varmin <- quantile(X3$PC2, c(0.01))
varmax <- quantile(X3$PC2, c(0.99))
cols <- round(10*((X3$PC2-varmin)/(varmax-varmin)))
cols <- ifelse(cols > 10, 10, ifelse(cols < 1, 1, cols))
with(X3, 
    plot(PC1, Hc, pch=3, cex=0.01, tcl=0,
         xlim=fivenum(PC1)[c(2,4)],
         ylim=range(Hc),
         xlab=expression(PC[1]),
         ylab=expression(H[c]),
         xaxs='i', yaxs='i', xaxt='n', yaxt='n', 
         col=col.pal[cols]))
magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=c(2, 0.25, 0))
dev.off()



#cairo_pdf('plots/corr-pca.pdf', width=8, height=7, family='Palatino')
#par(mar=c(0, 0, 0, 0), mgp=c(2, 0.25, 0), oma=c(0,0,0,0), cex=1, cex.lab=1)
#a <- corrplot(M3, diag=1, tl.col='white', type='lower', order="FPC", 
#    tl.cex = 0.3, cl.cex=1, tl.srt=90, 
#    p.mat = res2, sig.level = sig.level)
#    #hclust.method="median"
#pos <- as.numeric(sapply(colnames(a), function(x) which(x==names(X3))))
#colors <- ifelse(grepl('M|Y|Z|alpha', names(X3)[pos]), '#800080', 'black')
#text(1:ncol(X3)-0.3, (ncol(X3)+1):2-0.4, col=colors, x3_labs[pos], pos=4, srt=90)
#text(0.8, ncol(X3):1, x3_labs[pos], col=colors, pos=2)
#dev.off()

#pca <- prcomp(X[-1:-10], center=TRUE, scale.=TRUE)
