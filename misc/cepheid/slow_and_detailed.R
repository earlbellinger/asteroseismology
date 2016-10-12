library(magicaxis)
library(parallel)
library(parallelMap)
library(RColorBrewer)

log.dir <- file.path('work', 'LOGS')

hist.DF <- read.table(file.path(log.dir, 'history.data'),
    header=1, skip=5)

log.files <- list.files(log.dir)
profile.files <- log.files[grepl('profile\\d+\\.data$', log.files)]

pro.file <- profile.files[1]
early <- as.numeric(sub('.data', '', sub('profile', '', profile.files)))>0
h1.vals <- c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
m.vals <- seq(0.01, 1, 0.01)

parallelStartMulticore(16)

### Kippenhahn of eps_nuc/age/mass
eps <- do.call(rbind, parallelMap(function(pro.file) {
        pro.head <- read.table(file.path(log.dir, pro.file), header=1, nrow=1,
            skip=1)
        age <- pro.head$star_age/10**6
        pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
        eps_nuc <- pro.DF$eps_nuc
        cutoff <- eps_nuc > 0.1
        eps_nuc <- log10(eps_nuc[cutoff])
        mass <- pro.DF$mass[cutoff]
        #eps_nuc <- log10(pro.DF$eps_nuc)
        #mass <- pro.DF$mass
        contours <- approx(mass, eps_nuc, m.vals, rule=2)$y
        data.frame(age, rbind(contours))
    }, profile.files[early]))
eps <- eps[order(eps$age),]
conv <- do.call(rbind, parallelMap(function(pro.file) {
        pro.head <- read.table(file.path(log.dir, pro.file), header=1, nrow=1,
            skip=1)
        age <- pro.head$star_age/10**6
        pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
        grad.diff <- pro.DF$gradr - pro.DF$grada
        mass <- pro.DF$mass
        contours <- approx(mass, grad.diff, m.vals, rule=2)$y
        data.frame(age, rbind(contours))
    }, profile.files[early]))
conv <- conv[order(conv$age),]
cairo_pdf('slow_and_detailed-eps_nuc.pdf', family='Palatino')
par(mar=c(4, 5, 1, 7))
filled.contour(eps$age, log10(m.vals), as.matrix(eps[,-1]),
    color=colorRampPalette(brewer.pal(11, "Spectral")),
    xaxs='i', yaxs='i',
    key.axes={
        axis(4, tcl=0, line=0)
        mtext(expression("Nuclear Energy"~log[10](epsilon/erg/g/s)), 
            side=4, las=3, line=3)
    },
    plot.axes={
        contour(conv$age, log10(m.vals), as.matrix(conv[,-1]), levels=c(0),
            drawlabels=F, col='white', add=TRUE)
        magaxis(side=1:4, family='Palatino', tcl=0.25, labels=c(1,1,0,0),
           unlog='y')
    },
    plot.title={
        title(xlab=expression("Star age"~tau/"Myr"), line=2)
        title(ylab=expression(m/M["*"]), line=2)
    })
dev.off()


m.vals <- seq(0, 1, 0.001)
### Kippenhahn of X/age/mass
h1s <- do.call(rbind, parallelMap(function(pro.file) {
        pro.head <- read.table(file.path(log.dir, pro.file), header=1, nrow=1,
            skip=1)
        age <- pro.head$star_age/10**6
        pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
        mass <- pro.DF$mass
        h1 <- pro.DF$h1
        contours <- approx(mass, h1, m.vals, rule=2)$y#h1.vals)$y
        #contours[is.na(contours)] <- 0
        data.frame(age, rbind(contours))
    }, profile.files[early]))
h1s <- h1s[order(h1s$age),]
cairo_pdf('slow_and_detailed-H1.pdf', family='Palatino')
par(mar=c(4, 5, 1, 7))
filled.contour(h1s$age, m.vals, as.matrix(h1s[,-1]),
    color=colorRampPalette(brewer.pal(11, "Spectral")),
    xaxs='i', yaxs='i',
    key.axes={
        axis(4, tcl=0, line=0)
        mtext("Hydrogen Abundance X", side=4, las=3, line=3)
    },
    plot.axes={
        contour(h1s$age, m.vals, as.matrix(h1s[,-1]), 
            add=TRUE, labcex=0.5, levels=h1.vals)
        magaxis(side=1:4, family='Palatino', tcl=0.25, labels=c(1,1,0,0))
    },
    plot.title={
        title(xlab=expression("Star age"~tau/"Myr"), line=2)
        title(ylab=expression(m/M["*"]), line=2)
    })
dev.off()

### Kippenhahn of X/age/radius
m.vals <- seq(0.01, 1, 0.01)
h1s <- do.call(rbind, parallelMap(function(pro.file) {
        pro.head <- read.table(file.path(log.dir, pro.file), header=1, nrow=1,
            skip=1)
        age <- pro.head$star_age/10**6
        pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
        radius <- pro.DF$radius
        h1 <- pro.DF$h1
        contours <- approx(radius, h1, m.vals, rule=2)$y#h1.vals)$y
        #contours[is.na(contours)] <- 0
        data.frame(age, rbind(contours))
    }, profile.files[early]))
h1s <- h1s[order(h1s$age),]
cairo_pdf('slow_and_detailed-H1_r.pdf', family='Palatino')
par(mar=c(4, 5, 1, 2))
contour(h1s$age, log10(m.vals), as.matrix(h1s[,-1]), levels=h1.vals, axes=F)
magaxis(side=1:4, family='Palatino', tcl=0.25, labels=c(1,1,0,0),
       unlog='y')
title(xlab=expression("Star age"~tau/"Myr"), line=2)
title(ylab=expression(r/R["*"]), line=2)
dev.off()

### Kippenhahn of X/model/radius
h1s <- do.call(rbind, parallelMap(function(pro.file) {
        pro.head <- read.table(file.path(log.dir, pro.file), header=1, nrow=1,
            skip=1)
        model <- pro.head$model_number
        pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
        radius <- pro.DF$radius
        h1 <- pro.DF$h1
        contours <- approx(radius, h1, m.vals, rule=2)$y#h1.vals)$y
        #contours[is.na(contours)] <- 0
        data.frame(model, rbind(contours))
    }, profile.files[early]))
h1s <- h1s[order(h1s$model),]
cairo_pdf('slow_and_detailed-H1_r_model.pdf', family='Palatino')
par(mar=c(4, 5, 1, 7))
filled.contour(h1s$model, m.vals, as.matrix(h1s[,-1]),
    color=colorRampPalette(brewer.pal(11, "Spectral")),
    xaxs='i', yaxs='i',
    key.axes={
        axis(4, tcl=0, line=0)
        mtext("Hydrogen Abundance X", side=4, las=3, line=3)
    },
    plot.axes={
        contour(h1s$model, m.vals, as.matrix(h1s[,-1]), 
            add=TRUE, labcex=0.5, levels=h1.vals)
        magaxis(side=1:4, family='Palatino', tcl=0.25, labels=c(1,1,0,0))
    },
    plot.title={
        title(xlab=expression("Model number"), line=2)
        title(ylab=expression(r/R["*"]), line=2)
    })
dev.off()




plot(NA, axes=F, xaxs='i', yaxs='i', #log='xy',
    xlim=c(1, max(h1s$age)),#range(h1s$age), 
    ylim=range(h1s[,-1], na.rm=1),
    xlab=expression("Star age"~tau/"Myr"), 
    ylab=expression(m/M["*"]))
magaxis(1:4, labels=c(1,1,0,0), tcl=0.25, family='Palatino')
cnames <- colnames(h1s[,-1])
for (colname_i in 1:length(cnames)) {
    #lty <- as.numeric(substring(colname, 2))+1
    lines(h1s$age, h1s[[cnames[colname_i]]], lty=colname_i)
}
legend('topleft', lty=1:(ncol(h1s)-1), bty='n',
    legend=as.expression(Map(function(x) bquote(X==.(x)), h1.vals)))


pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)

