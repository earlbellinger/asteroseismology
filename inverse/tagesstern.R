library(magicaxis)
library(randomForest)

source('seismology.R')

sun <- read.table('../data/Sun-freqs.dat', header=1)
cygA <- read.table('../data/16CygA-freqs.dat', header=1)
cygB <- read.table('../data/16CygB-freqs.dat', header=1)

sun$nu <- sun$nu - 3108
cygA$nu <- cygA$nu - 2201
cygB$nu <- cygB$nu - 2552

cairo_pdf('plots/tagesstern-sigma.pdf', width=6, height=4.5, family='Palatino')
par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)
plot(sun$nu, sun$dnu, pch=20, 
     ylim=range(sun$dnu, cygA$dnu),
     xlab=expression(nu-nu[max]~"["*mu*Hz*"]"),
     ylab=expression("uncertainty"~sigma[nu]~"["*mu*Hz*"]"),
     tcl=0)
points(cygA$nu, cygA$dnu, pch=20, col='darkred')

abline(v=0, lty=2)
magaxis(1:4, tcl=0.25, labels=0)
legend("topleft", 
       legend=c(expression("("*nu["16CygA"]*", "*sigma["16CygA"]*")"),
                expression("("*nu["16CygA"]*", "*hat(sigma)["16CygA"]*")"),
                expression("("*nu["Sun"]*", "*sigma["Sun"]*")"),
                expression("("*hat(nu)["Sun"]*", "*hat(sigma)["Sun"]*")"),
                expression(nu[max]),
                expression(nu["cutoff"])),
       pch=c(20, 1, 20, 1, NA, NA), 
       col=c("darkred", "darkred", "black", "black", "black", "black"),
       lty=c(NA, NA, NA, NA, 2, 3),
       pt.cex=c(1, 0.75, 1, 0.75),
       bty='n')

relation <- randomForest(dnu~n+l+nu, data=cygA)
y <- predict(relation)
points(cygA$nu, y, pch=1, col='darkred', cex=0.75)
segments(cygA$nu, y, cygA$nu, cygA$dnu, 
         col=adjustcolor('darkred', alpha.f=0.25))

freq_range <- range(cygA$nu, cygB$nu)*1.025
sun <- sun[sun$nu >= freq_range[1] & sun$nu <= freq_range[2],]
abline(v=freq_range, lty=3)
#abline(v=min(sun$nu), lty=3)
#abline(v=max(sun$nu), lty=3)

y <- predict(relation, newdata=sun)
new_nu <- rnorm(nrow(sun), sun$nu, y)
points(sun$nu, y, pch=1, cex=0.75)
segments(new_nu, y, sun$nu, sun$dnu, col=adjustcolor('black', alpha.f=0.25))
dev.off()

sun$nu <- new_nu + 3108
sun$dnu <- y

seismology(sun, nu_max=3108, outf='Tagesstern')

write.table(sun, '../data/Tagesstern-freqs.dat', quote=FALSE, row.names=FALSE)

