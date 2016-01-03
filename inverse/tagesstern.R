#### Perturb solar frequencies to the level of Kepler uncertainty 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source('../scripts/seismology.R')
library(randomForest)

plot_dir <- file.path('plots', 'Tagesstern')

sun <- read.table(file.path('data', 'Sun-freqs.dat'), header=1)
cygA <- read.table(file.path('data', '16CygA-freqs.dat'), header=1)
cygB <- read.table(file.path('data', '16CygB-freqs.dat'), header=1)

sun$nu <- sun$nu - 3090
cygA$nu <- cygA$nu - 2201
cygB$nu <- cygB$nu - 2552

relation <- randomForest(dnu~n+l+nu, data=cygA)
y <- predict(relation)

freq_range <- range(cygA$nu, cygB$nu)*1.025
tag <- sun[sun$nu >= freq_range[1] & sun$nu <= freq_range[2],]

y2 <- predict(relation, newdata=tag)
new_nu <- rnorm(nrow(tag), tag$nu, y2)

#cairo_pdf('plots/tagesstern-sigma.pdf', width=6, height=4.5, family='Palatino')
#par(mar=c(3, 4, 1, 1), mgp=c(2, 0.25, 0), cex.lab=1.3)

plot_tag <- function(text.cex, ...) {
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
    
    points(cygA$nu, y, pch=1, col='darkred', cex=0.75)
    segments(cygA$nu, y, cygA$nu, cygA$dnu, 
             col=adjustcolor('darkred', alpha.f=0.25))

    
    abline(v=freq_range, lty=3)
    #abline(v=min(sun$nu), lty=3)
    #abline(v=max(sun$nu), lty=3)

    
    points(tag$nu, y2, pch=1, cex=0.75)
    segments(new_nu, y2, tag$nu, tag$dnu, col=adjustcolor('black',alpha.f=0.25))
}
#dev.off()

make_plots(plot_tag, 'Tagesstern', filepath=plot_dir)

tag$nu <- new_nu + 3108
tag$dnu <- y2

seismology(tag, nu_max=3108, outf='Tagesstern', filepath=plot_dir)
write.table(tag, file.path('data', 'Tagesstern-freqs.dat'), 
    quote=FALSE, row.names=FALSE)
