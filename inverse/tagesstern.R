#### Perturb solar frequencies to the level of Kepler uncertainty 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'seismology.R'))
library(randomForest)

plot_dir <- file.path('plots', 'Tagesstern')

sun.freqs <- read.table(file.path('data', 'Sun-freqs.dat'), header=1)
cygA.freqs <- read.table(file.path('data', '16CygA-freqs.dat'), header=1)
cygB.freqs <- read.table(file.path('data', '16CygB-freqs.dat'), header=1)

sun.obs <- read.table(file.path('data', 'Sun-obs.dat'), header=1)
cygA.obs <- read.table(file.path('data', '16CygA-obs.dat'), header=1)
cygB.obs <- read.table(file.path('data', '16CygB-obs.dat'), header=1)

sun.freqs$nu <- sun.freqs$nu - sun.obs[sun.obs$name == 'nu_max',]$value
cygA.freqs$nu <- cygA.freqs$nu - cygA.obs[cygA.obs$name == 'nu_max',]$value
cygB.freqs$nu <- cygB.freqs$nu - cygB.obs[cygB.obs$name == 'nu_max',]$value

ab.freqs <- rbind(cygA.freqs, cygB.freqs)
relation <- randomForest(dnu~l+nu+n, data=ab.freqs)
y <- predict(relation)
plot(ab.freqs$dnu, 100*(ab.freqs$dnu-y)/ab.freqs$dnu)

tag.freqs <- sun.freqs[-1:-(nrow(sun.freqs)-55),]

y2 <- predict(relation, newdata=tag.freqs)
new_nu <- rnorm(nrow(tag.freqs), tag.freqs$nu, y2)

plot_tag <- function(..., text.cex=1) {
    plot(sun.freqs$nu, sun.freqs$dnu, pch=20, 
         ylim=range(sun.freqs$dnu, cygA.freqs$dnu, y2),
         #xlab=expression((nu-nu[max])/mu*Hz),
         #ylab=expression("uncertainty"~sigma[nu]/mu*Hz),
         tcl=0)
    magaxis(1:4, tcl=0.25, labels=0)
    points(cygA.freqs$nu, cygA.freqs$dnu, pch=20, col='darkred')
    points(cygB.freqs$nu, cygB.freqs$dnu, pch=20, col='blue', cex=0.75)
    points(tag.freqs$nu, y2, pch=1, cex=0.75)
    segments(new_nu, y2, tag.freqs$nu, tag.freqs$dnu, 
        col=adjustcolor('black',alpha.f=0.25))
    abline(v=0, lty=2)
    legend("topleft", 
           legend=c(expression("("*nu["16CygA"]*", "*sigma["16CygA"]*")"),
                    expression("("*nu["16CygB"]*", "*sigma["16CygB"]*")"),
                    #expression("("*nu["16CygA"]*", "*hat(sigma)["16CygA"]*")"),
                    expression("("*nu["Sun"]*", "*sigma["Sun"]*")"),
                    expression("("*hat(nu)["Sun"]*", "*hat(sigma)["Sun"]*")"),
                    expression(nu[max])),
           pch=c(20, 20, 20, 1, NA), 
           col=c("darkred", "blue", "black", "black", "black"),
           lty=c(NA, NA, NA, NA, 2),
           pt.cex=c(1, 0.75, 1, 0.75),
           bty='n')
}

make_plots(plot_tag, 'Tagesstern', filepath=plot_dir)

tag.freqs$nu <- new_nu + sun.obs[sun.obs$name == 'nu_max',]$value
tag.freqs$dnu <- y2

## Generate uncertainties for classical observations
tag.obs <- sun.obs
rel.unc <- unlist(Map(mean, abs(cygA.obs$uncertainty / cygB.obs$value), 
                            abs(cygB.obs$uncertainty / cygB.obs$value) ))
tag.obs$uncertainty <- rel.unc * tag.obs$value
zeros <- tag.obs$uncertainty == 0
tag.obs$uncertainty[zeros] <- log10(10**rel.unc[zeros] * 
    10**tag.obs$value[zeros])

# Perturb the value within uncertainty one time 
tag.obs$value <- rnorm(nrow(tag.obs), tag.obs$value, tag.obs$uncertainty)

write.table(tag.obs, file.path('data', 'Tagesstern-obs.dat'), 
    quote=FALSE, row.names=FALSE)
seismology(tag.freqs, nu_max=tag.obs[tag.obs$name=="nu_max",]$value, 
    outf='Tagesstern', filepath=plot_dir)
write.table(tag.freqs, file.path('data', 'Tagesstern-freqs.dat'), 
    quote=FALSE, row.names=FALSE)

