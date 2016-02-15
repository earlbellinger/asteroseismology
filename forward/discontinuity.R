#### Mesh and scatterplot analysis of evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

suppressMessages(source(file.path('..', '..', '..', 'scripts', 'utils.R')))

suppressMessages(library(magicaxis))

DF <- read.table('LOGS/history.data', skip=5, header=1)
decreasing_L <- which(diff(DF$log_L) < 0 & DF$star_age[-1] < 0.25*10**9)
if (any(decreasing_L)) {
    goes_back_up <- diff(decreasing_L) > 1
    pms <- ifelse(any(goes_back_up), 
                  which(goes_back_up)[1] + 1, 
                  max(decreasing_L))
    DF <- DF[-1:-pms,]
}
DF <- DF[floor(nrow(DF)*0.05):(nrow(DF)-floor(nrow(DF)*0.05)),] 
#DF <- DF[-1:-floor(nrow(DF)*0.01),] 

y <- DF$log_Teff
x <- DF$star_age

new_xs <- seq(min(x), max(x), length=100000)
spl <- splinefun(x, y)(new_xs)
roc <- splinefun(x, y)(new_xs, 1)

outliers <- boxplot.stats(roc, coef=15)$out
locs <- which(roc %in% outliers)

plot_ts <- function(text.cex=1, font="Palatino", mgp=c(2, 0.25, 0), ...) {
    plot(x, y, cex=0.01, pch=1,
         tcl=0, axes=FALSE, 
         ylab=expression('Temperature' ~ 'lg'*(T[eff]/K)),
         xlab=expression('Star age' ~ tau/'yr'))
    magaxis(side=1:4, family=font, tcl=0.25, mgp=mgp, las=1, 
            cex.axis=text.cex, labels=c(1,1,0,0))
    points(new_xs[locs], spl[locs], col=adjustcolor("red", alpha=0.25), cex=5)
}

inlist <- readLines("inlist_1.0")
line_num <- grep('mesh_delta_coeff = ', inlist)
mesh_delta <- sub(".+= ", "", inlist[line_num])

invisible(make_plots(plot_ts, 
            paste0(basename(getwd()), "-discontinuity-", mesh_delta),
            filepath=file.path('..', '..', 'plots', 'discontinuity')))

cat(as.numeric(length(outliers) > 0))
