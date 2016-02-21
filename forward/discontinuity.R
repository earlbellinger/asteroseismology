#### Discontinuity detection in MESA evolutionary tracks 
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
DF <- DF[-1:-floor(nrow(DF)*0.05),] # in case more PMS leaks in

y <- DF$surface_h1
x <- DF$star_age

outliers <- abs(diff(y)) >= 0.01
locs <- which(outliers)
locs <- locs[c(1, which(diff(locs)>1)+1)]

plot_ts <- function(text.cex=1, font="Palatino", mgp=utils.mgp, ...) {
    plot(x, y, cex=0.01, pch=1,
         tcl=0, axes=FALSE, 
         ylab=expression('Surface hydrogen' ~ "X"["surf"]),
         xlab=expression('Star age' ~ tau/'yr'))
    magaxis(side=1:4, family=font, tcl=0.25, mgp=mgp, las=1, 
            cex.axis=text.cex, labels=c(1,1,0,0))
    points(x[locs], y[locs], col=adjustcolor("red", alpha=0.75), cex=5)
}

inlist <- readLines("inlist_1.0")
line_num <- grep('mesh_delta_coeff = ', inlist)
mesh_delta <- sub(".+= ", "", inlist[line_num])
line_num <- grep('max_years_for_timestep = ', inlist)
time_step <- sub(".+= ", "", inlist[line_num])

invisible(make_plots(plot_ts, 
            paste0(basename(getwd()), "-discontinuity-", 
                time_step, "_", mesh_delta),
            filepath=file.path('..', '..', 'plots', 'discontinuity')))

cat(as.numeric(length(which(outliers)) > 0))

