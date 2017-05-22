source('../../scripts/utils.R')

nu.cols <- c('l', 'n', 'nu')

dat <- read.table(file.path('diffusion3', 'LOGS_MS', 'profile1-freqs', 'profile1-freqs.dat'),
    col.names=nu.cols)
var <- read.table(file.path('diffusion3', 'LOGS_MS', 'profile1-freqs', 'profile1-freqs.var'),
    col.names=nu.cols)
comb <- merge(dat, var, by=nu.cols[1:2])

plot_var_rich <- function(..., text.cex=1) {
    with(comb, plot(nu.x, nu.y-nu.x, axes=F, tcl=0, pch=l+1,
	    xlab=expression("Richardson's extrapolated frequencies"~nu['R']/mu*Hz),
		ylab=expression("Difference with variational frequencies"~(nu['v']-nu['R'])/mu*Hz)))
    magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25, cex.axis=text.cex)
	legend('topleft', pch=1:4, legend=c("l=0", "l=1", "l=2", "l=3"))
}

make_plots(plot_var_rich, "var-rich")

dat <- read.table(file.path('diffusion', 'LOGS_MS', 'profile1-freqs', 'profile1-freqs.dat'),
    col.names=nu.cols)
var <- read.table(file.path('diffusion', 'LOGS_MS', 'profile1-freqs', 'profile1-freqs.var'),
    col.names=nu.cols)
comb <- merge(dat, var, by=nu.cols[1:2])

make_plots(plot_var_rich, "var-rich2")
