source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)

n <- 50
sigma <- 3
x <- sort(runif(n, 0, 10))
y <- 2*x+1 + rnorm(n, 0, sigma)

interp_vs_reg <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot(x,y, axes=0, tcl=0)
    lines(approx(x,y), col=red)
    abline(lm(y~x), lty=2, col=blue)
    legend("topleft", pch=c(1, NA, NA), lty=c(NA, 1, 2), 
        col=c('black', red, blue), bty='n',
        legend=c("Data", "Linear interpolation", "Linear regression"))
    magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0), mgp=mgp,
        las=1, cex.axis=text.cex)
}

make_plots(interp_vs_reg, "interpolation_vs_regression")

