#### Mesh and scatterplot analysis of evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source(file.path('..', 'scripts', 'utils.R'))

library(ggplot2)
library(GGally)
library(scales)

## Load data
combos <- read.table(file.path('..', 'forward', 'initial_conditions.dat'), 
    col.names=c("M", "Y", "Z", "alpha", "overshoot", "diffusion"))#[1:1000,]

log_vars <- c(3, 5, 6)

# Make inputs diagram
X <- 1-combos$Y-combos$Z

# set up corner
p = ggpairs(data=combos, axisLabels="show", upper="blank",
        aes(color=X, size=0.01), 
        lower = list(continuous = wrap("points", alpha = 0.75, size=0.1)),
        columnLabels=sapply(names(combos)[1:6], 
            function (name) as.expression(get_label_nameless(name))))

# log plot alpha_mlt, alpha_ov, and D
for (col_j in 1:6) {
    for (row_i in 1:6) {
        if (row_i == col_j) {
            p <- putPlot(p, 
                ggplot(combos) + geom_line(aes_string(names(combos)[row_i]), 
                   stat="density"), row_i, row_i)
            next
        }
        pp <- getPlot(p, row_i, col_j)
        if (row_i %in% log_vars)
            pp <- pp + scale_y_log10(limits=signif(range(combos[,row_i]), 1))
        if (col_j %in% log_vars)
            pp <- pp + scale_x_log10(limits=signif(range(combos[,col_j]), 1))
        if (row_i %in% log_vars || col_j %in% log_vars) {
            pp$subtype <- 'logpoints'
            pp$type <- 'logcontinuous'
        }
        p <- putPlot(p, pp, row_i, col_j)
    }
}

# change number lines
number_ticks <- function(xs) {
    repr <- seq(min(xs), max(xs), length.out=1000)
    signif(as.numeric(quantile(repr, c(0.25, 0.75))), 2)
} 
log_number_ticks <- function(xs) {
    repr <- seq(log10(min(xs)), log10(max(xs)), length.out=1000)
    signif(10**as.numeric(quantile(repr, c(0.25, 0.75))), 2)
} 
minor_ticks <- function(xs) {
    repr <- seq(min(xs), max(xs), length.out=1000)
    signif(as.numeric(quantile(repr, 0:10/10), 2))
} 
for (ii in 1:6) {
    for (jj in 1:6) {
        pp <- getPlot(p, ii, jj)
        if (jj %in% log_vars) {
            pp <- pp + scale_x_log10(breaks=log_number_ticks)
        } else {
            pp <- pp + scale_x_continuous(breaks=number_ticks)
        }
        if (ii %in% log_vars && ii != jj) {
            pp <- pp + scale_y_log10(breaks=log_number_ticks)
        } else {
            pp <- pp + scale_y_continuous(breaks=number_ticks)
        }
        #pp <- pp + scale_size_area(range=c(1, 1))
        p <- putPlot(p, pp, ii, jj)
    }
}

# remove the labels on the top left plot 
pp <- getPlot(p, 1, 1)
pp <- pp + theme(axis.line=element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank(),
                 text=element_blank())
p <- putPlot(p, pp, 1, 1)

# turn into a function so that make_plots can call it
inputs_plot <- function(..., text.cex) {
    print(p,
    #leftWidthProportion=0.6, bottomHeightProportion=0.5
    leftWidthProportion=0.4, bottomHeightProportion=0.3
    )
}

# save!
make_plots(inputs_plot, "inputs", filepath=file.path("plots", "inputs"),
    short=FALSE, thin=FALSE, make_png=TRUE,
        paper_pdf_width=1.3*6.97522, 
        paper_pdf_height=1.3*4.17309)

