#### Utility functions for forward & inverse asteroseismology 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Plotting values
utils.mgp <<- c(2, 0.25, 0)

font <- "Palatino"
png_res <- 400
slide_cex.lab <- 1.3

paper_pdf_width <- 6.97522 # inches
paper_pdf_height <- 4.17309

latex_pt_per_in <- 5 * 72.27
paper_png_width <- paper_pdf_width * latex_pt_per_in
paper_png_height <- paper_pdf_height * latex_pt_per_in

slide_pdf_width <- 6.22665 
slide_pdf_height <- 4.1511

slide_png_width <- slide_pdf_width * latex_pt_per_in
slide_png_height <- slide_pdf_height * latex_pt_per_in

## Make the same plot as a pdf and a png suitable for papers and for slides
# takes a plotting function plot_f that calls `plot` 
make_plots <- function(plot_f, filename, ..., 
        filepath='plots', mar=c(3, 4, 1, 1), mgp=utils.mgp, 
        wide=TRUE, thin=TRUE, 
        make_png=TRUE, make_pdf=TRUE, 
        paper=TRUE, slide=TRUE) {
    
    if (paper & make_pdf & wide) {
        directory <- file.path(filepath, 'paper', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '.pdf')),
                  width=paper_pdf_width, height=paper_pdf_height, 
                  family=font)
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=1, family=font, ...)
        plot_f(text.cex=1, ...)
        dev.off()
    }
    
    if (paper & make_pdf & thin) {
        directory <- file.path(filepath, 'paper', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '-thin.pdf')),
                  width=paper_pdf_width/2, height=paper_pdf_height, 
                  family=font)
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=1, family=font)
        plot_f(text.cex=1, ...)
        dev.off()
    }
    
    if (paper & make_png & wide) {
        directory <- file.path(filepath, 'paper', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '.png')),
                  width=paper_png_width, height=paper_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=1, family=font)
        plot_f(text.cex=1, ...)
        dev.off()
    }
    
    if (paper & make_png & thin) {
        directory <- file.path(filepath, 'paper', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '-thin.png')),
                  width=paper_png_width/2, height=paper_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=1, family=font)
        plot_f(text.cex=1, ...)
        dev.off()
    }
    
    if (slide & make_pdf & wide) {
        directory <- file.path(filepath, 'slides', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '-slides.pdf')),
                  width=slide_pdf_width, height=slide_pdf_height, 
                  family=font)
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=slide_cex.lab, family=font)
        plot_f(text.cex=slide_cex.lab, ...)
        dev.off()
    }
    
    if (slide & make_pdf & thin) {
        directory <- file.path(filepath, 'slides', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        cairo_pdf(file.path(directory, paste0(filename, '-slides-thin.pdf')),
                  width=slide_pdf_width/2, height=slide_pdf_height, 
                  family=font)
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=slide_cex.lab, family=font)
        plot_f(text.cex=slide_cex.lab, ...)
        dev.off()
    }
    
    if (slide & make_png & wide) {
        directory <- file.path(filepath, 'slides', 'wide')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '-slides.png')),
                  width=slide_png_width, height=slide_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=slide_cex.lab, family=font)
        plot_f(text.cex=slide_cex.lab, ...)
        dev.off()
    }
    
    if (slide & make_png & thin) {
        directory <- file.path(filepath, 'slides', 'thin')
        dir.create(directory, showWarnings=FALSE, recursive=TRUE)
        png(file.path(directory, paste0(filename, '-slides-thin.png')), 
                  width=slide_png_width/2, height=slide_png_height, 
                  family=font, res=png_res, type='cairo')
        par(mar=mar, mgp=c(2, 0.25, 0), cex.lab=slide_cex.lab, family=font)
        plot_f(text.cex=slide_cex.lab, ...)
        dev.off()
    }
}
