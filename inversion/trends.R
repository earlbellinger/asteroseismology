#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

source('../scripts/utils.R') 
library(RColorBrewer)

DF <- read.table('inv-res.dat', header=1)
DF <- DF[order(DF$M),]
#col.pal <- colorRampPalette(brewer.pal(11, 'Spectral'))(nrow(DF))
#col.pal <- colorRampPalette(brewer.pal(12, 'Set1'))(nrow(DF))
#col.pal <- colorRampPalette(c(red, 'darkgray', blue))(nrow(DF))
col.pal <- colorRampPalette(c('#34738f', '#122f3d', '#be3e2b', 
    '#ed8a45', '#eec70e'))(nrow(DF))

plot_params <- function(varname, tr, make_xlab=T, make_ylab=T, ..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3,-0.4,-0.4,-0.4), mgp=mgp+c(0, 0.4, 0), lwd=1.5)
    
    plot(NA, axes=F, xaxs='i', yaxs='i',
        xlim=if (varname=="age")   c(0, 16) else
             if (varname=="M")     c(0.7, 1.5) else
             if (varname=="R")     c(0.9, 1.8) else
             if (varname=="FeH")   c(-0.6, 0.6) else #c(-1.1,  0.5) else
             if (varname=="alpha") c(1.2, 2.2),
        ylim=c(-0.25, 0.25),
        xlab="", ylab="")
    
    abline(h=0, lty=2, lwd=par()$lwd)
    
    #col <- if (varname=="age") 1 else
    #       if (varname=="M")   red else
    #       if (varname=="R")   blue else
    #       if (varname=="FeH") orange else
    #       if (varname=="alpha") "#723d46"
    col <- col.pal
    
    #with(DF, points(get(varname), get(tr), pch=20, col=1))
    with(DF,
        arrows(get(varname) - get(paste0('e_', varname)), 
           get(tr), #get(paste0('e_', tr)),
           get(varname) + get(paste0('e_', varname)), 
           get(tr), #get(paste0('e_', tr)),
           length=0, #.05, 
           code=3, lwd=par()$lwd*2, angle=90, col=col))
           #adjustcolor(col, alpha.f=0.9)))
    with(DF,
        arrows(get(varname), 
           get(tr) + get(paste0('e_', tr)),
           get(varname), 
           get(tr) - get(paste0('e_', tr)),
           length=0, #.05, 
           code=3, lwd=par()$lwd*2, angle=90, col=col))
           #adjustcolor(col, alpha.f=0.9)))
    
    #with(DF, points(get(varname), get(tr), pch=20, cex=1.5))
    #with(DF, points(get(varname), get(tr), pch=20, col=1, cex=0.66))
    with(DF, points(get(varname), get(tr), pch=21, col=1, lwd=1.5, 
        cex=0.66, bg=col.pal))
    
    
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, labels=make_xlab, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=make_ylab, lwd.ticks=par()$lwd)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.55, 0, 0))
    if (make_xlab) title(xlab=if (varname=="age") "Age [Gyr]" else
               if (varname=="M")   "Mass" else
               if (varname=="R")   "Radius" else
               if (varname=="FeH") "[Fe/H]" else
               if (varname=="alpha") expression(alpha["MLT"]))
    par(mgp=mgp+c(0.8, 0, 0))
    target_radius <- if (tr=='tr005') '0.05' else
                     if (tr=='tr01')  '0.1' else
                     if (tr=='tr015') '0.15' else
                     if (tr=='tr02')  '0.2' else
                     if (tr=='tr025') '0.25' else
                     if (tr=='tr03')  '0.3'
    if (make_ylab) title(ylab=bquote(delta*u/u*"("*.(target_radius)*")"))
}

for (varname in c('age', 'M', 'R', 'FeH', 'alpha')) {
    for (tr in c('tr005', 'tr01', 'tr015', 'tr02', 'tr025', 'tr03')) {
        make_plots(plot_params, paste0(varname, tr), 
            filepath=file.path('plots', 'trends'),
            wide=F, tall=F, slides=F, make_png=F,
            varname=varname,
            tr=tr,
            cex.paper=1.16,
            make_xlab=tr=='tr03',
            make_ylab=varname=='age',
            paper_pdf_height=4.5, #6.75206*2/3*5/6,#*12/11,
            paper_pdf_width=5) #4.17309*6/4*4/5)
    }
}


if (F) {
#varname <- 'alpha'
with(DF, plot(get(varname),   tr005, pch=20, col=1))
with(DF, points(get(varname), tr01,  pch=20, col=2))
with(DF, points(get(varname), tr015, pch=20, col=3))
with(DF, points(get(varname), tr02,  pch=20, col=4))
with(DF, points(get(varname), tr025, pch=20, col=5))
with(DF, points(get(varname), tr03,  pch=20, col=6))
abline(h=0, lty=2)


varname <- 'e_M'
op <- par(ask=TRUE)
for (tr in c('tr005', 'tr01', 'tr015', 'tr02', 'tr025', 'tr03')) {
    with(DF, 
        plot(get(varname), get(tr), pch=20, col=1, ylim=c(-0.15, 0.15),
            xlab=varname, ylab=tr))
        
    with(DF,
        arrows(get(varname) - get(paste0('e_', varname)), 
               get(tr), #get(paste0('e_', tr)),
               get(varname) + get(paste0('e_', varname)), 
               get(tr), #get(paste0('e_', tr)),
               code=3, length=0.05, angle=90, col=adjustcolor(1, alpha.f=0.5)))
    with(DF,
        arrows(get(varname), 
               get(tr) + get(paste0('e_', tr)),
               get(varname), 
               get(tr) - get(paste0('e_', tr)),
               code=3, length=0.05, angle=90, col=adjustcolor(1, alpha.f=0.5)))
        
    with(DF, points(get(varname), get(tr), pch=20))
    abline(h=0, lty=2)
}


# alpha at 0.2 


varname <- 'e_R'
with(DF, plot(get(varname),   e_tr005, pch=20, col=1))
with(DF, points(get(varname), e_tr01,  pch=20, col=2))
with(DF, points(get(varname), e_tr015, pch=20, col=3))
with(DF, points(get(varname), e_tr02,  pch=20, col=4))
with(DF, points(get(varname), e_tr025, pch=20, col=5))
with(DF, points(get(varname), e_tr03,  pch=20, col=6))
abline(h=0, lty=2)
}



