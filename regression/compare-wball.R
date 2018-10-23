source(file.path('..', 'scripts', 'utils.R'))

plot_comparison <- function(x, y, e_xu, e_xl, e_y, mod, color, lim, cutoff, 
        xlab, ylab, laboff=0, leg=NULL, ...,
        font=utils.font, mgp=utils.mgp, mar=utils.mar, tcl=utils.tcl) {
    par(mar=mar+c(0.5,-0.2,-0.3,-0.3), mgp=mgp+c(0, 0.4, 0), pty="s", lwd=1.66)
    #lim <- range(x-e_x, x+e_x, y-e_y, y+e_y)
    plot(NA, axes=F, xaxs='i', yaxs='i', 
        pch=20, cex=1, 
        xlim=lim, ylim=lim,
        xlab="", ylab="")
    segments(0,0,20,20, lty=2, lwd=par()$lwd)
    abline(v=cutoff, lty=3, lwd=2)
    abline(h=cutoff, lty=3, lwd=2)
    
    #arrows(x, y+e_y, x, y-e_y, lwd=2.5, col=1, angle=90, length=0.01, code=3)
    #arrows(x+e_xu, y, x+e_xl, y, lwd=2.5, col=1, angle=90, length=0.01, code=3)
    arrows(x, y+e_y, x, y-e_y, lwd=2, col=color, angle=90, length=0.01, code=3)
    arrows(x+e_xu, y, x+e_xl, y, lwd=2, col=color, angle=90, length=0.01, code=3)
    #points(x, y, pch=20, cex=0.6)#, col=color)
    points(x, y, pch=20, cex=0.66, lwd=1)
    #points(x, mod, pch=1, cex=1)
    
    if (!is.null(leg)) {
        par(family='Helvetica')# LT Std Light')
        legend('topleft', bty='n', inset=c(-0.12, -0.04), #pch=NULL, lty=NULL,
            legend=c(leg), cex=text.cex)
        par(family=font)
    }
    
    #magaxis(1:4, tcl=0, labels=F)
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, majorn=3, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, lwd.ticks=par()$lwd)
    
    box(lwd=par()$lwd)
    
    par(mgp=mgp+c(0.6+laboff, 0, 0))
    title(xlab=as.expression(xlab))
    par(mgp=mgp+c(0.8-laboff, 0, 0))
    title(ylab=as.expression(ylab))
}

dirs <- c('inv', 'inv-ov', 'inv-D', 'inv-Dov')
#description <- c('No diffusion, no overshoot', 
#                 'No diffusion, overshoot',
#                 'Diffusion, no overshoot',
#                 'Diffusion, overshoot')
description <- c('Vanilla', 'Overshoot', 'Diffusion', 'Both')
pipelines <- c('BASTA') #c('GOE', 'BASTA')
for (pipeline in pipelines) {
    res <- read.table(file.path('data', 
            paste0(pipeline, '_LEGACY_TableC4-1.txt')), 
        header=1)
    res <- with(res, 
        data.frame(Name=KIC, M=Mass, e_M=(sMassP-sMassM)/2, 
                   e_Mu=sMassP, e_Ml=sMassM,
                   age=Age, e_age=(sAgeP-sAgeM)/2, 
                   e_ageu=sAgeP, e_agel=sAgeM))
    
    for (dir_i in 1:length(dirs)) {
        name <- dirs[dir_i]
        inv <- read.table(file.path(paste0('learn-', name), 
                paste0('tables-', name), 'inversions.dat'), 
            header=1, fill=T)
        DF <- merge(res, inv, by='Name')
        
        #modes <- read.table(file.path(paste0('learn-', name), 
        #        paste0('tables-', name), 'inversions_modes.dat'), 
        #    header=1)
        #DF <- merge(DF, modes, by='Name')
        
        color <- c('darkgray', orange, blue, red)[dir_i]
        
        x <- DF$M.x
        y <- DF$M.y
        e_xu <- DF$e_Mu
        e_xl <- DF$e_Ml
        e_y <- DF$e_M.y
        mod <- DF$M.y #DF$M
        make_plots(plot_comparison, paste0(name, '-mass'), 
            filepath=file.path('plots', paste0('compare-', pipeline)),
            x=x, y=y, e_xu=e_xu, e_xl=e_xl, e_y=e_y, mod=mod,
            color=color, 
            lim=c(0.6, 1.99141), #c(0.5, 2.00736), 
            cutoff=1.8, leg=description[dir_i],
            xlab="BASTA Mass", #bquote("GOE Mass"~M/M['ʘ']), 
            ylab="SPI Mass", #bquote("ML Mass"~M/M['ʘ']),
            cex.paper=1.16, wide=F, tall=F, make_png=F,
            paper_pdf_height=5.5,
            paper_pdf_width=5.5)
            #paper_pdf_height=6.97522)
            #paper_pdf_width=5.5, paper_pdf_height=5.5)
            #paper_pdf_width=4.17309)#, paper_pdf_height=6.97522)
        
        x <- DF$age.x
        y <- DF$age.y
        e_xu <- DF$e_ageu
        e_xl <- DF$e_agel
        e_y <- DF$e_age.y
        mod <- DF$age.y #DF$age
        make_plots(plot_comparison, paste0(name, '-age'), 
            filepath=file.path('plots', paste0('compare-', pipeline)),
            x=x, y=y, e_xu=e_xu, e_xl=e_xl, e_y=e_y, mod=mod,
            color=color, 
            lim=c(0, 16), cutoff=13.799, laboff=0.2,
            leg=description[dir_i],
            xlab="BASTA Age", #bquote("GOE Age"~tau/Gyr), 
            ylab="SPI Age", #bquote("ML Age"~tau/Gyr),
            cex.paper=1.16, wide=F, tall=F, make_png=F,
            paper_pdf_height=5.5,
            paper_pdf_width=5.5)
            #paper_pdf_height=6.97522)
            #paper_pdf_width=5.7, paper_pdf_height=5.7)
            #paper_pdf_width=4.17309)#, paper_pdf_height=6.97522)
    }
}


if (F) {
wball <- read.table(file.path('data', 'GOE_LEGACY_TableC4-1.txt'), header=1)
wball <- with(wball, 
    data.frame(Name=KIC, M=Mass, e_M=(sMassP-sMassM)/2, 
               e_Mu=sMassP, e_Ml=sMassM,
               age=Age, e_age=(sAgeP-sAgeM)/2, 
               e_ageu=sAgeP, e_agel=sAgeM))

basta <- read.table(file.path('data', 'BASTA_LEGACY_TableC4-1.txt'), header=1)
basta <- with(basta, 
    data.frame(Name=KIC, M=Mass, e_M=(sMassP-sMassM)/2, 
               e_Mu=sMassP, e_Ml=sMassM,
               age=Age, e_age=(sAgeP-sAgeM)/2, 
               e_ageu=sAgeP, e_agel=sAgeM))
}


if (F) {
x <- DF$M.x
y <- DF$M.y
e_x <- DF$e_M.x
e_y <- DF$e_M.y
lim <- range(x-e_x, x+e_x, y-e_y, y+e_y)
plot(x, y, 
    pch=20, cex=1, 
    xlim=lim, ylim=lim,
    xlab="GOE Mass",
    ylab="ML Mass")
segments(x, y+e_y, x, y-e_y)
segments(x+e_x, y, x-e_x, y)
segments(0,0,20,20, lty=2)
dev.off()

x <- DF$age.x
y <- DF$age.y
e_x <- DF$e_age.x
e_y <- DF$e_age.y
lim <- range(x-e_x, x+e_x, y-e_y, y+e_y)
plot(x, y, 
    pch=20, cex=1, 
    xlim=lim, ylim=lim,
    xlab="GOE Age",
    ylab="ML Age")
segments(x, y+e_y, x, y-e_y)
segments(x+e_x, y, x-e_x, y)
segments(0,0,20,20, lty=2)
dev.off()
}
