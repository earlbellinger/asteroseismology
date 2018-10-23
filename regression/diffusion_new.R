#### New diffusion diagram 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 
#library(emojifont)
library(parallelMap)
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

directory <- file.path('feh', 'learn-feh')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    print(filename)
    
    covs <- read.table(filename, header=T)
    with(covs, 
        data.frame(KIC=KIC,
                   M=median(M),         e_M=mad(M),
                   alpha=median(alpha), e_alpha=mad(alpha),
                   age=median(age),     e_age=mad(age),
                   diffusion=median(diffusion), e_diffusion=mad(diffusion)))
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(file.path(directory))

track.col <- "#222222"#'black'
point.col <- 'darkblue'
point.border <- 'white'#'black'

plot_diffusion <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    par(mar=mar+c(0.3, -0.5, 0.1, 0.1), lwd=1.5, las=1, cex.axis=text.cex)
    
    xlim <- c(0.7, 1.8)
    ylim <- c(0.01, 10)
    
    plot(NA, axes=F, log='y', 
        xaxs='i', yaxs='i', 
        xlim=xlim, ylim=ylim, 
        xlab="", ylab="")
    
    with(DF, segments(M-e_M, diffusion, M+e_M, diffusion, lwd=0.5, lty=1))
    with(DF, segments(M, diffusion-e_diffusion, 
                      M, diffusion+e_diffusion, lwd=0.5, lty=1))
    
    with(DF,
        points(M, diffusion, 
            pch=21, 
            cex=1.1, 
            bg=adjustcolor(point.col, alpha.f=0.75),
            lwd=0.75, 
            col=point.border))
    
    ## solar symbol
    points(1, 1, pch=20, cex=1.1/2,   lwd=1.5)
    points(1, 1, pch=1,  cex=1.1, lwd=1.5)
    
    #par(xpd=NA)
    #rect(xlim[2], ylim[1]*0.1, xlim[2]-100, ylim[2]*10, col='white', border=NA)
    #rect(xlim[1]*1.1, ylim[1], xlim[2]*0.9, ylim[1]*0.9, col='white', border=NA)
    #par(xpd=F)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.4, 0), cex.axis=text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    box(lwd=par()$lwd)
    
    mtext('Diffusion', 2, 2, outer=F, las=0, cex=text.cex)
    mtext('Mass', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_diffusion, 
    paste0('diffusion'),
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, 
    wide=F, tall=F, 
    paper_pdf_height=4.17309*1.385,
    cex.paper=0.95,
    use.cairo=T, font='Palatino Linotype') 


