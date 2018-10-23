#### Assessing the impact of systematic errors in [Fe/H]
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

options(scipen=100000)
source('../scripts/utils.R') 
library(parallelMap)
num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
parallelStartMulticore(num_procs)

draw.grid = F

solar_Dnu = 135.1
solar_e_Dnu = 0.1
solar_Teff = 5772
solar_e_Teff = 0.8
solar_nu_max = 3090
solar_e_nu_max = 30

directory <- file.path('learn-final')

get_star <- function(filename) {
    if (!length(grep('^\\d+\\.dat$', basename(filename)))) return(NULL)
    KIC <- as.numeric(strsplit(basename(filename), '.dat')[[1]][1])
    if (KIC == 5774694) return(NULL)
    print(filename)
    
    covs <- read.table(filename, header=T)
    covs <- covs[complete.cases(covs),]
    cov.DF <- with(covs, 
        data.frame(KIC=KIC,
                   R=median(radius),    e_R=mad(radius),
                   M=median(M),         e_M=mad(M),
                   rho=median(density), e_rho=mad(density),
                   age=median(age),     e_age=mad(age)))
    
    if ('L' %in% names(covs))
        cov.DF <- cbind(cov.DF,
            with(covs, data.frame(L=median(L), e_L=mad(L))))
    
    perturb <- read.table(file.path('perturb', 'feh', 
        sub('\\.', '_perturb.', basename(filename))), header=1)
    perturb.DF <- with(perturb,
            data.frame(Teff=median(Teff), e_Teff=mad(Teff),
                       nu_max=median(nu_max), e_nu_max=mad(nu_max),
                       FeH=median(Fe.H),  e_FeH=mad(Fe.H)))
    
    if ('Dnu0' %in% names(perturb)) 
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(Dnu=median(Dnu0),  e_Dnu=mad(Dnu0))))
    
    if ('dnu02' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF, 
            with(perturb, data.frame(dnu=median(dnu02), e_dnu=mad(dnu02))))
    
    if ('L' %in% names(perturb))
        perturb.DF <- cbind(perturb.DF,
            with(perturb, data.frame(L=median(L), e_L=mad(L))))
    
    cbind(cov.DF, perturb.DF)
}

get_DF <- function(directory) do.call(plyr:::rbind.fill, 
    Map(get_star, list.files(directory, recursive=T, full.names=T)))

DF <- get_DF(file.path(directory))
DF <- DF[order(DF$KIC),]

scaling_R <- with(DF, 
    nu_max/solar_nu_max * (Dnu/solar_Dnu)**-2 * (Teff / solar_Teff)**0.5)
set.seed(0)
scaling_dR <- apply(do.call(rbind, Map(function(ii) with(DF, 
        (rnorm(nrow(DF), nu_max, e_nu_max) / 
            rnorm(1, solar_nu_max, solar_e_nu_max)) * 
        (rnorm(nrow(DF), Dnu, e_Dnu) / 
            rnorm(1, solar_Dnu, solar_e_Dnu))**-2 * 
        (rnorm(nrow(DF), Teff, e_Teff) / 
            rnorm(1, solar_Teff, solar_e_Teff))**0.5), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))
set.seed(0)
R.sigma <- apply(do.call(rbind, Map(function(ii) 
        rnorm(nrow(DF), DF$R, DF$e_R) 
        - 
        rnorm(nrow(DF), scaling_R, scaling_dR), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))

scaling_M <- with(DF, 
    (nu_max / solar_nu_max)**3 * 
       (Dnu / solar_Dnu)**-4 * 
      (Teff / solar_Teff)**(3/2))
set.seed(0)
scaling_dM <- apply(do.call(rbind, Map(function(ii) with(DF, 
        (rnorm(nrow(DF), nu_max, e_nu_max) / 
            rnorm(1, solar_nu_max, solar_e_nu_max))**3 * 
        (rnorm(nrow(DF), Dnu, e_Dnu) / 
            rnorm(1, solar_Dnu, solar_e_Dnu))**-4 * 
        (rnorm(nrow(DF), Teff, e_Teff) /
            rnorm(1, solar_Teff, solar_e_Teff))**(3/2)), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))
set.seed(0)
M.sigma <- apply(do.call(rbind, Map(function(ii) 
        rnorm(nrow(DF), DF$M, DF$e_M) 
        - 
        rnorm(nrow(DF), scaling_M, scaling_dM), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))

scaling_rho <- with(DF,  
       (Dnu / solar_Dnu)**2 * 
      (Teff / solar_Teff))
set.seed(0)
scaling_drho <- apply(do.call(rbind, Map(function(ii) with(DF, 
        (rnorm(nrow(DF), nu_max, e_nu_max) / 
            rnorm(1, solar_nu_max, solar_e_nu_max))**0 * 
        (rnorm(nrow(DF), Dnu, e_Dnu) / 
            rnorm(1, solar_Dnu, solar_e_Dnu))**2 * 
        (rnorm(nrow(DF), Teff, e_Teff) /
            rnorm(1, solar_Teff, solar_e_Teff))), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))
set.seed(0)
rho.sigma <- apply(do.call(rbind, Map(function(ii) 
        rnorm(nrow(DF), DF$rho, DF$e_rho) 
        - 
        rnorm(nrow(DF), scaling_rho, scaling_drho), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))


###### GAIA #######
gaia_filename <- 'gaia.dat'
con  <- file(gaia_filename, open = "r")
first_line <- readLines(con, n=1, warn=F)
close(con)
widths <- nchar(unlist(regmatches(first_line, gregexpr(".+?\\|", first_line))))
gaia.DF <- read.fwf(file=gaia_filename, widths=widths, skip=4)
colnames(gaia.DF) <- strsplit(
    gsub('\\s+', '', substring(first_line, 2)), '\\|')[[1]]

#LUM_VAL_GAIA
#ERROR_LUM_PERCENTILE_UPPER_GAIA
#ERROR_LUM_PERCENTILE_LOWER_GAIA
#RADIUS_VAL_GAIA
#ERROR_RADIUS_PERCENTILE_UPPER_GAIA
#ERROR_RADIUS_PERCENTILE_LOWER_GAIA

gaia.DF2 <- with(gaia.DF, data.frame(KIC=KIC_ID,
    L=LUM_VAL_GAIA, 
    e_L=ERROR_LUM_PERCENTILE_UPPER_GAIA-LUM_VAL_GAIA,
    R=RADIUS_VAL_GAIA,
    e_R=(ERROR_RADIUS_PERCENTILE_UPPER_GAIA
           -ERROR_RADIUS_PERCENTILE_LOWER_GAIA)/2))
gaia.DF2 <- gaia.DF2[!(gaia.DF2$L < -100),]

extras <- read.table('twomissingstars.csv', header=1, sep=',')
gaia.DF2 <- rbind(gaia.DF2,
    with(extras, data.frame(
        KIC=c(8379927, 10514430),
        L=lum_val,
        e_L=lum_percentile_upper-lum_percentile_lower,
        R=radius_val,
        e_R=(radius_percentile_upper - radius_percentile_lower)/2)))

new.DF <- merge(DF, gaia.DF2, by='KIC')
new.DF <- new.DF[order(new.DF$KIC),]

set.seed(0)
gaia.sigma.R <- apply(do.call(rbind, Map(function(ii) 
            rnorm(nrow(new.DF), new.DF$R.x, new.DF$e_R.x) 
            - 
            rnorm(nrow(new.DF), new.DF$R.y, new.DF$e_R.y), 
        ii=1:10000)),
    2, function(x) mad(x, na.rm=T))

set.seed(0)
gaia.sigma.L <- apply(do.call(rbind, Map(function(ii) 
        rnorm(nrow(new.DF), new.DF$L.x, new.DF$e_L.x) 
        - 
        rnorm(nrow(new.DF), new.DF$L.y, new.DF$e_L.y), 
    ii=1:10000)),
2, function(x) mad(x, na.rm=T))

#outliers <- unique(c( 
#    with(new.DF, KIC[abs(L.y-L.x)/gaia.sigma.L > 5
#                   | abs(R.y-R.x)/gaia.sigma.R > 5]),
#    with(DF,     KIC[abs(scaling_M-M)/M.sigma > 5 
#                   | abs(scaling_R-R)/R.sigma > 5])))
#out.cols <- adjustcolor(c('#111111', 'darkgreen', 'purple'), alpha.f=0.75)
#out.pch  <- 22:25
out.cols <- adjustcolor('#111111', alpha.f=0.75)
#out.cols <- adjustcolor(c('#161032', '#2D2D2A', 'darkgreen'), alpha.f=0.75)
out.pch <- 24

threshold <- 3
outliers.M <- with(DF, KIC[abs(scaling_M-M)/M.sigma > threshold]) 
outliers.R <- with(DF, KIC[abs(scaling_R-R)/R.sigma > threshold]) 
outliers.rho <- with(DF, KIC[abs(scaling_rho-rho)/rho.sigma > threshold]) 
outliers.R.gaia <- with(new.DF, KIC[abs(L.y-L.x)/gaia.sigma.L > threshold])
outliers.L <- with(new.DF, KIC[abs(R.y-R.x)/gaia.sigma.R > threshold])

outliers <- Reduce(union, list(
    intersect(outliers.M, outliers.R), 
    intersect(outliers.M, outliers.rho), 
    intersect(outliers.M, outliers.R.gaia), 
    intersect(outliers.M, outliers.L), 
    intersect(outliers.R, outliers.rho), 
    intersect(outliers.R, outliers.R.gaia), 
    intersect(outliers.R, outliers.L), 
    intersect(outliers.rho, outliers.R.gaia), 
    intersect(outliers.rho, outliers.L), 
    intersect(outliers.R.gaia, outliers.L)))
first <- outliers %in% outliers.M
second <- !outliers %in% outliers.M & outliers %in% outliers.R
third <- !outliers %in% outliers.M & !outliers %in% outliers.R &
    outliers %in% outliers.rho
fourth <- !outliers %in% outliers.M & !outliers %in% outliers.R &
    !outliers %in% outliers.rho & outliers %in% outliers.R.gaia
fifth <- !outliers %in% outliers.M & !outliers %in% outliers.R &
    !outliers %in% outliers.rho & !outliers %in% outliers.R.gaia
outliers <- c(
    outliers[first][order(sapply(outliers[first],
            function(x) DF$M[DF$KIC==x]))], 
    outliers[second][order(sapply(outliers[second],
            function(x) DF$R[DF$KIC==x]))],
    outliers[third][order(sapply(outliers[third],
            function(x) DF$rho[DF$KIC==x]))],
    outliers[fourth][order(sapply(outliers[fourth],
            function(x) new.DF$R.x[DF$KIC==x]))],
    outliers[fifth][order(sapply(outliers[fifth],
            function(x) new.DF$L.x[DF$KIC==x]))])
#outliers <- outliers[order(sapply(outliers, function(x) DF$M[DF$KIC==x]))]
#outlier.labs <- strsplit("abcdefghijklmnopqrstuvwxyz", '')[[1]]
#outlier.labs <- c(1:9, 'A', 'B', 'C', 'D', 'E')
#outlier.labs <- strsplit("ABCDEFGHIJKLMNOPQRSTUVWXYZ", '')[[1]]
outlier.labs <- 1:length(outliers)


formatter <- function(val, unc) {
    unc <- toString(signif(unc, 2))
    if (nchar(unc) == 1) {
        unc <- paste0(unc, '.0')
    }
    if (substr(unc, nchar(unc)-1, nchar(unc)-1) == '0' ||
        substr(unc, nchar(unc)-1, nchar(unc)-1) == '.' && 
        substr(unc, nchar(unc)-2, nchar(unc)-2) == '0') 
            unc <- paste0(unc, '0')
    dec_places <- nchar(strsplit(toString(unc), '\\.')[[1]][2])
    paste0('$',
        sprintf(val, fmt=paste0('%.', dec_places, 'f')),
        '\\pm',
        unc,
        '$'
    )
}

refs.DF <- read.table('KIC-planets.dat', header=T, stringsAsFactors=F)

## Print table 
sink(file='tab-parameters.tex')
refs <- c()
for (ii in c(1:nrow(DF))[order(as.numeric(DF$KIC))]) {#[1:10]) {
    star <- DF[ii,]
    planets <- if (star$KIC %in% refs.DF$KIC) {
        planet <- refs.DF[star$KIC == refs.DF$KIC,]
        output <- paste0(planet$Planets)
        planet.refs <- if (grepl(',', planet$References))
            strsplit(planet$References, ',')[[1]] else planet$References 
        
        planet.order <- c()
        for (planet.ref in planet.refs) {
            planet.order <- c(planet.order, if (planet.ref %in% refs) {
                which(planet.ref == refs)
            } else 99)
        }
        planet.refs <- planet.refs[order(planet.order)]
        
        for (planet.ref_i in 1:length(planet.refs)) {
            planet.ref <- planet.refs[planet.ref_i]
            output <- if (!planet.ref %in% refs) {
                refs <- c(refs, planet.ref)
                paste0(output, 
                    '\\footnotemark[', length(refs), ']',
                    if (planet.ref_i < length(planet.refs)) '$^,$' else '',
                    '\\footnotetext[', length(refs), ']{\\citealt{', 
                    planet.ref, '}}')
            } else {
                paste0(output, 
                    '\\footnotemark[', which(planet.ref == refs), ']')
            }
        }
        output
    } else ''
    with(star, cat(paste0(
        if (KIC == 8349582) '\\pagebreak ' else '',
        KIC, 
        if (KIC %in% outliers) {
            paste0('$^{\\text{\\hyperref[fig:scaling]{',
                outlier.labs[which(KIC == outliers)],
                '}}}$')
        } else '',
        '   \t&\t',
        formatter(age, e_age),
        '   \t&\t',
        formatter(M, e_M),
        '   \t&\t',
        formatter(R, e_R),
        '   \t&\t',
        formatter(rho, e_rho),
        '   \t&\t',
        formatter(L, e_L),
        '   \t&\t',
        planets,
        ' \\\\\n'
    )))
}
sink()

cat(paste0(
    formatter(mean(scaling_M - DF$M), sd(scaling_M - DF$M)), 
    '~$\\text{M}_\\odot$\n',
    formatter(mean(scaling_R - DF$R), sd(scaling_R - DF$R)), 
    '~$\\text{R}_\\odot$\n',
    formatter(mean(scaling_rho - DF$rho), sd(scaling_rho - DF$rho)), 
    '~$\\rho_\\odot$\n',
    formatter(mean(new.DF$R.x - new.DF$R.y), sd(new.DF$R.x - new.DF$R.y)), 
    '~$\\text{R}_\\odot$\n',
    formatter(mean(new.DF$L.x - new.DF$L.y), sd(new.DF$L.x - new.DF$L.y)), 
    '~$\\text{L}_\\odot$\n'
))


plot_R <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5, cex=text.cex)
    
    xlim <- c(0.4, 2.8)
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    # draw grid
    if (draw.grid) {
    for (xtick in xticks) 
        for (ytick in pretty(ylim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    
    for (ytick in yticks) 
        for (xtick in pretty(xlim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    }
    
    abline(a=0, b=1, lwd=1.5, lty=2)
    
    selection <- !DF$KIC %in% outliers.R
    with(DF[selection,], 
        points(R, scaling_R[selection], 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.R
    with(DF[selection,], 
        points(R, scaling_R[selection], 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.R) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        if (which(outlier == outliers) == 6) {
            text(0.992*DF$R[index], 
                 1*scaling_R[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        } else {
            text(0.992*DF$R[index], 0.985*scaling_R[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    x.mid <- 0.95 * xlim[2]#8
    y.mid <- 0.085 * (xlim[2] - xlim[1]) + xlim[1]#0.7
    
    segments(x.mid+median(DF$e_R), y.mid, 
             x.mid-median(DF$e_R), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(scaling_dR, na.rm=T), 
             x.mid, y.mid-median(scaling_dR, na.rm=T), 
        lwd=1, col='black')
    
    ## solar symbol
    points(1, 1, pch=20, cex=0.75, lwd=1.5)
    points(1, 1, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=2, labels=T, lwd.ticks=par()$lwd)
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Scaling Radius', 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    ## lower panel
    
    ylim <- c(-0.3, 0.3)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    abline(a=0, b=0, lwd=1.5, lty=2)
    
    for (ii in 1:3)
        rect(xlim[1],  ii*median(R.sigma, na.rm=T), 
             xlim[2], -ii*median(R.sigma, na.rm=T),
            col=adjustcolor('black', alpha.f=0.15), border=NA)
    
    selection <- !DF$KIC %in% outliers.R
    with(DF[selection,], 
        points(R, scaling_R[selection] - R, 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.R
    with(DF[selection,], 
        points(R, scaling_R[selection] - R, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.R) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        text(0.992*DF$R[index], 0.985*(scaling_R[index] - DF$R[index]), 
            labels=outlier.labs[which(outlier == outliers)],
            cex=0.7*text.cex, pos=3)
    }
    par(family=font)
    
    ## solar symbol
    points(1, 0, pch=20, cex=0.75,   lwd=1.5)
    points(1, 0, pch=1,  cex=1.5, lwd=1.5)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex,
        family=font, majorn=2, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Residuals', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('SPI Radius', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_R, 'R', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,#1.3,
    paper_pdf_width=4.17309*0.9,#1.1,
    cex.paper=0.75) 




plot_M <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5, cex=text.cex)
    
    xlim <- c(0.5, 1.9)
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)#, n=5)
    yticks <- xticks
    
    # draw grid 
    if (draw.grid) {
    for (xtick in xticks) 
        for (ytick in pretty(ylim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    
    for (ytick in yticks) 
        for (xtick in pretty(xlim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    }
    
    abline(a=0, b=1, lwd=1.5, lty=2)
    
    #scaling_M <- with(DF, 
    #    (nu_max / solar_nu_max)**3 * 
    #       (Dnu / solar_Dnu)**-4 * 
    #      (Teff / solar_Teff)**(3/2))
    
    #points(DF$M, scaling_M,
    #    #pch=21, 
    #    pch=sapply(DF$KIC, 
    #        function(x) if (!x %in% outliers.M) 21 
    #            else out.pch[which(x == outliers.M) %% length(out.pch) + 1]), 
    #    cex=1.1, 
    #    bg=sapply(DF$KIC, function (x) if (x %in% outliers.M) 
    #        out.cols[which(x == outliers.M) %% length(out.cols) + 1] 
    #        else adjustcolor(orange, alpha.f=0.75)), 
    #    lwd=0.75, col='white')
    
    selection <- !DF$KIC %in% outliers.M
    with(DF[selection,], 
        points(M, scaling_M[selection], 
        pch=21, cex=1.1, 
        bg=adjustcolor(orange, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.M
    with(DF[selection,], 
        points(M, scaling_M[selection], 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.M) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        text(0.992*DF$M[index], 0.985*scaling_M[index], 
            labels=outlier.labs[which(outlier == outliers)],
            cex=0.7*text.cex, pos=3)
    }
    par(family=font)
    
    #scaling_dM <- with(DF, 
    #    #(e_nu_max/3090)**3 * (e_Dnu/135)**-4 * (e_Teff / 5777)**(3/2))
    #    abs((nu_max/3090)**2 * 3 * (e_nu_max/3090))
    set.seed(0)
    scaling_dM <- apply(do.call(rbind, Map(function(ii) with(DF, 
                (rnorm(nrow(DF), nu_max, e_nu_max) / 
                    rnorm(1, solar_nu_max, solar_e_nu_max))**3 * 
                (rnorm(nrow(DF), Dnu, e_Dnu) / 
                    rnorm(1, solar_Dnu, solar_e_Dnu))**-4 * 
                (rnorm(nrow(DF), Teff, e_Teff) /
                    rnorm(1, solar_Teff, solar_e_Teff))**(3/2)), 
            ii=1:1000)),
        2, function(x) mad(x, na.rm=T))
    
    
    x.mid <- 0.95 * xlim[2]#8
    y.mid <- 0.085 * (xlim[2] - xlim[1]) + xlim[1]#0.7
    
    segments(x.mid+median(DF$e_M), y.mid, 
             x.mid-median(DF$e_M), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(scaling_dM, na.rm=T), 
             x.mid, y.mid-median(scaling_dM, na.rm=T), 
        lwd=1, col='black')
    
    ## solar symbol
    points(1, 1, pch=20, cex=0.75, lwd=1.5)
    points(1, 1, pch=1,  cex=1.5,  lwd=1.5)
    
    #magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=text.cex, 
    #    family=font, majorn=3, labels=T, lwd.ticks=par()$lwd)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=4, labels=F, lwd.ticks=par()$lwd)
    
    #axis(1, yticks, labels=F, tick=F, 
    #    cex.axis=text.cex, tcl=0, las=1, 
    #    mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, pretty(xlim, n=3)[-1], labels=pretty(xlim, n=3)[-1], tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Scaling Mass', 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    
    
    
    ## lower panel
    
    ylim <- c(-0.7, 0.7)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    #set.seed(0)
    #sigma <- apply(do.call(rbind, Map(function(ii) 
    #            rnorm(nrow(DF), DF$M, DF$e_M) 
    #            - 
    #            rnorm(nrow(DF), scaling_M, scaling_dM), 
    #        ii=1:1000)),
    #    2, function(x) mad(x, na.rm=T))
    
    for (ii in 1:3)
        rect(xlim[1],  ii*median(M.sigma, na.rm=T), 
             xlim[2], -ii*median(M.sigma, na.rm=T),
            col=adjustcolor('black', alpha.f=0.15), border=NA)
    
    abline(a=0, b=0, lwd=1.5, lty=2)
    
    #segments(x.mid+median(DF$e_M), -0.35, 
    #         x.mid-median(DF$e_M), -0.35,
    #    lwd=1.5, col='black')
    #segments(x.mid, -0.35+median(scaling_dM, na.rm=T), 
    #         x.mid, -0.35-median(scaling_dM, na.rm=T), 
    #    lwd=1.5, col='black')
    
    #points(DF$M, scaling_M - DF$M,
    #    #pch=21, 
    #    pch=sapply(DF$KIC, 
    #        function(x) if (!x %in% outliers.M) 21 
    #            else out.pch[which(x == outliers.M) %% length(out.pch) + 1]), 
    #    cex=1.1, 
    #    bg=sapply(DF$KIC, function (x) if (x %in% outliers.M) 
    #        out.cols[which(x == outliers.M) %% length(out.cols) + 1] 
    #        else adjustcolor(orange, alpha.f=0.75)), 
    #    lwd=0.75, col='white')
    
    selection <- !DF$KIC %in% outliers.M
    with(DF[selection,], 
        points(M, scaling_M[selection] - M, 
        pch=21, cex=1.1, 
        bg=adjustcolor(orange, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.M
    with(DF[selection,], 
        points(M, scaling_M[selection] - M, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.M) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        text(0.992*DF$M[index], 0.985*(scaling_M[index] - DF$M[index]), 
            labels=outlier.labs[which(outlier == outliers)],
            cex=0.7*text.cex, pos=3)
    }
    par(family=font)
    
    ## solar symbol
    points(1, 0, pch=20, cex=0.75, lwd=1.5)
    points(1, 0, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex, 
        family=font, majorn=4, labels=T, lwd.ticks=par()$lwd,)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Residuals', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('SPI Mass', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_M, 'M', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,#1.3,
    paper_pdf_width=4.17309*0.9,#1.1,
    cex.paper=0.75) 




plot_rho <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5, cex=text.cex)
    
    xlim <- c(0, 2)
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim, n=6)
    yticks <- xticks #pretty(ylim)
    
    # draw grid 
    if (draw.grid) {
    for (xtick in xticks) 
        for (ytick in pretty(ylim, n=18))
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    
    for (ytick in yticks) 
        for (xtick in pretty(xlim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    }
    
    abline(a=0, b=1, lwd=1.5, lty=2)
    
    #points(DF$rho, scaling_rho,
    #    pch=sapply(DF$KIC, 
    #        function(x) if (!x %in% outliers.rho) 21 
    #            else out.pch[which(x == outliers.rho) %% length(out.pch) + 1]), 
    #    cex=1.1, 
    #    bg=sapply(DF$KIC, function (x) if (x %in% outliers.rho) 
    #        out.cols[which(x == outliers.rho) %% length(out.cols) + 1] 
    #        else adjustcolor(blue, alpha.f=0.75)), 
    #    lwd=0.75, col='white')
    
    selection <- !DF$KIC %in% outliers.rho
    with(DF[selection,], 
        points(rho, scaling_rho[selection], 
        pch=21, cex=1.1, 
        bg=adjustcolor(blue, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.rho
    with(DF[selection,], 
        points(rho, scaling_rho[selection], 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.rho) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        text(0.992*DF$rho[index], 0.985*scaling_rho[index], 
            labels=outlier.labs[which(outlier == outliers)],
            cex=0.7*text.cex, pos=3)
        if (which(outlier == outliers) == 1) {
            text(0.992*DF$rho[index], 
                1.05*scaling_rho[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else {
            text(0.992*DF$rho[index], 0.985*scaling_rho[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    x.mid <- 0.95 * xlim[2]#8
    y.mid <- 0.085 * (xlim[2] - xlim[1]) + xlim[1]#0.7
    
    segments(x.mid+median(DF$e_rho), y.mid, 
             x.mid-median(DF$e_rho), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(scaling_drho, na.rm=T), 
             x.mid, y.mid-median(scaling_drho, na.rm=T), 
        lwd=1, col='black')
    
    ## solar symbol
    points(1, 1, pch=20, cex=0.75, lwd=1.5)
    points(1, 1, pch=1,  cex=1.5,  lwd=1.5)
    
    par(xpd=NA)
    rect(xlim[2], ylim[1]*0.1, xlim[2]*1.1, ylim[2]*1.1, col='white', border=NA)
    rect(xlim[1]*0.9, ylim[1], xlim[2]*1.1, ylim[1]*0.9, col='white', border=NA)
    par(xpd=F)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=4, labels=F, lwd.ticks=par()$lwd)
    
    axis(2, yticks[-c(1, length(yticks))], 
        labels=yticks[-c(1, length(yticks))], 
        tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Scaling Density', 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    ## lower panel
    
    ylim <- c(-0.3, 0.3)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    #xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    for (ii in 1:3)
        rect(xlim[1],  ii*median(rho.sigma, na.rm=T), 
             xlim[2], -ii*median(rho.sigma, na.rm=T),
            col=adjustcolor('black', alpha.f=0.15), border=NA)
    
    abline(a=0, b=0, lwd=1.5, lty=2)
    
    selection <- !DF$KIC %in% outliers.rho
    with(DF[selection,], 
        points(rho, scaling_rho[selection] - rho, 
        pch=21, cex=1.1, 
        bg=adjustcolor(blue, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- DF$KIC %in% outliers.rho
    with(DF[selection,], 
        points(rho, scaling_rho[selection] - rho, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.rho) {
        if (!outlier %in% outliers || !outlier %in% DF$KIC) next
        index <- which(outlier == DF$KIC)
        if (which(outlier == outliers) == 1) {
            text(0.992*DF$rho[index], 
                0.85*(scaling_rho[index] - DF$rho[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else {
            text(0.992*DF$rho[index], 0.985*(scaling_rho[index] - DF$rho[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    ## solar symbol
    points(1, 0, pch=20, cex=0.75, lwd=1.5)
    points(1, 0, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex, 
        family=font, majorn=4, labels=T, lwd.ticks=par()$lwd,)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, yticks, labels=F, tick=F, 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Residuals', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('SPI Density', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_rho, 'rho', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,#1.3,
    paper_pdf_width=4.17309*0.9,#1.1,
    cex.paper=0.75) 







########## GAIA ###########
#outliers <- unique(outliers.R.gaia, outliers.L)

plot_R <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5, cex=text.cex)
    
    xlim <- c(0.4, 2.8)
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    abline(a=0, b=1, lwd=1.5, lty=2)
    
    selection <- !new.DF$KIC %in% outliers.R.gaia
    with(new.DF[selection,], 
        points(R.x, R.y, 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- new.DF$KIC %in% outliers.R.gaia
    with(new.DF[selection,], 
        points(R.x, R.y, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.R.gaia) {
        if (!outlier %in% outliers || !outlier %in% new.DF$KIC) next
        index <- which(outlier == new.DF$KIC)
        #if (which(outlier == outliers) == 6) {
        #    text(0.992*new.DF$R.x[index], 
        #         1.02*new.DF$R.y[index], 
        #        labels=outlier.labs[which(outlier == outliers)],
        #        cex=0.7*text.cex, pos=1)
        #} else {
        if (which(outlier == outliers) == 7) {
            text(0.98*new.DF$R.x[index], 
                 0.985*new.DF$R.y[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        } else if (which(outlier == outliers) == 10) {
            text(1.01*new.DF$R.x[index], 
                 0.985*new.DF$R.y[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        } else {
            text(0.992*new.DF$R.x[index], 
                 0.985*new.DF$R.y[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    x.mid <- 0.95 * xlim[2]#8
    y.mid <- 0.085 * (xlim[2] - xlim[1]) + xlim[1]#0.7
    
    segments(x.mid+median(new.DF$e_R.x), y.mid, 
             x.mid-median(new.DF$e_R.x), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(new.DF$e_R.y), 
             x.mid, y.mid-median(new.DF$e_R.y), 
        lwd=1, col='black')
    
    ## solar symbol
    points(1, 1, pch=20, cex=0.75, lwd=1.5)
    points(1, 1, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=2, labels=T, lwd.ticks=par()$lwd)
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Gaia Radius', 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    
    
    
    ## lower panel
    
    ylim <- c(-0.7, 0.7)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    abline(a=0, b=0, lwd=1.5, lty=2)
    
    
    
    for (ii in 1:3)
        rect(xlim[1],  ii*median(gaia.sigma.R, na.rm=T), 
             xlim[2], -ii*median(gaia.sigma.R, na.rm=T),
            col=adjustcolor('black', alpha.f=0.15), border=NA)
    
    selection <- !new.DF$KIC %in% outliers.R.gaia
    with(new.DF[selection,], 
        points(R.x, R.y - R.x, 
        pch=21, cex=1.1, 
        bg=adjustcolor(red, alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- new.DF$KIC %in% outliers.R.gaia
    with(new.DF[selection,], 
        points(R.x, R.y - R.x, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    par(family="Helvetica")
    for (outlier in outliers.R.gaia) {
        if (!outlier %in% outliers || !outlier %in% new.DF$KIC) next
        index <- which(outlier == new.DF$KIC)
        if (which(outlier == outliers) == 7) {
            text(0.98*new.DF$R.x[index], 
                 1.1*(new.DF$R.y[index] - new.DF$R.x[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else if (which(outlier == outliers) == 10) {
            text(1.01*new.DF$R.x[index], 
                 1.1*(new.DF$R.y[index] - new.DF$R.x[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else {
            text(0.992*new.DF$R.x[index], 
                 0.985*(new.DF$R.y[index] - new.DF$R.x[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    ## solar symbol
    points(1, 0, pch=20, cex=0.75,   lwd=1.5)
    points(1, 0, pch=1,  cex=1.5, lwd=1.5)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex, 
        cex=text.cex,
        family=font, majorn=2, labels=T, lwd.ticks=par()$lwd,
        cex=text.cex,)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Residuals', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('SPI Radius', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_R, 'R-gaia', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    paper_pdf_width=4.17309*0.9,
    cex.paper=0.75) 




plot_L <- function(..., 
        text.cex=1, mgp=utils.mgp, mar=utils.mar, font=utils.font, 
        tcl=utils.tcl) {
    
    layout(rbind(1,1,2), respect=F)
    
    par(oma=mar+c(0.2, -0.6, -0.3, -0.3), mar=c(0,0,0,0),
        mgp=mgp+c(0, 0.4, 0), lwd=1.5, cex=text.cex)
    
    xlim <- c(0, 9)
    ylim <- xlim
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    # draw grid 
    if (draw.grid) {
    for (xtick in xticks) 
        for (ytick in pretty(ylim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    
    for (ytick in yticks) 
        for (xtick in pretty(xlim, n=18)) 
            points(xtick, ytick, pch=19, lwd=0.1, cex=0.2, 
                bg='darkgray', col='darkgray')
    }
    
    abline(a=0, b=1, lwd=1.5, lty=2)
    
    #points(new.DF$L.x, new.DF$L.y,
    #    #pch=21, 
    #    pch=sapply(new.DF$KIC, 
    #        function(x) if (!x %in% outliers.L) 21 
    #            else out.pch[which(x == outliers.L) %% length(out.pch) + 1]), 
    #    cex=1.1, 
    #    bg=sapply(new.DF$KIC, function (x) if (x %in% outliers.L) 
    #        out.cols[which(x == outliers.L) %% length(out.cols) + 1] 
    #        else adjustcolor('darkblue', alpha.f=0.75)), 
    #    lwd=0.75, col='white')
    
    selection <- !new.DF$KIC %in% outliers.L
    with(new.DF[selection,], 
        points(L.x, L.y, 
        pch=21, cex=1.1, 
        bg=adjustcolor('#509E6A', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- new.DF$KIC %in% outliers.L
    with(new.DF[selection,], 
        points(L.x, L.y, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    
    par(family="Helvetica")
    for (outlier in outliers.L) {
        if (!outlier %in% outliers || !outlier %in% new.DF$KIC) next
        index <- which(outlier == new.DF$KIC)
        #text(0.992*new.DF$L.x[index], 0.985*new.DF$L.y[index], 
        #    labels=outlier.labs[which(outlier == outliers)],
        #    cex=0.7*text.cex, pos=3)
        if (which(outlier == outliers) == 7) {
            text(0.992*new.DF$L.x[index], 
                 1.2*new.DF$L.y[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else {
            text(0.992*new.DF$L.x[index], 
                 0.985*new.DF$L.y[index], 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    x.mid <- 0.95 * xlim[2]#8
    y.mid <- 0.085 * (xlim[2] - xlim[1]) + xlim[1]#0.7
    
    segments(x.mid+median(new.DF$e_L.x), y.mid, 
             x.mid-median(new.DF$e_L.x), y.mid,
        lwd=1, col='black')
    segments(x.mid, y.mid+median(new.DF$e_L.y), 
             x.mid, y.mid-median(new.DF$e_L.y), 
        lwd=1, col='black')
    
    ## solar symbol
    #points(1, 1, pch=20, cex=0.75, lwd=2.4, col='white')
    #points(1, 1, pch=1,  cex=1.1,  lwd=2.4, col='white')
    points(1, 1, pch=20, cex=0.75, lwd=1.5)
    points(1, 1, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=4, labels=F, lwd.ticks=par()$lwd)
    
    axis(2, xticks[-1], tick=F, #labels=xticks[-1], 
        cex.axis=1.3*text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Gaia Luminosity', 2, 2.5, outer=F, las=0, cex=text.cex)
    
    
    
    
    
    ## lower panel
    
    ylim <- c(-3, 3)
    
    plot(NA, axes=F, 
        xaxs='i', yaxs='i', 
        xlab="", ylab="", 
        xlim=xlim, ylim=ylim)
    
    xticks <- pretty(xlim)
    yticks <- pretty(ylim)
    
    #set.seed(0)
    #sigma <- apply(do.call(rbind, Map(function(ii) 
    #            rnorm(nrow(new.DF), new.DF$L.x, new.DF$e_L.x) 
    #            - 
    #            rnorm(nrow(new.DF), new.DF$L.y, new.DF$e_L.y), 
    #        ii=1:1000)),
    #    2, function(x) mad(x, na.rm=T))
    
    for (ii in 1:3)
        rect(xlim[1],  ii*median(gaia.sigma.L), 
             xlim[2], -ii*median(gaia.sigma.L),
            col=adjustcolor('black', alpha.f=0.15), border=NA)
    
    abline(a=0, b=0, lwd=1.5, lty=2)
    
    #segments(x.mid+median(DF$e_M), -0.35, 
    #         x.mid-median(DF$e_M), -0.35,
    #    lwd=1.5, col='black')
    #segments(x.mid, -0.35+median(scaling_dM, na.rm=T), 
    #         x.mid, -0.35-median(scaling_dM, na.rm=T), 
    #    lwd=1.5, col='black')
    
    #points(new.DF$L.x, new.DF$L.y - new.DF$L.x,
    #    #pch=21, 
    #    pch=sapply(new.DF$KIC, 
    #        function(x) if (!x %in% outliers.L) 21 
    #            else out.pch[which(x == outliers.L) %% length(out.pch) + 1]), 
    #    cex=1.1, 
    #    bg=sapply(new.DF$KIC, function (x) if (x %in% outliers.L) 
    #        out.cols[which(x == outliers.L) %% length(out.cols) + 1] 
    #        else adjustcolor('darkblue', alpha.f=0.75)), 
    #    lwd=0.75, col='white')
    
    selection <- !new.DF$KIC %in% outliers.L
    with(new.DF[selection,], 
        points(L.x, L.y - L.x, 
        pch=21, cex=1.1, 
        bg=adjustcolor('#509E6A', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    selection <- new.DF$KIC %in% outliers.L
    with(new.DF[selection,], 
        points(L.x, L.y - L.x, 
        pch=24, cex=1.1, 
        bg=adjustcolor('#111111', alpha.f=0.75), 
        lwd=0.75, col='white'))
    
    
    par(family="Helvetica")
    for (outlier in outliers.L) {
        if (!outlier %in% outliers || !outlier %in% new.DF$KIC) next
        index <- which(outlier == new.DF$KIC)
        if (which(outlier == outliers) == 7) {
            text(0.992*new.DF$L.x[index], 
                 1.1*(new.DF$L.y[index] - new.DF$L.x[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=1)
        } else {
            text(0.992*new.DF$L.x[index], 
                 0.985*(new.DF$L.y[index] - new.DF$L.x[index]), 
                labels=outlier.labs[which(outlier == outliers)],
                cex=0.7*text.cex, pos=3)
        }
    }
    par(family=font)
    
    ## solar symbol
    #points(1, 0, pch=20, cex=0.75, lwd=2.4, col='white')
    #points(1, 0, pch=1,  cex=1.1,  lwd=2.4, col='white')
    points(1, 0, pch=20, cex=0.75, lwd=1.5)
    points(1, 0, pch=1,  cex=1.5,  lwd=1.5)
    
    magaxis(1, tcl=tcl, mgp=mgp+c(0, 0.25, 0), cex.axis=1.3*text.cex, 
        family=font, majorn=4, labels=T, lwd.ticks=par()$lwd,)
    magaxis(2, tcl=tcl, mgp=mgp+c(0, 0.35, 0), cex.axis=1.3*text.cex, 
        family=font, las=1, majorn=3, labels=T, lwd.ticks=par()$lwd)
    
    axis(1, yticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.25, 0), lwd=1.5) 
    
    axis(2, xticks, labels=F, tick=F, 
        cex.axis=text.cex, tcl=0, las=1, 
        mgp=mgp+c(0, 0.35, 0), lwd=1.5) 
    
    box(lwd=par()$lwd)
    
    mtext('Residuals', 2, 2.5, outer=F, las=0, cex=text.cex)
    mtext('SPI Luminosity', 1, 2, outer=F, cex=text.cex)
}

make_plots(plot_L, 'L', 
    filepath=file.path('plots', 'feh'), 
    slides=F, make_png=F, wide=F, tall=F,
    paper_pdf_height=4.17309*1.1,
    paper_pdf_width=4.17309*0.9,
    cex.paper=0.75) 


