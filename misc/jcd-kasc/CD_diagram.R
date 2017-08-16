#### Plot JCD Diagram for KASC proceedings
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Department of Astronomy, Yale University;
#### Stellar Ages & Galactic Evolution Group,
#### Max-Planck-Institut fur Sonnensystemforschung

library(magicaxis)
library(parallelMap)
#library(akima)
#library(RColorBrewer)
source(file.path('..', '..', 'scripts', 'seismology.R'))

separations <- function(filename) {
    print(filename)
    prof.filename <- paste0(strsplit(filename, '-')[[1]][1], '.data')
    prof <- read.table(prof.filename, header=1, skip=1, nrow=1)

    nu_max <- scaling_nu_max(prof$photosphere_r, prof$star_mass, prof$Teff)

    freqs <- parse_freqs(filename, gyre=T)
    #freqs <- read.table(filename, #select=1:3,
    #    col.names=c('l', 'n', 'nu', 'inertia'))
    ##freqs <- fread(filename, select=1:3,
    ##     col.names=c('l', 'n', 'nu'))
    seps <- seismology(freqs, nu_max)

    if ("Dnu0" %in% names(seps) & "dnu02" %in% names(seps))
        data.frame(age=prof$star_age/10**9, center_h1=prof$center_h1,
                   mass=prof$star_mass, Dnu=seps$Dnu0, dnu=seps$dnu02,
                   Teff=prof$Teff, L=prof$photosphere_L)
    else data.frame()
}

parse_dir <- function(fgong.dir) {
    all.files <- list.files(fgong.dir)
    fgong.files <- all.files[grepl('*.dat$', all.files)]

    DF <- do.call(plyr:::rbind.fill,
                  Map(separations, file.path(fgong.dir, fgong.files)))
    DF <- DF[order(DF[,1]),]
    #DF <- DF[DF$age > 0.01,]
    #DF <- DF[DF$center_h1 >= 0.001,]
    #DF <- DF[DF$dnu > -1,]

    # clip PMS
    decreasing_L <- which(diff(DF$L) < 0 & DF$center_h1[-1] > 0.6)
    if (any(decreasing_L)) {
        goes_back_up <- diff(decreasing_L) > 1
        pms <- max(decreasing_L)
        print(paste(fgong.dir, "Clipping", pms, "points"))
        DF <- DF[-1:-pms,]
    }

    #while (any(diff(DF$dnu) > 1) | any(diff(DF$dnu) < -10)) {
    #    DF <- DF[c(1, 1+which(diff(DF$dnu) <= 1)),]
    #    DF <- DF[c(1, 1+which(diff(DF$dnu) >= -10)),]
    #}

    #while (any(diff(DF$dnu) > 0)) DF <- DF[c(1, 1+which(diff(DF$dnu) <= 0)),]

    #print(DF)

    DF[DF$center_h1 > 0.01 & DF$dnu > -1,]
}



###############################################
### Plot the JCD diagram for just one track ###
###############################################
fgong.dir <- file.path('grid', 'mass', '1.0', 'LOGS')
DF <- parse_dir(fgong.dir)
cairo_pdf('CD_diagram.pdf', family='Palatino')
par(mar=c(4, 5, 1, 2))
plot(DF$Dnu, DF$dnu, type='l', lty=2, axes=F,
    xlab=expression("Large frequency separation"~"<"*Delta*nu[0]*">"/mu*Hz),
    ylab=expression("Small separation"~
        "<"*delta*nu[0*","*2]*">"/mu*Hz))
magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1, family='Palatino')
dev.off()



######################################################
### Plot the JCD diagram for a grid varied in mass ###
######################################################
h1s <- c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
masses <- c(0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0,
            #1.3, 1.4, 1.6,
            3.0, 5.0)#, 10.0)

parallelStartMulticore(min(length(masses), 16))

fgong.dirs <- file.path('grid', 'mass', sprintf('%.1f', masses), 'LOGS')
DF <- do.call(plyr:::rbind.fill, parallelMap(parse_dir, fgong.dir=fgong.dirs))
        #list.dirs(file.path('grid', 'mass'), recursive=F), 'LOGS')))
DF <- DF[order(DF$mass),]

save(DF, file='DF7')

# A helper function to obtain the interpolated large and small
# separations at a given hydrogen abundance
Dnu_dnu <- function(track, h1s) {
   data.frame(Dnu=splinefun(track$center_h1, track$Dnu, method="monoH.FC")(h1s),
              dnu=splinefun(track$center_h1, track$dnu, method="monoH.FC")(h1s))
}

# # Obtain all of the Dnus and dnus at each central H1 value
# # for which we want to plot an iso"chrone"
Dnus <- list()
dnus <- list()
for (h1 in h1s) {
    tmp_Dnus <- c()
    tmp_dnus <- c()
    for (mass in masses) {
        track <- DF[DF$mass==paste(mass),]
        vals <- Dnu_dnu(track, h1)
        tmp_Dnus <- c(tmp_Dnus, vals$Dnu)
        tmp_dnus <- c(tmp_dnus, vals$dnu)
    }
    has_val <- !is.na(tmp_Dnus) & !is.na(tmp_dnus)
    Dnus[[paste(h1)]] <- tmp_Dnus[has_val]
    dnus[[paste(h1)]] <- tmp_dnus[has_val]
}

exclude <- c(5774694) # fake Sun
pred.dir <- file.path('..', '..', 'inverse', 'learn',
                      'covs-simulations', 'legacy')
for (pred.fname in list.files(pred.dir)) {
    pred.data <- read.table(file.path(pred.dir, pred.fname), header=1)
    if (any(pred.data$X_c < 0.01))
        exclude <- c(exclude, strsplit(pred.fname, '.dat')[[1]])
}
print(exclude)


legacy.dir <- file.path('..', '..', 'regression', 'perturb', 'legacy')
legacy.fnames <- setdiff(list.files(legacy.dir),
                         paste0(exclude, '_perturb.dat'))

FeHs <- c()
for (legacy.fname in legacy.fnames) {
    legacy.data <- read.table(file.path(legacy.dir, legacy.fname), header=1,
                              nrows=1)
    FeHs <- c(FeHs, range(legacy.data$Fe))
}
print(range(FeHs))

plot_CD_diagram <- function(DF,
        col.pal=colorRampPalette(
            #c("#FDAE61", "#D7191C", "#000000", "#2B83BA", "#ABDDA4"))(15),
            c("#FC89AC", red, "#000000", blue, "#ADD8E6"))(15),
        #brewer.pal(11, 'Spectral'))(15),
        #c(blue, "black", red))(15),
        ..., text.cex=1, mgp=utils.mgp, font="Times", mar=utils.mar) {
    par(mar=mar+c(0,1,0,4), mgp=mgp+c(0, 0.2, 0))
    plot(NA, type='l', axes=F,
        xlim=range(DF$Dnu)+c(0, 30), #33), 
        ylim=c(0.5, max(DF$dnu)),#range(DF$dnu), #+c(0, 0.5),
        xlab=expression("Large frequency separation"~"<"*Delta*nu[0]*">"/mu*Hz),
        ylab=expression("Small separation"~
            "<"*delta*nu[0*","*2]*">"/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,0,0,0), las=1, family=font,
            cex.axis=text.cex, mgp=mgp+c(0, 0.25, 0))
    magaxis(side=2, tcl=-0.25, labels=T, las=1, family=font,
            cex.axis=text.cex, mgp=mgp+c(0, 0.25, 0))
    #legend("topleft", bty='n', cex=0.66, lty=1, pch=20, col="blue",
    #       legend=c("LEGACY Data"))
    #       #lty=c(1,2,1), pch=c(20,NA,NA), col=c("blue", "black", "black"),
    #       #legend=expression('Legacy Data', M, X['c']))

    # Add colorbar
    var1range <- diff(par()$usr)[1]
    color.legend(par()$usr[2],#+0.05*var1range,
                 par()$usr[3],
                 par()$usr[2]+0.05*var1range,
                 par()$usr[4], #+0.10*var1range, par()$usr[4],
                 pretty(range(FeHs)),
                 cex=text.cex,
                 col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(expression("[Fe/H]"), 4, cex=text.cex, line=4)#2)

    for (h1_i in 1:length(h1s)) {
        h1 <- h1s[h1_i]
        Dnu <- Dnus[[paste(h1)]]
        dnu <- dnus[[paste(h1)]]
        if (length(Dnu) > 3) {
            lines(dnu ~ Dnu, lty=2, lwd=2, col='#999999')
        }
    }

    for (mass in masses) {
        #with(DF[DF$mass==paste(mass),], lines(Dnu, dnu))
        track <- Dnu_dnu(DF[DF$mass==paste(mass),], h1s)
        with(track, lines(Dnu, dnu, lwd=2, col='#999999'))
        if (mass==5)
            text(track$Dnu[1]-15, #-9, 
                 track$dnu[1]-0.15, #-0.5, 
                 pos=3, labels=expression(M==5), cex=text.cex)
        else
            text(track$Dnu[1]-0.5, #+0,
                 track$dnu[1]-0.14, #-0.9, 
                 pos=3, labels=mass, cex=text.cex)
    }

    for (h1 in h1s) {
        vals <- Dnu_dnu(DF[DF$mass=="0.7",], h1)
        if (h1 == 0.7) {
            text(vals$Dnu-3, 
                 vals$dnu-0.2, 
                 labels="ZAMS", pos=4, cex=text.cex)
        } else if (h1 == 0.1) {
            text(vals$Dnu-1.5, 
                 vals$dnu, 
                 labels=expression(X["c"]==0.1), pos=4, cex=text.cex)
        } else {
            text(vals$Dnu-2, 
                 vals$dnu-0.05, 
                 labels=h1, pos=4, cex=text.cex)
        }
    }

    points(137, 9, pch=1, cex=1, lwd=2)
    points(137, 9, pch=20, cex=0.1,
           col=col.pal[floor((-min(FeHs)) / (max(FeHs)-min(FeHs))
                             * (length(col.pal)-1)) + 1])
    
    for (legacy.fname in legacy.fnames) {
        legacy.data <- read.table(file.path(legacy.dir, legacy.fname), header=1)
        Dnu <- fivenum(legacy.data[['Dnu0']])
        dnu <- fivenum(legacy.data[['dnu02']])
        FeH <- legacy.data$Fe[1]
        color <- col.pal[floor((FeH-min(FeHs)) / (max(FeHs)-min(FeHs))
                               * (length(col.pal)-1)) + 1]
        arrows(Dnu[2], dnu[3], Dnu[4], dnu[3],
            length=0.01, lwd=2, angle=90, code=3,
            col=color#'blue'
            )
        arrows(Dnu[3], dnu[2], Dnu[3], dnu[4],
            length=0.01, lwd=2, angle=90, code=3,
            col=color#'blue'
            )
        points(Dnu[3], dnu[3], pch=20, cex=text.cex, lwd=2, 
               col=color#'blue'
               )
    }
}
make_plots(plot_CD_diagram, 'CD', DF=DF, #make_png=F, short=F, thin=F,
           #paper=F, 
           #tall=F, wide=F,
           mar=utils.mar+c(0,-1,0,3.5))#2.5))



legacy.Teffs <- c()
legacy.Dnus <- c()
for (legacy.fname in legacy.fnames) {
    legacy.data <- read.table(file.path(legacy.dir, legacy.fname), header=1)
    legacy.Teffs <- c(legacy.Teffs, fivenum(legacy.data$Teff)[c(2,4)])
    legacy.Dnus <- c(legacy.Dnus, fivenum(legacy.data$Dnu0)[c(2,4)])
}

plot_Dnu_Teff <- function(DF, ..., text.cex=1, mgp=utils.mgp, font="Times") {
    plot(NA, type='l', axes=F,
         xlim=rev(range(c(legacy.Teffs, 5000, 6800))),
         ylim=rev(range(legacy.Dnus)+c(0,5)),
         xlab=expression("Effective Temperature"~T["eff"]/K),
         ylab=expression("Large separation"~"<"*Delta*nu[0]*">"/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,0,0,0), las=1, family=font,
            cex.axis=text.cex, mgp=mgp)
    magaxis(side=2, tcl=-0.25, labels=T, las=1, family=font,
            cex.axis=text.cex, mgp=mgp+c(0, 0.25, 0))
    legend("bottomleft", bty='n', cex=text.cex,
           lty=c(1,1), pch=c(20,NA), col=c("blue", "black"),
           legend=expression('Legacy data', 'Evolutionary tracks'))

    points(5777, 135, pch=1, cex=1)
    points(5777, 135, pch=20, cex=0.1)

    for (mass in unique(DF$mass)) {
        track <- DF[DF$mass==paste(mass),]
        lines(track$Dnu ~ track$Teff)
        text(DF[DF$mass==paste(mass),]$Teff[1],
             DF[DF$mass==paste(mass),]$Dnu[1]+6, cex=0.6,
             labels=if (mass==1.0) "M=1" else mass)
    }

    for (legacy.fname in legacy.fnames) {
        legacy.data <- read.table(file.path(legacy.dir, legacy.fname), header=1)
        xs <- fivenum(legacy.data[['Teff']])
        ys <- fivenum(legacy.data[['Dnu0']])
        arrows(xs[2], ys[3], xs[4], ys[3],
               length=0.01, lwd=2, angle=90, code=3, col='blue')
        arrows(xs[3], ys[2], xs[3], ys[4],
               length=0.01, lwd=2, angle=90, code=3, col='blue')
        points(xs[3], ys[3], col='blue', pch=20, cex=0.66)
    }
}
make_plots(plot_Dnu_Teff, 'Dnu_Teff', DF=DF, make_png=F, tall=F, wide=F)

#dev.off()
