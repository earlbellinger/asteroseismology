#### Animate asteroseismological quantities across an evolutionary track
#### Author: Earl Bellinger ( bellinger@mps.mpg.de )
#### Department of Astronomy, Yale University;
#### Stellar Ages & Galactic Evolution Group,
#### Max-Planck-Institut fur Sonnensystemforschung

library(parallel)
library(parallelMap)
library(data.table)
library(pracma)
library(akima)
library(Bolstad)
library(RColorBrewer)
source(file.path('/', 'scratch', 'seismo', 'bellinger',
    'asteroseismology', 'scripts', 'seismology.R'))

setwd('gemma')#nosmooth')#

parallelStartMulticore(16)#max(1,as.integer(Sys.getenv()[['OMP_NUM_THREADS']])))

transblack <- adjustcolor('black', alpha.f=0.66)

ells <- c(10, 11, 9, 3)

log_dir <- 'LOGS'
ev.DF <- read.table(file.path(log_dir, 'history.data'), header=TRUE, skip=5)
logs <- list.files(log_dir)

freq_files <- logs[grep('profile.+-freqs.dat$', logs)]
prof_files <- sub('-freqs.dat', '.data', freq_files)

mdl_nums <- as.numeric(parallelMap(function(pro_file)
        read.table(pro_file, header=TRUE, nrows=1, skip=1)$model_number,
    pro_file=file.path(log_dir, prof_files)))
ages <- as.numeric(parallelMap(function(pro_file)
        read.table(pro_file, header=TRUE, nrows=1, skip=1)$star_age,
    pro_file=file.path(log_dir, prof_files))) / 10**9
min.age <- min(ages)
ages <- ages[order(mdl_nums)] - min.age
mdl_nums. <- sort(mdl_nums)

burn.max <- ceil(max(as.numeric(parallelMap(function(pro_file)
        max(read.table(pro_file, header=TRUE, skip=5)$eps_nuc),
    pro_file=file.path(log_dir, prof_files)))))

mdl.selector <- T

get_freqs <- function(model_number) {
    idx <- which(mdl_nums == model_number)
    freq_file <- freq_files[idx]
    parse_freqs(file.path(log_dir, freq_file))#, gyre=T)
}

get_closest <- function(freqs, nu_max, ells) {
    new.freqs <- data.frame()
    for (ell in 1:length(ells)) {
        freqs.ell <- freqs[freqs$l==ell-1,]
        for (ii in 1:ells[ell]) {
            closest <- find_closest2(freqs.ell$nu, nu_max)
            new.freqs <- plyr:::rbind.fill(new.freqs, freqs.ell[closest,])
            freqs.ell <- freqs.ell[-closest,]
        }
    }
    new.freqs[order(new.freqs$l, new.freqs$n),]
}

collapse_freqs <- function(model_number) {
    print(model_number)
    DF <- get_freqs(model_number)
    if (nrow(DF) <= 0) return(data.frame(NA))
    #DF <- DF[DF$l == 0,]
    nu_max <- ev.DF[ev.DF$model_number==model_number,]$nu_max
    DF <- get_closest(freqs=DF, nu_max=nu_max, ells=ells)
    freqs <- rbind(DF$nu)
    colnames(freqs) <- paste0('l', DF$l, '.', 'n', DF$n)
    as.data.frame(freqs)
}

track <- ev.DF#[ev.DF$model_number %in% mdl_nums[mdl.selector],]

track_freqs <- do.call(plyr:::rbind.fill, 
    parallelMap(collapse_freqs, sort(mdl_nums[mdl.selector])))



k.pair <- 'u-Y' #'psi' #'Gamma1-rho'#'c2-rho'#
mode.i <- 'l.1_n.11'

for (k.pair in c('psi', 'Gamma1-rho', 'c2-rho', 'u-Y')) {
for (mode.i in c('l.1_n.11', 'l.0_n.11', 'l.1_n.15')) {

out_dir <- paste0('animate-', k.pair, '-', mode.i)
if (k.pair == 'psi') {
    k.fname1 <- 'psi.dat'
} else {
    k.f1 <- strsplit(k.pair, '-')[[1]][1]
    k.f2 <- strsplit(k.pair, '-')[[1]][2]
    k.fname1 <- paste0('E_K_', k.pair, '.dat') #'E_K_u-Y.dat'
    k.fname2 <- paste0('E_K_', k.f2, '-', k.f1, '.dat')
}

add_convection <- function(DF, ylim, conv.col="#D9D6BE",
        conv.core.col="#FFE1A8") {
    yrange <- c(ylim[1]-0.5*abs(ylim[1]), ylim[2]+0.5*abs(ylim[2]))
    radius <- DF$radius / max(DF$radius)
    bounds <- which(diff(DF$gradr - DF$grada < 0) != 0)
    conv <- DF$gradr > DF$grada
    conv.core <- conv[length(conv)]
    if (conv.core) {
        rect(10**-3, yrange[1], radius[bounds[length(bounds)]], yrange[2],
            lwd=2, col=conv.core.col, border=transblack, lty=3)
        bounds <- bounds[-length(bounds)]
    }
    bounds <- rev(bounds)
    while (length(bounds)>1) {
        rect(radius[bounds[1]], yrange[1], 
             radius[bounds[2]], yrange[2],
            lwd=2, col=conv.col, border=transblack, lty=3)
        bounds <- bounds[-c(1,2)]
    }
    if (length(bounds) == 1) {
       rect(radius[bounds[1]], yrange[1], 1, yrange[2],
            lwd=2, col=conv.col, border=transblack, lty=3)
    }
}

add_burning <- function(DF, ylim) {
    yrange <- c(ylim[1]-0.5*abs(ylim[1]), ylim[2]+0.5*abs(ylim[2]))
    radius <- DF$radius / max(DF$radius)
    
    rs <- seq(0.001, 1, 0.001)
    #burn.max <- 200 
    #burning <- floor(splinefun(radius, DF$eps_nuc)(rs)/burn.max)+1
    eps.nuc <- splinefun(radius, DF$eps_nuc)(rs)
    eps.nuc[eps.nuc<=0] <- 10**-100
    num.cols <- 1001
    burning <- ceil(num.cols*log10(eps.nuc)/log10(burn.max))
    burning[burning<1] <- 1
    #log.burn.max <- 100*log10(burn.max) / log10(burn.max) + 1
    cols <- sapply(1:num.cols, function(burn)
        adjustcolor(colorRampPalette(c(red, '#780000'))(num.cols)[burn], 
            alpha.f=ifelse(burn>1, 1.0025**burn / 1.0025**num.cols, 0)))
            #log10(burn/num.cols)/log10(num.cols)))#exp(burn/num.cols)/exp(1)))
    #burning[burning>burn.max] <- burn.max
    #cols <- sapply(1:burn.max-1, function(alpha.f)
    #    adjustcolor(colorRampPalette(c(red, '#780000'))(burn.max)[alpha.f+1], 
    #        alpha.f=alpha.f/burn.max))
    #
    #    adjustcolor(
    #    colorRampPalette(c('white', red, 'darkred'))(25)[1:25],#[burning],
    #    alpha.f=0.75)
    
    burn <- 1
    start.r <- rs[1]
    #burn <- F
    #start.r <- -1
    for (ii in 1:length(burning)) { 
        if (burning[ii] == burn) next 
        end.r <- rs[ii] #if (ii < length(burning)) rs[ii+1] else rs[ii]
        if (burn>1) rect(start.r, yrange[1], end.r, yrange[2], 
            col=cols[burn], border=NA)
        burn <- burning[ii]
        start.r <- rs[ii]
    }
    if (F) {
        if (!burn) { #& burning[ii]>1) {
            burn <- burning[ii]
            start.r <- rs[ii]
        } else { #if (burn) {# & burning[ii]>1) {
            if (burning[ii] == burn) next
            end.r <- if (ii < length(burning)) rs[ii+1] else rs[ii]
            if (burn>1)
                rect(start.r, yrange[1], end.r, yrange[2], 
                    col=cols[burn], border=NA)
            burn <- burning[ii]
            start.r <- rs[ii]
        }
    }
    #r.min <- rs[which.min(burning>1)] # 0.001
    #r.max <- rs[which.max(burning>1)] # 1
    #gradient.rect(r.min, yrange[1], r.max, yrange[2],
    #    col=colorRampPalette(c('white', red, 'darkred'))(25)[burning])
}

hr_echelle <- function(model_number) {
    par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(3,5,1,1), cex=1.5, pch=20)
    
    selector <- track$model_number == model_number
    model <- track[selector,]
    idx <- which(mdl_nums == model_number)
    idx. <- which(mdl_nums. == model_number)
    DF <- read.table(file.path(log_dir, prof_files[idx]), 
        skip=5, header=1)
    radius <- read.table(file.path(log_dir, prof_files[idx]), 
        skip=1, nrow=1, header=1)$photosphere_r
    freqs <- get_freqs(model_number)
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freqs <- get_closest(freqs=freqs, nu_max=nu_max, ells=ells)
    large_sep <- avg(Dnu, freqs)[[1]]
    freq_range <- range(freqs$nu) 
    age <- model$star_age/10**9 - min.age
    
    
    
    ## hr
    with(ev.DF, 
        plot(log_Teff, log_L, axes=F, tcl=0, type='l', lwd=2, 
            xlim=rev(range(log_Teff)), 
            xlab="",
            ylab=expression("Luminosity" ~ log(L/L[Sun]))))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)#,
    #magaxis(side=1:4, tcl=-0.25, labels=c(0,0,0,0), las=1)#,
        #ylab="Luminosity log(L/L\\sb\\SO\\eb)", hersh=T)
    #axis(1, at=axTicks(1), labels=signif(10**axTicks(1), 4), tick=F)
    #axis(2, at=axTicks(2), labels=signif(10**axTicks(2), 3), tick=F, las=1)
    par(mgp=c(3, 0.15, 0))
    axis(3, at=log10(c(41000, 31000, 9500, 7240, 5920, 5300, 3850)),
        labels=c("O", "B", "A", "F", "G", "K", "M"), tcl=0)
    par(mgp=c(3, 1, 0))
    with(model, points(log_Teff, log_L))
    points(log10(5777.7), 0)
    points(log10(5777.7), 0, pch=20, cex=0.1)
    legend("bottomright", bty='n', #inset=c(-0.1, -0.05),
        legend=c(paste("Age:", signif(age, 4), "/", 
                       signif(max(ages), 4), "Gyr"),
                 paste("Model:", 
                       model$model_number-min(track$model_number)+1, "/",
                       max(track$model_number)-min(track$model_number)+1)))
    par(mgp=c(2, 1, 0))
    title(xlab=expression("Temperature" ~ log(T["eff"]/K)))
    par(mgp=c(3, 1, 0))
    
    
    
    ## u kernel
    lowest.l1 <- freqs[freqs$l==1&freqs$n==11,]#freqs[freqs$l==1,][1,]
    #k.fname1 <- paste0('E_K_', k.pair, '.dat') #'E_K_u-Y.dat'
    #k.fname2 <- paste0('E_K_', 
    #    strsplit(k.pair, '-')[[1]][2], '-', strsplit(k.pair, '-')[[1]][1], 
    #    '.dat')
    u.Y <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]), k.fname1), header=1)
    if (k.pair == 'psi') {
        ylab <- bquote(psi)
        Y.u <- u.Y$x * 0
    } else {
        ylab <- expression("Kernel" ~ K)#^(u*","~Y))
        Y.u <- read.table(file.path('LOGS', 
              sub('.data', '-freqs', prof_files[idx]), k.fname2), 
          header=1)[[mode.i]]
    }
    K <- u.Y[[mode.i]]
    x <- u.Y$x
    x[1] <- 10**-10
    if (is.null(K)) K <- x
    inside <- K[x<0.98]
    ylim <- range(abs(inside), -abs(inside))
    plot(NA, axes=F, tcl=0, log='x',
        xaxs='i',
        xlab="",
        ylab=ylab,
        xlim=c(10**-3, 1),#c(0, 1),#
        ylim=ylim)
    add_convection(DF, ylim)
    #abline(h=0, col='darkgray', lwd=1.5, lty=2)
    #lines(x[-1], K[-1], lwd=2, col='darkred')#col="#11151C")#
    #lines(x[-1], Y.u[[mode.i]][-1], lwd=2, lty=2, col='black') #"#31572C") 
    lines(x, K, lwd=2, col='darkred')#col="#11151C")#
    lines(x, Y.u, lwd=2, lty=2, col='black') #"#31572C") 
    legend("bottom", bty='n', 
        legend=gsub('_', ', ', gsub('\\.', '=', mode.i)))
        #bquote(l==1*","~n==11))
        #bquote(l==.(lowest.l1$l)*","~n==.(lowest.l1$n)))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    par(mgp=c(2, 1, 0))
    title(xlab=expression("Radius" ~ r/R["star"]))
    par(mgp=c(3, 1, 0))
    
    
    
    ## freq_evolv_plot
    l.0.col <- "#2B2D42" #'black'
    l.1.col <- "#BF211E" #red
    hilight <- "#E09F3E" #"#FFF45F"
    radials <- track_freqs[grepl('l0', colnames(track_freqs))]
    dipoles <- track_freqs[grepl('l1', colnames(track_freqs))]
    age.range <- which(ages>(age-0.15) & ages<(age+0.15))
    
    plot(NA, axes=F, xaxs='i',
        xlim=c(max(min(ages), age-.15),
               min(max(ages), age+.15)), 
        #ylim=range(radials[selector,], na.rm=1), 
        ylim=range(radials[age.range,], dipoles[age.range,], freq_range, 
            na.rm=T),
        xlab="",
        ylab=expression("Frequency"~nu/mu*Hz))
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    lines(ages, track_freqs[['l1.n11']], lwd=4, col=hilight, lty=1)
    for (col_i in 1:ncol(dipoles)) 
        lines(ages, dipoles[,col_i], lty=1, lwd=2, col=l.1.col)
    for (col_i in 1:ncol(radials)) 
        lines(ages, radials[,col_i], lty=2, lwd=2, col=l.0.col)
    #lines(ages, track_freqs[['l1.n11']], lwd=1, col="#FFF45F", lty=1)
    points(age, track_freqs[['l1.n11']][idx.], cex=1.5, col=hilight)
    for (row_i in 1:ncol(radials)) 
        points(age, radials[idx., row_i], cex=0.5, col=l.0.col)
    for (row_i in 1:ncol(dipoles)) 
        points(age, dipoles[idx., row_i], cex=0.5, col=l.1.col)
    points(age, track_freqs[['l1.n11']][idx.], cex=1.5, col=hilight)
    points(age, track_freqs[['l1.n11']][idx.], cex=0.5, col=l.1.col)
    legend("topleft", lty=c(2,1), pch=20, cex=0.8, 
           col=c(l.0.col, l.1.col), 
           inset=c(0.02, 0.02), legend=c(paste0("\u2113=", 0:1)))
    par(mgp=c(2, 1, 0))
    title(xlab=expression("Star age"~tau/Gyr))
    par(mgp=c(3, 1, 0))
    
    
    
    ## propagation
    pprop.col <- adjustcolor("#F9E27F", alpha.f=0.75) #"#F8F4A6"
    lambs.col <- "#3C6E71" #"#70AE6E" #"#E08E45" #blue
    brunt.col <- "black" #"#3C6E71" # 'black'
    gprop.col <- adjustcolor("#E2F0FB", alpha.f=0.75)
    
    ell <- 1
    S <- DF$csound**2 / (DF$radius * solar_radius)**2
    lambs <- 10**6/(2*pi) * sqrt(ell*(ell+1) * S)
    
    brunt_N2 <- DF$brunt_N2
    brunt <- 10**6/(2*pi) * sqrt(brunt_N2)
    brunt[is.na(brunt)] <- 1
    
    ylim <- range(10, 10**5)
    plot(NA, axes=F, tcl=0, type='l', 
        log='xy', xaxs='i', #cex=0.5,
        #ylim=freq_range,
        ylim=ylim,
             #c(max(1, nu_max-20*large_sep), 
             #  min(10**5, nu_max+20*large_sep)),
        #range(1, 10**5),#range(lambs),#, 0.01, 10**6), 
        #ylim=c(1, 10**5),
        #xlim=range(DF$radius/max(DF$radius)),#range(0.01, max(DF$radius)), 
        xlim=c(10**-3, 1),#c(0, 1),#
        xlab="",
        ylab=expression("Frequency" ~ nu/mu*Hz))
    
    radius <- DF$radius/max(DF$radius)
    yrange <- c(ylim[1]-0.5*abs(ylim[1]), ylim[2]+0.5*abs(ylim[2]))
    #nus <- 10**seq(log10(yrange[1]), log10(yrange[2]), length=1000)
    g.bot <- sapply(1:length(radius), function(ii) 1)
    g.top <- sapply(1:length(radius), function(ii) min(brunt[ii], lambs[ii]))
    p.bot <- sapply(1:length(radius), function(ii) max(brunt[ii], lambs[ii]))
    p.top <- sapply(1:length(radius), function(ii) yrange[2])
    
    g.x <- c(radius[g.top>1], rev(radius[g.top>1]))
    g.y <- c(g.top[g.top>1], g.bot[g.top>1])
    polygon(x = g.x, y = g.y, 
        col=gprop.col, border=NA)
    polygon(x = c(radius, rev(radius)), y = c(p.bot, p.top), 
        col=pprop.col, border=NA)
    
    add_burning(DF, ylim)
    
    abline(h=freq_range[1], lty=3, lwd=2, col=transblack)
    abline(h=freq_range[2], lty=3, lwd=2, col=transblack)
    magaxis(side=1:4, tcl=-0.25, labels=c(1,1,0,0), las=1)
    
    #for (ell in 1:3) lines(DF$radius/max(DF$radius), lambs[[ell]], lty=ell+1)
    bounds <- which(diff(DF$gradr - DF$grada < 0) != 0)
    abline(v=DF$radius[bounds]/max(DF$radius), lty=3, lwd=2, col=transblack)
    lines(DF$radius/max(DF$radius), lambs, lty=2, lwd=2, col=lambs.col)
    lines(DF$radius/max(DF$radius), brunt, lwd=2, col=brunt.col)
    #legend("top", lty=c(1,2), lwd=2, col=c('black', blue), #pch=c(1, NA),
    #    legend=c("Brunt", "Lamb"), cex=0.8, inset=c(0.05, 0))
    #legend("topright", cex=0.8, 
    par(mgp=c(2, 1, 0))
    title(xlab=expression("Radius" ~ r/R["star"]))
    par(mgp=c(3, 1, 0))
}

model_number <- mdl_nums.[1]
options(bitmapType='cairo')
png('Rplots.png', width=1600, height=900, res=140)
hr_echelle(model_number)
dev.off()

dir.create(file.path(out_dir), showWarnings = FALSE)
options(bitmapType='cairo')
#for (model_number in sort(mdl_nums[mdl.selector])) {
#    png(file.path('animate', 
#            paste0(formatC(model_number, width=6, format='d', flag=0), '.png')),
#        width=800, height=600)
#    hr_echelle(model_number)
#    dev.off()
#}
parallelMap(function(model_number) {
    png(file.path(out_dir, 
            paste0(formatC(model_number, width=6, format='d', flag=0), '.png')),
        width=1600, height=900, res=140)
    hr_echelle(model_number)
    dev.off()
}, model_number=sort(mdl_nums[mdl.selector]))







for (k.fname in c(k.fname1, k.fname2)) {

if (k.fname == "E_K_Gamma1-rho.dat") {
    plot.fname <- 'Gamma1-rho'
    k1 <- bquote(Gamma[1])
    k2 <- bquote(rho)
} else if (k.fname == "E_K_rho-Gamma1.dat") {
    plot.fname <- 'rho-Gamma1'
    k1 <- bquote(rho)
    k2 <- bquote(Gamma[1])
} else if (k.fname == "E_K_u-Y.dat") {
    plot.fname <- 'u-Y'
    k1 <- bquote(u)
    k2 <- bquote(Y)
} else if (k.fname == "E_K_Y-u.dat") {
    plot.fname <- 'Y-u'
    k1 <- bquote(Y)
    k2 <- bquote(u)
} else if (k.fname == "E_K_c2-rho.dat") {
    plot.fname <- 'c2-rho'
    k1 <- bquote(c^2)
    k2 <- bquote(rho)
} else if (k.fname == "E_K_rho-c2.dat") {
    plot.fname <- 'rho-c2'
    k1 <- bquote(rho)
    k2 <- bquote(c^2)
}

#k1 <- gsub('-.+', '', gsub('.+_', '', k.fname1))
#k2 <- gsub('-.+', '', gsub('.+_', '', k.fname1))

#for (mode.i in c('l.1_n.11', 'l.1_n.15', 'l.0_n.11')) {
xs <- seq(0.001, 1, 0.001)
r.max <- do.call(rbind, parallelMap(function(model_number) {
    selector <- track$model_number == model_number
    model <- track[selector,]
    idx <- which(mdl_nums == model_number)
    freqs <- get_freqs(model_number)
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freqs <- get_closest(freqs=freqs, nu_max=nu_max, ells=ells)
    
    u.Y <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]), k.fname), header=1)
    if (ncol(u.Y)<=1) {
        print(model_number)
        return(data.frame())
    }
    
    K_i <- u.Y[[mode.i]]
    #K_i <- abs(K_i) / sintegral(u.Y$x, abs(K_i))$value
    K_i <- splinefun(u.Y$x/max(u.Y$x), K_i)(xs)
    K_i / max(abs(K_i))
}, model_number=sort(mdl_nums[mdl.selector])))
age.selector <- T#ages>=3.9
x.selector <- T
r.max. <- r.max#[,-1]

plot_contour <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(0.5, 0.5, 0, 4))
    
    filled.contour(10**(ages[age.selector]), log10(xs[x.selector]), 
        r.max.[age.selector,x.selector], 
        log='y',
        nlevels=22,
        xlim=c(0, 10**max(ages[age.selector])),
        color=colorRampPalette(c('darkblue', 'white', 'darkred')),
        xaxs='i', yaxs='i',
        key.axes={
            #par(fg='black')
            axis(4, tcl=0, line=0)
            mtext(bquote("Normalized Kernel"~K^(.(k1)*','~.(k2))), 
                side=4, las=3, line=3, family=font, cex.axis=text.cex)
        },
        plot.axes={
            abline(v=10**5.589, lty=2, lwd=1)
            
            magaxis(side=c(2,4), tcl=0.25, labels=c(1,0),
                unlog='y', family=font, cex.axis=text.cex)
            
            #yticks <- pretty(10**(ages[age.selector]))
            labs <- signif(log10(pretty(10**ages[age.selector])), 2)
            labs[1] <- 0
            yticks <- 10**labs
            #10**c(0, 4, 4.3, 4.6, max(ages[age.selector]))
            #yticks[1] <- 1
            #labs <- signif(log10(yticks), 3)
            
            axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
                labels=labs, tcl=0.25)
            axis(1, tick=T, at=10**seq(0, max(ages), 0.05), tcl=0.125, labels=F)
            
            axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
                labels=F, tcl=0.25)
            axis(3, tick=T, at=10**seq(0, max(ages), 0.05), tcl=0.125, labels=F)
        },
        plot.title={
            title(xlab=expression("Age"~tau/Gyr), line=2)
            title(ylab=expression("Radius"~r/R["*"]), line=2)
        })
    
    #par(mgp=mgp+c(1.5, 0, 0))
    #title(ylab=bquote(''))
}
make_plots(plot_contour, 
    paste0("kernel_evolution-", plot.fname, '-', mode.i), 
    short=F, thin=F)

}


}


}













#for (mode.i in c('l.1_n.11', 'l.1_n.15', 'l.0_n.11')) {
xs <- seq(0.001, 1, 0.001)
r.max <- do.call(rbind, parallelMap(function(model_number) {
    selector <- track$model_number == model_number
    model <- track[selector,]
    idx <- which(mdl_nums == model_number)
    
    dgam1 <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]), 'gamder.out'), 
      col.names=c('x', 'p', 'rho', 'y'))
    
    splinefun(dgam1[['x']], dgam1[['P']])(xs)
}, model_number=sort(mdl_nums[mdl.selector])))
age.selector <- T#ages>=3.9
x.selector <- T
r.max. <- r.max#[,-1]

plot_contour <- function(..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
                           font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(0.5, 0.5, 0, 4))
    
    filled.contour(10**(ages[age.selector]), #log10(xs[x.selector]), 
        xs[x.selector],
        r.max.[age.selector,x.selector], 
        #log='y',
        nlevels=22,
        xlim=c(0, 10**max(ages[age.selector])),
        color=colorRampPalette(c('darkblue', 'white', 'darkred')),
        xaxs='i', yaxs='i',
        key.axes={
            #par(fg='black')
            axis(4, tcl=0, line=0)
            mtext(bquote(
                bgroup('(', 
                    frac(delta~ln~Gamma[1], delta~ln~P),#rho),#Y),#
                ')')[rho*","*Y]
            ), 
                side=4, las=3, line=3, family=font, cex.axis=text.cex)
        },
        plot.axes={
            abline(v=10**5.589, lty=2, lwd=1)
            
            magaxis(side=c(2,4), tcl=0.25, labels=c(1,0),
                #unlog='y', 
                family=font, cex.axis=text.cex)
            
            labs <- signif(log10(pretty(10**ages[age.selector])), 2)
            labs[1] <- 0
            yticks <- 10**labs
            
            axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
                labels=labs, tcl=0.25)
            axis(1, tick=T, at=10**seq(0, max(ages), 0.05), tcl=0.125, labels=F)
            
            axis(1, tick=T, at=yticks, cex.axis=text.cex, las=1,
                labels=F, tcl=0.25)
            axis(3, tick=T, at=10**seq(0, max(ages), 0.05), tcl=0.125, labels=F)
        },
        plot.title={
            title(xlab=expression("Age"~tau/Gyr), line=2)
            title(ylab=expression("Radius"~r/R["*"]), line=2)
        })
    
    #par(mgp=mgp+c(1.5, 0, 0))
    #title(ylab=bquote(''))
}
make_plots(plot_contour, paste0("dgamma1_P"), short=F, thin=F)
#make_plots(plot_contour, paste0("dgamma1_rho"), short=F, thin=F)
#make_plots(plot_contour, paste0("dgamma1_Y"), short=F, thin=F)














if (F) {


#ages <- r.max[,1]
#ages <- ages - min(ages)

image(ages[age.selector], xs[x.selector], 
      r.max.[age.selector,x.selector], log='y',
      xaxs='i', yaxs='i',
    xlab=expression("Age"~tau/Gyr),
    ylab=expression("Radius"~r/R["star"]),
    col=colorRampPalette(c(blue, 'white', red))(21))
dev.off()









r.max <- do.call(plyr:::rbind.fill, parallelMap(function(model_number) {
    selector <- track$model_number == model_number
    model <- track[selector,]
    idx <- which(mdl_nums == model_number)
    freqs <- get_freqs(model_number)
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freqs <- get_closest(freqs=freqs, nu_max=nu_max, ells=ells)
    
    u.Y <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]),
      'E_K_u-Y.dat'), header=1)
    if (ncol(u.Y)<=1) {
        print(model_number)
        return(data.frame())
    }
    
    info <- do.call(cbind, parallelMap(function(freq_i) {
        freq.i <- freqs[freq_i,]
        K_i <- u.Y[[paste0('l.', freq.i$l, '_n.', freq.i$n)]]
        cumulative <- cumtrapz(u.Y$x, K_i / sintegral(u.Y$x, K_i)$value)
        #u.Y$x[which.max()]
        u.Y$x[which.min(cumulative<=0.5)]
    }, freq_i=1:nrow(freqs)))
    colnames(info) <- paste0('l.', freqs$l, '_n.', freqs$n)
    
    age <- model$star_age/10**9
    as.data.frame(cbind(age, info))
}, model_number=sort(mdl_nums[mdl.selector])))

plot_r.max <- function(r.max, ..., text.cex=1, mgp=utils.mgp, mar=utils.mar,
        font="Times") {
    text.cex <- text.cex*1.25
    par(cex.lab=text.cex, mar=mar+c(-0.5, 0.5, 0, 0))
    
    ages <- r.max[,1]
    r.max. <- r.max[,-1]
    plot(NA, axes=F, 
         xlim=range(ages),#range(3, max(ages)),
         ylim=range(r.max., na.rm=T),
         xlab=expression("Star age"~tau/"Gyr"),
         ylab="")
    
    for (ii in 1:ncol(r.max.)) {
        #ell <- if (ii <= cumsum(ells)[1]) 1 else if (ii <= cumsum(ells)[2])
        #    2 else if (ii <= cumsum(ells)[3]) 3 else 4
        ell <- if ('l.0' %in% colnames(r.max.)[ii]) 1 else if 
            ('l.1' %in% colnames(r.max.)[ii]) 2 else if
            ('l.2' %in% colnames(r.max.)[ii]) 3 else 4
        lines(ages, r.max.[,ii], lwd=2, lty=ell, col=ell)
    }
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0),
            las=1, mgp=mgp-c(0.5,0.15,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=1, 
            las=1, mgp=mgp+c(1,0,0),
            family=font, cex.axis=text.cex)
    
    #legend("topleft", lty=c(2, 1), col=c('darkgray', 'darkred'), lwd=2, 
    #       inset=c(0.01, 0.015),
    #       pch=c(NA, 20), legend=c("Actual", "Inverted"))
    
    par(mgp=mgp+c(1.5, 0, 0))
    title(ylab=bquote("Radius of"~max(K^(u*","~Y))))

}


xs <- seq(0.01, 1, 0.001)
r.max <- do.call(rbind, parallelMap(function(model_number) {
    selector <- track$model_number == model_number
    model <- track[selector,]
    idx <- which(mdl_nums == model_number)
    freqs <- get_freqs(model_number)
    nu_max <- with(model, nu_max_scaling(star_mass, 10**log_R, 10**log_Teff))
    freqs <- get_closest(freqs=freqs, nu_max=nu_max, ells=ells)
    
    u.Y <- read.table(file.path('LOGS', 
      sub('.data', '-freqs', prof_files[idx]),
      'E_K_u-Y.dat'), header=1)
    if (ncol(u.Y)<=1) {
        print(model_number)
        return(data.frame())
    }
    
    info <- do.call(cbind, parallelMap(function(freq_i) {
        freq.i <- freqs[freq_i,]
        K_i <- u.Y[[paste0('l.', freq.i$l, '_n.', freq.i$n)]]
        K_i <- abs(K_i) / sintegral(u.Y$x, abs(K_i))$value
        splinefun(u.Y$x, K_i)(xs)
    }, freq_i=1:nrow(freqs)))
    c(model$star_age/10**9, rowSums(info))
}, model_number=sort(mdl_nums[mdl.selector])))

ages <- r.max[,1]
r.max. <- r.max[,-1]

age.selector <- ages>=4
image(ages[age.selector], xs[1:900], r.max.[age.selector,1:900], log='y',
    xlab=expression("Age"~tau/Gyr),
    ylab=expression("Radius"~r/R["star"]),
    col=colorRampPalette(c('white', blue, red))(21))
dev.off()


levelplot(K~age*x, 
    data=list(age=ages[age.selector], 
              x=log10(xs[x.selector]), 
              K=r.max.[age.selector,x.selector]),
  xlab = "X Coordinate (feet)", ylab = "Y Coordinate (feet)",
  main = "Surface elevation data",
  col.regions = terrain.colors(100)
)



image.plot(ages[age.selector], xs[x.selector], 
      r.max.[age.selector,x.selector], log='y',
      nlevel=100,
    xlab=expression("Age"~tau/Gyr),
    ylab=expression("Radius"~r/R["star"]),
    col=colorRampPalette(c('darkblue', 'white', 'darkred'))(21))
dev.off()


drape.plot(ages[age.selector], xs[x.selector], 
      r.max.[age.selector,x.selector],
      theta=0, phi=45)
dev.off()


library(lattice)
g <- expand.grid(x=ages[age.selector], y=xs[x.selector])
g$z <- c(r.max.[age.selector,x.selector])
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
wireframe(z ~ x * y, g, drape = TRUE,
           perspective = FALSE,
           aspect = c(3,1), colorkey = FALSE)

        
g <- expand.grid(x=ages[age.selector], y=xs[x.selector])
g$z <- c(r.max.[age.selector,x.selector])        
mesh <- interp(g$x, g$y, g$z)
filled.contour(mesh, 
    log='y',
    nlevels=22,
    color=colorRampPalette(c('darkblue', 'white', 'darkred')),
    xaxs='i', yaxs='i',
    key.axes={
        axis(4, tcl=0, line=0)
        mtext(expression("Normalized Kernel"~K^(u*','~Y)), 
            side=4, las=3, line=3)
    },
    plot.axes={
        magaxis(side=1:4, family=utils.font, tcl=0.25, labels=c(1,1,0,0),
            unlog='y')
    },
    plot.title={
        title(xlab=expression("Age"~tau/Gyr), line=2)
        title(ylab=expression("Radius"~r/R["star"]), line=2)
    })
dev.off()

x.selector <- c(T, F, F, F, F)
persp(ages[age.selector], xs[x.selector], 
      r.max.[age.selector,x.selector], phi=45, theta=0)
dev.off()

}
