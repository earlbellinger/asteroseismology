library(magicaxis)
library(lomb)
library(RColorBrewer)
source(file.path('..', 'scripts', 'utils.R'))
source(file.path('..', 'scripts', 'seismology.R'))

plot_dir <- file.path("plots", "16CygB")

freqs <- read.table(file.path('..', 'inverse', 'data', '16CygB-freqs.dat'), 
    header=1)

filenames <- c('kplr100002742-2011334093404_slc_wg1.dat',
               'kplr100002742-2011303113607_slc_wg1.dat',
               'kplr100002742-2012004120508_slc_wg1.dat')
time.series <- do.call(rbind, Map(function(filename) {
    data <- read.table(file.path('data', '16CygB', filename))
    data <- data[,c(1,4)]
    colnames(data) <- c("time", "flux")
    data <- data[is.finite(data$flux),]
    data$flux <- (data$flux - median(data$flux))/mad(data$flux)
    data
}, filename=filenames))
#time.series <- time.series[!time.series$flux%in%boxplot.stats(time.series$flux,
#                                                              coef=10)$out,]

time.series <- rbind(read.table('kplr100002742-2011334093404_slc_wg1.dat'),
                     read.table('kplr100002742-2011303113607_slc_wg1.dat'),
                     read.table('kplr100002742-2012004120508_slc_wg1.dat'))
time.series <- time.series[,c(1,4)]
colnames(time.series) <- c("time", "flux")
time.series <- time.series[is.finite(time.series$flux),]


attach(time.series)
flux <- (flux - mean(flux))/sqrt(var(flux))

#par(mfrow=c(1,2))
plot_time_series <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot((time-min(time))*24, time.series$flux, pch=4, cex=1, 
         axes=FALSE, xaxs='i', yaxs='i', tcl=0, 
         xlab=expression("Time"~tau/hours), ylab="Flux",
         xlim=c(0, 24), ylim=c(-3, 3))
    magaxis(c(2,4), labels=c(1,0), prettybase=1.5, cex=text.cex, 
            mgp=mgp, family=font, las=1)
    magaxis(c(1,3), labels=c(1,0), prettybase=4, cex=text.cex, 
            mgp=mgp, family=font, las=1)
}
make_plots(plot_time_series, "time_series", filepath=plot_dir)
make_plots(plot_time_series, "time_series-half", filepath=plot_dir, 
    slides_pdf_height=4.1511/2)

#plot(runif(1469, 0, 24), rnorm(1469, 0, 1), pch=4, cex=1, 
#     axes=FALSE, xaxs='i', yaxs='i', tcl=0, 
#     xlab="Time (hours)", ylab="Standardized Flux",
#     xlim=c(0, 24), ylim=c(-2, 2))
#magaxis(1:4, labels=c(1,1,0,0), prettybase=2)

pgram <- lsp(time.series$flux,#flux, 
             times=time.series$time*(60*60*24)/10**6, 
             from=1500, to=3500)
attach(pgram)

plot_raw_spectrum <- function(..., text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot(scanned, power/max(power),
         col="#dbdbdb", lwd=0.5,
         ylim=c(0, 1.05),
         ylab="Power", 
         xlab=expression("Frequency"~nu/mu*Hz),
         type='h', xaxs='i', yaxs='i', axes=FALSE)
    magaxis(2:4, labels=c(1,0,0), tcl=0.25, cex=text.cex, mgp=mgp, 
        family=font, prettybase=0.5)
    magaxis(1, labels=1, tcl=-0.25, cex=text.cex, mgp=mgp, family=font)
}

diff.scan <- diff(scanned)[1]
nu_cols <- rep(0, length(scanned))
for (freq_ii in 1:nrow(freqs)) {
    freq <- freqs[freq_ii,]
    unc <- min(3*freq$dnu, 2)
    colorize <- which(scanned >= freq$nu - unc &
                      scanned <= freq$nu + unc)
    colorize <- colorize[which.max(power[colorize])]
    print(length(colorize))
    if (length(colorize) == 0) print("Error!")
    nu_cols[colorize] <- freq$l+1
}

#col.pal <- brewer.pal(4, "Set1")
col.pal <- c("#ca0020", "#f4a582", "#0571b0", "#800080")

plot_power_spectrum <- function(colors=0, text.cex=1, mgp=utils.mgp, 
                                font=utils.font) {
    #plot_raw_spectrum(text.cex=text.cex, mgp=mgp, font=font)
    #points(scanned[nu_cols>0], 
    
    plot(scanned[nu_cols>0], power[nu_cols>0]/max(power[nu_cols>0]), 
         ylim=c(0, 1.05),
         ylab="Power", 
         xlim=c(1501, 3499),
         xlab=expression("Frequency"~nu/mu*Hz),
         type='h', lwd=ifelse(colors, 1.5, 1), 
         xaxs='i', yaxs='i', axes=FALSE, 
         col=if (colors) col.pal[nu_cols[nu_cols>0]] else 1)
    magaxis(2:4, labels=c(1,0,0), tcl=0.25, cex=text.cex, mgp=mgp, 
        family=font, las=1)
    magaxis(1, labels=1, tcl=-0.25, cex=text.cex, mgp=mgp, family=font)
    if (colors) {
        legend("topleft", bty="n", col=col.pal, lty=1, legend=c(
            expression("ℓ"==0),
            expression("ℓ"==1),
            expression("ℓ"==2),
            expression("ℓ"==3)
        ))
    }
}

plot_mode <- function(ells=0, annotate=FALSE,
        text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot_power_spectrum(text.cex=text.cex, mgp=mgp, font=font)
    for (ell in ells) {
        points(scanned[nu_cols==(ell+1)], 
               power[nu_cols==(ell+1)]/max(power[nu_cols>0]), 
               type='h', col=col.pal[ell+1], lwd=2)
    }
    if (annotate) {
        left <- scanned[nu_cols==(ell+1)][7]
        right <- scanned[nu_cols==(ell+1)][8]
        pwr <- power[nu_cols==(ell+1)][7]/max(power[nu_cols>0])
        arrows(left+1, 0.8*pwr, right-1, length=0.1)
        arrows(right-1, 0.8*pwr, left+1, length=0.1)
        text((left+right)/2, pwr, expression(Delta*nu))
    }
    magaxis(2:4, labels=FALSE, tcl=0.25, cex=text.cex, mgp=mgp, family=font)
    magaxis(1, labels=FALSE, tcl=-0.25, cex=text.cex, mgp=mgp, family=font)
    legend("topleft", bty="n", col=col.pal[ells+1], lty=1, 
        legend=as.expression(Map(function(ell) bquote("ℓ"==.(ell)), ell=ells)))
}

plot_nu_max <- function(text.cex=1, mgp=utils.mgp, font=utils.font) {
    plot_power_spectrum(text.cex=text.cex, mgp=mgp, font=font)
    x <- seq(1500,3500,length=10000)
    y <- dnorm(x, mean=2552, sd=(0.66*2552**0.88)/(2*sqrt(2*log(2))))
    lines(x, 1.025*y/max(y), lty=2)
}

make_plots(plot_power_spectrum, 'power_spectrum', filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l0_mode', filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l0_mode-Dnu', annotate=TRUE, filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l1_mode', ells=1, filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l2_mode', ells=2, filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l3_mode', ells=3, filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_mode, 'l02_modes', ells=c(0,2), filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_power_spectrum, 'all_modes', colors=1, filepath=plot_dir,
    slides_pdf_height=4.1511/2)
make_plots(plot_nu_max, 'nu_max', filepath=plot_dir,
    slides_pdf_height=4.1511/2)

