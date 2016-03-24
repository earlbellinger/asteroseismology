#### Plot random forest OOB score from learn.py 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(magicaxis)
library(RColorBrewer)
source(file.path('..', 'scripts', 'utils.R'))

data <- read.table(header=0, col.names=c("n_trees", "oob_estimate", "time"), 
    textConnection(
"4 -1.59092383389 1.0694997310638428
5 -0.699668004518 1.0913584232330322
6 -0.0770409535235 1.1505346298217773
7 0.275452343631 1.1438889503479004
8 0.524621385099 1.338620662689209
9 0.676525305658 1.4479658603668213
10 0.763739335917 1.2791359424591064
11 0.830172261431 1.233341932296753
12 0.870283188661 1.2746436595916748
13 0.899253443168 1.1950063705444336
14 0.91282219679 1.3169989585876465
15 0.92588248411 1.3896033763885498
16 0.93393885205 1.4302539825439453
18 0.942025728178 1.5908150672912598
20 0.949021334034 1.4341449737548828
24 0.954036496768 1.73848557472229
32 0.958056673017 2.1601779460906982
64 0.964926578017 3.4933571815490723
128 0.967952290197 6.176987648010254
256 0.969232403824 13.008396625518799
512 0.96991333332 25.358351945877075
1024 0.970277333503 53.30591607093811"))
#2048 0.970496444205 118.6785409450531"))

data2 <- read.table(header=0, col.names=c("n_trees", "oob_estimate", "time"), 
    textConnection(
"4 -8.21757538345 1.0032823085784912
5 -4.90218872624 0.9420266151428223
6 -2.80996270248 0.9958856105804443
7 -1.37211318706 1.0066289901733398
8 -0.48200043081 1.0346086025238037
9 0.00230861943514 1.0545179843902588
10 0.353116747938 1.0928094387054443
11 0.586527068599 1.1130430698394775
12 0.714980369202 1.1459684371948242
13 0.779206078954 1.2733943462371826
14 0.85724069986 1.409426212310791
15 0.897342452183 1.4091274738311768
16 0.908976841376 1.2863829135894775
18 0.939183440084 1.5827152729034424
20 0.950468547049 1.4342658519744873
24 0.955553546978 1.7327349185943604
32 0.960371165712 2.0115349292755127
64 0.965803430183 3.467203378677368
128 0.968348916194 6.066576719284058
256 0.969417524348 12.017822504043579
512 0.970192674268 24.709973335266113
1024 0.970414603705 44.51620101928711"))
#2048 0.970666178286 94.13048672676086"))

data3 <- read.table(header=0, col.names=c("n_trees", "oob_estimate", "time"), 
    textConnection(
"4 -7.52520545308 0.9224388599395752
5 -4.42919701112 0.9668724536895752
6 -2.41010250401 0.9904253482818604
7 -1.2427633602 1.0408079624176025
8 -0.356889368714 1.1731603145599365
9 0.0994669956148 1.0059773921966553
10 0.379441879423 1.1420676708221436
11 0.593313836189 1.1808910369873047
12 0.737979647719 1.194148063659668
13 0.807072433289 1.3135206699371338
14 0.854028833923 1.2673249244689941
15 0.901904000605 1.4008500576019287
16 0.908824743319 1.3584856986999512
18 0.93139850311 1.4456632137298584
20 0.948582862039 1.5135846138000488
24 0.955156767122 1.7025775909423828
32 0.960174956648 2.0814244747161865
64 0.965277064616 3.602566957473755
128 0.967526986742 6.694079875946045
256 0.968645144161 11.87510895729065
512 0.969431999446 25.76741886138916
1024 0.969733055455 45.883561849594116"))
#2048 0.969820548189 88.62634110450745"))

oob_plot <- function(max_n=1024, plotboth=0, y='oob_estimate', ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino", plotlegend=T) {
    plot(data$n_trees, if (y=="oob_estimate") 1-data[[y]] else data[[y]], 
         axes=FALSE, pch=20, cex=1, tcl=0,
         xlab="Number of Decision Trees in Forest",
         ylab=if (y=='oob_estimate') "Out-of-Bag Error Rate" else 
              "Training time t/s",
         xaxs='i', yaxs='i', 
         log='xy',#if (y=='oob_estimate') 'x' else 'xy',
         ylim=if (y=='oob_estimate') range(10**-2, 1) else
             range(1, 60, data[[y]]))
    lines(data$n_trees, 
        if (y=='oob_estimate') 1-data[[y]] else data[[y]], type='l')
    minor <- 32
    axis(1, at=seq(minor, max_n, minor), labels=FALSE, tcl=0.125)
    axis(3, at=seq(minor, max_n, minor), labels=FALSE, tcl=0.125)
    
    axis(1, at=2^(2:log2(minor)), labels=FALSE, tcl=0.25)
    axis(3, at=2^(2:log2(minor)), labels=FALSE, tcl=0.25)
    
    axis(1, at=2^(2:log2(max_n)), labels=2^(2:log2(max_n)), tcl=0.25,
        cex.axis=text.cex)
    axis(3, at=2^(2:log2(max_n)), labels=FALSE, tcl=0.25)
    
    if (y == 'oob_estimate')
        magaxis(side=c(2,4), tcl=0.25, labels=c(1,0), mgp=mgp,#mgp+c(0,0.25,0), 
            family=font, #las=2, 
            cex.axis=text.cex)
    else {
        axis(2, at=seq(1, 60, 6), labels=FALSE, tcl=0.125)
        axis(4, at=seq(1, 60, 6), labels=FALSE, tcl=0.125)
        axis(2, at=c(1, 10, 60), labels=c(1, 10, 60), tcl=0.25, 
            cex.axis=text.cex)
        axis(4, at=c(1, 10, 60), labels=FALSE, tcl=0.25)
    }
    
    #abline(h=0.99, lty=3, col='darkred')
    #abline(v=data$n_trees[min(which(data$oob_estimate>0.99))],
    #       lty=3, col='darkred')
    ##lines(data3$n_trees, data3$oob_estimate, lty=3)
    if (plotboth) {
        points(data2$n_trees, 
            if (y=='oob_estimate') 1-data2[[y]] else data2[[y]], 
            col=blue, cex=0.5)
        lines(data2$n_trees, 
            if (y=='oob_estimate') 1-data2[[y]] else data2[[y]], 
            lty=3, col=blue)
        points(data3$n_trees, 
            if (y=='oob_estimate') 1-data2[[y]] else data3[[y]], 
            pch=3, col=red, cex=0.5)
        lines(data3$n_trees, 
            if (y=='oob_estimate') 1-data2[[y]] else data3[[y]], 
            lty=4, col=red)
    }
    if (plotlegend) {
        legend("topright", bty='n', xjust=1, yjust=1,
               lty=c(2,3,4), pch=c(20, 1, 3),
               col=c("black", blue, red), cex=text.cex,
               legend=c("Fit using all observables",
                        expression("Fit without"~log~g*", "*
                                   delta*nu[1*","*3]* ", and"~r[1*","*3]),
                        expression("Fit without"~L*", "*log~g*", "*
                                   delta*nu[1*","*3]* ", and"~r[1*","*3])))
    }
}

make_plots(oob_plot, "oob", plotlegend=F)
make_plots(oob_plot, "oob-without", plotboth=1)
make_plots(oob_plot, "time-without", plotboth=1, y='time', plotlegend=F)


ppt.DF <- read.table(file.path('subsets', 'points_per_track.dat'), header=1)



averages <- NULL
for (n_tracks in unique(DF$n_tracks)) {
    for (m_points in unique(DF$m_points)) {
        for (variable in unique(DF$variable)) {
            ss <- DF[DF$n_tracks == n_tracks & 
                     DF$m_points == m_points &
                     DF$variable == variable,]
            new.row <- ss[1,]
            new.row$cv <- mean(ss$cv)
            new.row$oob <- mean(ss$oob)
            averages <- rbind(averages, new.row)
        }
    }
}







legend2 <- legend
body(legend2)[[c(38,4,6)]] <- quote({
    cidx <- col(matrix(1, ncol=ncol, nrow=n.legpercol))[1:n.leg]
    wc <- tapply(rep_len(w0,n.leg), cidx, max)
    w <- sum(wc) + 0.5 * xchar
    w0 <- rep(cumsum(c(0,wc[-length(wc)])), each=n.legpercol)[1:n.leg]
})
body(legend2)[[40]] <- quote(xt <- left + xchar + xextra + w0)

plot_accuracy <- function(DF, score='cv', m.or.n='n_tracks', ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino", plotlegend=T) {
    
    if ('log_g' %in% DF$variable) {
        DF <- DF[-which(DF$variable == 'log_g'),]
    }
    
    par(xpd=FALSE)
    plot(NA, axes=F, 
         xlim=range(DF[[m.or.n]]), 
         ylim=c(1, 0.001),#rev(range(1-DF[[score]])),#c(0,1), 
         log='xy', xaxs='i', yaxs='i', 
         xlab=ifelse(m.or.n == 'n_tracks', 
                     "Number of Evolutionary Tracks",
                     "Models Per Track"),
         ylab="Accuracy")
    
    yticks <- c(1, 0.316, 0.1, 0.03, 0.01, 0.003, 0.001)
    axis(2, tick=F, at=yticks,
         labels=c("0%", "68%", "90%", "97%", "99%", "99.7%", "99.9%"))
    
    minor <- 2 ^ (floor(log(sqrt(max(DF[[m.or.n]]))) / log(2)))
    tick.places <- 2^(log2(min(DF[[m.or.n]])):log2(max(DF[[m.or.n]])))
    axis(1, at=tick.places, labels=tick.places, 
         tcl=-0.25, cex.axis=text.cex, tick=F)
    
    #grid(0, 5, lty=6, col="cornsilk2", equilogs=F) 
    for (ytick in yticks[1:(length(yticks)-1)]) abline(h=ytick, col='cornsilk2')
    for (tick in tick.places) abline(v=tick, col='cornsilk2')
    
    par(xpd=TRUE)
    variables <- unique(DF[['variable']])
    col.pal <- colorRampPalette(brewer.pal(11, 'Spectral'))(length(variables))
    for (variable_i in 1:length(variables)) {
        average <- NULL
        for (n_tracks in unique(DF[[m.or.n]])) {
            ss <- DF[DF[[m.or.n]] == n_tracks &
                         DF$variable == variables[variable_i],]
            new.row <- ss[1,]
            new.row[[score]] <- median(ss[[score]])
            average <- rbind(average, new.row)
        }
        vals <- ifelse(average[[score]] < 0, 0, average[[score]])
        vals <- ifelse(vals > 0.999, 0.999, vals)
        lines(average[[m.or.n]], 1-vals, col='black', lwd=3)
        lines(average[[m.or.n]], 1-vals, col=col.pal[variable_i], lwd=1.5)
        points(average[[m.or.n]], 1-vals, col='black', 
               pch=21+(variable_i%%3), bg=col.pal[variable_i], cex=0.75)
    }
    
    if (plotlegend) {
        shapes <- 1:length(variables)%%3 + 1
        labels <- as.expression(seis.labs[variables])
        text.widths <- rep(0.1, length(labels))
        #text.widths[4] <- 0.15
        text.widths[which(grepl('alpha|Y_surf', variables))] <- 0.175
        legend2("topleft", inset=c(-0.025, -0.11), horiz=T, bty='n',
               legend=labels, text.col='white', pch=c(16, 15, 18)[shapes], 
               col=col.pal, text.width=text.widths)
        legend2("topleft", inset=c(-0.025, -0.11), horiz=T, bty='n',
               legend=labels, pch=c(1, 0, 5)[shapes], col='black', 
               text.width=text.widths)
        #legend("topleft", inset=c(0, -0.15), horiz=T, bty='n',
        #       legend=labels[1:6], text.col='white', pch=15+shapes, 
        #       col=col.pal, text.width=0.05)
        #legend("topleft", inset=c(0, -0.15), horiz=T, bty='n',
        #       legend=labels[1:6], pch=0+shapes, col='black', text.width=0.05)
        #legend("topleft", inset=c(0, -0.075), horiz=T, bty='n',
        #       legend=labels[-1:-6], text.col='white', pch=15+shapes, 
        #       col=col.pal[-1:-6], text.width=0.05)
        #legend("topleft", inset=c(0, -0.075), horiz=T, bty='n',
        #       legend=labels[-1:-6], pch=0+shapes, col='black', text.width=0.05)
    }
}

make_plots(plot_accuracy, "num_tracks", 
    DF=read.table(file.path('subsets', 'num_tracks.dat'), 
                  header=T, stringsAsFactors=F), mar=c(3, 3, 1.5, 1))
make_plots(plot_accuracy, "num_points", m.or.n='m_points',
    DF=read.table(file.path('subsets', 'points_per_track.dat'), 
                  header=T, stringsAsFactors=F), mar=c(3, 3, 1.5, 1))













# plot cv as a function of num tracks
plot_accuracy <- function(DF, score='cv', m.or.n='n_tracks', ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino", plotlegend=T) {
    plot(NA, axes=F, xlim=range(DF[[m.or.n]]), ylim=c(0,1), 
         log='x', xaxs='i', yaxs='i', 
         xlab="Number of Evolutionary Tracks",
         ylab="Accuracy")
    for (variable in unique(DF$variable)) {
        average <- NULL
        for (n_tracks in unique(DF[[m.or.n]])) {
            ss <- DF[DF[[m.or.n]] == n_tracks &
                     DF$variable == variable,]
            new.row <- ss[1,]
            new.row[[score]] <- mean(ss[[score]])
            average <- rbind(average, new.row)
        }
        points(average[[m.or.n]], average[[score]])
        lines(average[[m.or.n]], average[[score]])
    }
    magaxis(side=c(2,4), tcl=-0.25, labels=c(1,0), mgp=mgp,
            family=font, las=1, cex.axis=text.cex)
    minor <- 
    minor <- 2 ^ (floor(log(sqrt(max(DF[[m.or.n]]))) / log(2)))
    axis(1, at=seq(minor, max(DF[[m.or.n]]), minor), labels=FALSE, tcl=-0.125)
    axis(3, at=seq(minor, max(DF[[m.or.n]]), minor), labels=FALSE, tcl=-0.125)
    
    axis(1, at=2^(2:log2(minor)), labels=FALSE, tcl=-0.25)
    axis(3, at=2^(2:log2(minor)), labels=FALSE, tcl=-0.25)
    
    axis(1, at=2^(2:log2(max(DF[[m.or.n]]))), 
        labels=2^(2:log2(max(DF[[m.or.n]]))), 
        tcl=-0.25, cex.axis=text.cex)
    axis(3, at=2^(2:log2(max(DF[[m.or.n]]))), 
        labels=FALSE, tcl=-0.25)
}
make_plots(plot_accuracy, "num_tracks", DF=DF)
make_plots(plot_accuracy, "num_points", DF=DF, m.or.n='m_points')

# plot oob as a function of num tracks 
for (variable in unique(DF$variable)) {
    
}

# plot cv as a function of num points 
for (variable in unique(DF$variable)) {
    
}

# plot oob as a function of num points 
for (variable in unique(DF$variable)) {
    
}
