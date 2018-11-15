#### HR diagram for two tracks of M = 1, 2
#### Author: Earl Bellinger ( bellinger@phys.au.dk )
#### Stellar Astrophysics Centre Aarhus

source(file.path('..', 'scripts', 'utils.R'))

DF <-
    read.table(file.path('examples', '1', 'LOGS', 'history.data'),
        header = 1,
        skip = 5)
DF2 <-
    read.table(file.path('examples', '2', 'LOGS', 'history.data'),
        header = 1,
        skip = 5)
DF3 <-
    read.table(file.path('examples', '3', 'LOGS', 'history.data'),
        header = 1,
        skip = 5)

plot_HR <- function(...,
    text.cex = 1,
    mgp = utils.mgp,
    mar = utils.mar,
    font = utils.font,
    tcl = utils.tcl) {
    par(
        mar = mar + c(0.3, -0.5, 1.5, 0.1),
        lwd = 1.5,
        las = 1,
        cex.axis = text.cex
    )
    
    xlim <- c(6000, 4500)
    ylim <- c(0.1, 30)
    
    plot(
        NA,
        axes = F,
        log = 'xy',
        xaxs = 'i',
        yaxs = 'i',
        xlim = xlim,
        ylim = ylim,
        xlab = "",
        ylab = ""
    )
    
    spectral.divs <- c(30000, 10000, 7500, 6000, 5200, 3700, 2400)
    rs <- c(175 / 255, 199 / 255, 1, 1, 1, 1, 1, 1)
    gs <- c(201, 216, 244, 229, 217, 199, 166) / 255
    bs <-
        c(1, 1, 243 / 255, 207 / 255, 178 / 255, 142 / 255, 81 / 255)
    cols <- c(
        rgb(175 / 255, 201 / 255, 1),
        # O
        rgb(199 / 255, 216 / 255, 1),
        # B
        rgb(1,       244 / 255, 243 / 255),
        # A
        rgb(1,       229 / 255, 207 / 255),
        # F
        rgb(1,       217 / 255, 178 / 255),
        # G
        rgb(1,       199 / 255, 142 / 255),
        # K
        rgb(1,       166 / 255, 81 / 255)
    )  # M
    #if (F) {
    for (ii in 1:length(spectral.divs)) {
        div <- spectral.divs[ii]
        if (div >= xlim[1])
            next
        if (div < xlim[2])
            div <- xlim[2]
        if (ii == 1) {
            rect(xlim[1],
                ylim[1],
                div,
                ylim[2],
                col = cols[ii],
                border = NA)
        } else {
            prev <- spectral.divs[ii - 1]
            if (prev > xlim[1])
                prev <- xlim[1]
            rect(prev,
                ylim[1],
                div,
                ylim[2],
                col = cols[ii],
                border = NA)
        }
    }
    for (ii in 2:(length(spectral.divs) - 1)) {
        div <- spectral.divs[ii]
        gradient.rect(
            div + 0.0025,
            ylim[1],
            div - 0.0025,
            ylim[2],
            nslices = 10,
            border = NA,
            reds = c(rs[ii], rs[ii + 1]),
            greens = c(gs[ii], gs[ii + 1]),
            blues = c(bs[ii], bs[ii + 1])
        )
    }
    
    with(DF,  lines(10 ** log_Teff, 10 ** log_L, lwd = 1.5, col=blue,  lty=1))
    with(DF2, lines(10 ** log_Teff, 10 ** log_L, lwd = 1.5, col=red, lty=2))
    with(DF3, lines(10 ** log_Teff, 10 ** log_L, lwd = 1.5,           lty=3))
    
    with(DF,  points(10 ** log_Teff[1], 10 ** log_L[1], lwd = 1, pch=20))
    with(DF2, points(10 ** log_Teff[1], 10 ** log_L[1], lwd = 1, pch=20))
    with(DF3, points(10 ** log_Teff[1], 10 ** log_L[1], lwd = 1, pch=20))
    
    
    par(family='Helvetica LT Std Light')
    text(5772, 1,    'Sun', pos=3,                        cex=0.9*text.cex)
    text(5070, 10,  expression(beta == 0.01),  col=1,    cex=0.9*text.cex)
    text(4780, 6.5,  expression(beta == -0.01), col=red, cex=0.9*text.cex)
    #text(5270, 19.5, 'Pre-main sequence',                 cex=0.9*text.cex)
    text(5060, 19.6, expression(tau == 0),                cex=0.9*text.cex)
    text(5075, 1.85, 'Hayashi',                           cex=0.9*text.cex)
    text(5075, 1.4,  'track',                             cex=0.9*text.cex)
    text(4950, 0.23, 'Henyey track',                      cex=0.9*text.cex)
    text(5500, 0.23, 'ZAMS',                              cex=0.9*text.cex)
    text(5800, 0.64, expression(tau == 10~Gyr),           cex=0.9*text.cex)
    #text(5850, 0.65, 'TAMS',                              cex=0.9*text.cex)
    par(family=font)
    
    ## solar symbol
    points(5772,
        1,
        pch = 20,
        cex = 1.1 / 2,
        lwd = 1.5)
    points(5772,
        1,
        pch = 1,
        cex = 1.1,
        lwd = 1.5)
    
    par(xpd = NA)
    rect(xlim[2],
        ylim[1] * 0.1,
        xlim[2] - 100,
        ylim[2] * 10,
        col = 'white',
        border = NA)
    rect(xlim[1] * 1.1,
        ylim[1],
        xlim[2] * 0.9,
        ylim[1] * 0.9,
        col = 'white',
        border = NA)
    par(xpd = F)
    
    nxticks <- 4
    nyticks <- 4
    nxminor <- 4
    nyminor <- 4
    xticks <- pretty(xlim, n = nxticks)
    #yticks <- pretty(ylim, n=nyticks)
    xticks.minor <- pretty(xlim, n = nxticks * nxminor)
    #yticks.minor <- pretty(ylim, n=nyticks*nyminor)
    xticks.minor <- xticks.minor[!xticks.minor %in% xticks]
    #yticks.minor <- yticks.minor[!yticks.minor %in% yticks]
    par(mgp = mgp + c(0, 0.25, 0))
    #xpos <- seq(10**xlim[2], 10**xlim[1], 2000)
    #xpos2 <- seq(10**xlim[2], 10**xlim[1], 500)
    xpos <- xticks#seq(xlim[2], xlim[1], 1000)
    xpos2 <- xticks.minor#seq(xlim[2], xlim[1], 200)
    #axis(side=1, tcl=tcl/2, at=log10(xpos2), labels=F, lwd.ticks=par()$lwd)
    axis(
        side = 1,
        tcl = tcl / 2,
        at = xpos2,
        labels = F,
        lwd.ticks = par()$lwd
    )
    #axis(side=1, tcl=tcl, at=log10(xpos), labels=xpos, cex.axis=text.cex,
    axis(
        side = 1,
        tcl = tcl,
        at = xpos,
        labels = xpos,
        cex.axis = text.cex,
        lwd.ticks = par()$lwd
    )
    par(mgp = mgp + c(0, 0.43, 0))
    magaxis(
        2,
        tcl = tcl,
        mgp = mgp + c(0, 0.4, 0),
        cex.axis = text.cex,
        family = font,
        las = 1,
        majorn = 3,
        labels = T,
        lwd.ticks = par()$lwd
    )
    
    box(lwd = par()$lwd)
    
    mtext(
        expression(L / L["solar"]),
        2,
        2,
        outer = F,
        las = 0,
        cex = 1.15*text.cex
    )
    mtext(expression(T["eff"] / K),
        1,
        2,
        outer = F,
        cex = 1.15*text.cex)
    
    mtext(
        expression("Spectral Type"),
        side = 3,
        line = 1.3,
        cex = 1.15*text.cex
    )
    
    spectral.labs <- c("O", "B", "A", "F", "G", "K", "M")
    selector <- 1:(length(spectral.divs))
    spectral.Teffs <- sapply(selector,
        function(ii) {
            div <- spectral.divs[ii]
            if (div >= xlim[1])
                return(Inf)
            if (div < xlim[2])
                div <- xlim[2]
            if (ii == 1)
                return((xlim[1] + div) / 2)
            prev <- spectral.divs[ii - 1]
            if (prev > xlim[1])
                prev <- xlim[1]
            (div + prev) / 2
        })
    axis(
        3,
        at = spectral.divs,
        tcl = tcl,
        labels = F,
        cex.axis = text.cex,
        lwd.ticks = par()$lwd
    )
    par(mgp = mgp + c(0, -0.05, 0))
    axis(
        3,
        at = spectral.Teffs,
        labels = spectral.labs[selector],
        cex.axis = text.cex,
        tcl = 0
    )
}

make_plots(
    plot_HR,
    'HR',
    paper_pdf_height = 4.17309 * 0.8,
    #cex.paper = 0.95,
    use.cairo = T,
    font = 'Palatino Linotype'
)

