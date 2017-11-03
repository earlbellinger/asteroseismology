#### Plot kernel functions of a stellar model 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
source(file.path('..', 'scripts', 'utils.R')) 
source('models.R') 
source('frequencies.R') 
source('OLA_invert.R') 
models <- get_model_list()
perturb <- F
k.pair <- NULL
col.pal <- c(NA, "#283845", "#F29559")

ref.mod <- 'diffusion'
model <- get_model(freqs=NULL, model.name='diffusion', target.name=NULL, 
    k.pair=NULL, square.Ks=F) 

MDI.freqs <- get_freqs(target.name='MDI') 
BiSON.freqs <- get_freqs(target.name='BiSON') 
CygA.freqs <- get_freqs(target.name='BiSON', mode.set='CygA') 

freqs.list <- list(MDI.freqs, BiSON.freqs, CygA.freqs)

for (ii in 1:length(freqs.list)) {
    freqs <- freqs.list[[ii]]
    rs <- model$r * model$R * solar_radius
    freqs.list[[ii]]$r_ts <- apply(freqs, 1, function(freq) {
        ell <- freq[['l']]
        nu <- freq[['nu']]
        if (ell == 0) return(0)
        model$r[which.min((
                model$cs.spl(rs / (model$R * solar_radius))**2 / (rs**2) 
                - (10**-6*nu*(2*pi))**2/(ell*(ell+1))
            )**2)]
    })
}

plot_turning_points <- function(freqs.list, ..., 
        text.cex=1, mgp=utils.mgp, font="Palatino Linotype", 
        thin=F, mar=utils.mar) {
    
    par(mar=mar-c(1, 1, 0, 0))
    #par(mar=mar-c(0.5, 0, 0, 0), mgp=mgp-c(0.2, 0, 0))
    
    plot(NA, axes=F, #frame.plot=T, 
        xlim=range(sapply(freqs.list, function(freqs) range(freqs$nu))),
        ylim=c(0, 1),
        xlab="",
        ylab=expression("Lower turning point"~r[t]/R ))
    
    for (ii in 1:length(freqs.list)) {
        nus <- freqs.list[[ii]]
        
        for (ell in rev(sort(unique(nus$l)))) {
            nus.ell <- nus[nus$l==ell,]
            if (nrow(nus.ell) <= 1) next 
            if (any(nus.ell$r_ts < 0.999)) {
                nu.seq <- with(nus.ell, seq(min(nu), max(nu), 1))
                r_ts.seq <- splinefun(nus.ell$nu, nus.ell$r_ts)(nu.seq)
                lines(nu.seq, r_ts.seq, col="gray", lty=1, lwd=1)
            }
        }
        for (ell in rev(sort(unique(nus$l)))) {
            nus.ell <- nus[nus$l==ell,]
            cex <- ifelse(nus.ell$l < 4 | nus.ell$r_ts <= 0.2, 
                    0.6,  
                    0.6*log(nus.ell$r_ts) / log(0.2))
            points(nus.ell$nu, nus.ell$r_ts, 
                cex=cex,
                col=1, lwd=0.5, 
                bg=col.pal[ii], 
                pch=if (ell < 4) c(22,24,23,21)[ell+1] else 1)
        }
        
    }
    
    text(x=4500, y=0.85, labels='MDI', cex=text.cex)#, font=font)
    text(x=1650, y=0.042, labels='BiSON', cex=text.cex)#, font=font)
    text(x=4400, y=0.01, labels=substitute(paste(italic('Kepler'))), 
        cex=text.cex)#, font=font)
    
    text(x=972.613, y=0.038, labels='0', cex=text.cex)
    text(x=1010, y=0.139, labels='1', cex=text.cex)
    text(x=1220, y=0.191, labels='2', cex=text.cex)
    text(x=1430, y=0.225, labels='3', cex=text.cex)
    text(x=1195, y=0.28, labels='\u2113 \u2265 4', cex=text.cex)
    
    
    
    magaxis(side=c(1,3,4), tcl=0.25, labels=c(1,0,0), 
            las=thin, mgp=mgp-c(0.5,0.25,0), 
            family=font, cex.axis=text.cex)
    magaxis(side=2, tcl=0.25, labels=T, 
            las=thin, mgp=mgp-c(0,0,0), 
            family=font, cex.axis=text.cex)
    
    par(mgp=mgp-c(0.5, 0, 0))
    title(xlab=expression("Frequency"~nu/mu*Hz))
    
}

#plot_turning_points(freqs.list, font='Palatino Linotype')

make_plots(plot_turning_points, 
    paste0('turning_points'), 
    filepath=file.path('plots', 'echelle'), 
    wide=F, short=F, make_png=T, slides=F, 
    use.cairo=T, font='Palatino Linotype', 
    freqs.list=freqs.list, cex.paper=0.7)

