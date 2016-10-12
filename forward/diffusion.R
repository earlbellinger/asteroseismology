#### Plot and animate Kippenhahns of solar-calibrated stars varied by diffusion
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source(file.path('..', 'scripts', 'utils.R'))

library(magicaxis)
library(parallel)
library(parallelMap)
library(RColorBrewer)

sim.dirs <- file.path('simulations-diffusion')
dirs <- list.dirs(sim.dirs, recursive=F)

#pro.file <- profile.files[1]
#early <- as.numeric(sub('.data', '', sub('profile', '', profile.files)))<100
h1.vals <- c(0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
m.vals <- seq(0.01, 1, 0.01)

parallelStartMulticore(32)

get_diffusion <- function(directory) {
    params.DF <- NULL
    for (var in unlist(strsplit(basename(directory), '_'))) { 
        nameval <- unlist(strsplit(var, "=")) 
        params.DF[nameval[1]] <- as.numeric(nameval[2])
    }
    diffusion <- params.DF[['diffusion']]
    if (diffusion > 90) diffusion <- 100
    sprintf('%.1f', diffusion)
    #round(params.DF[['diffusion']], -1)
}

### Kippenhahn of eps_nuc/age/mass
nuclear_Kippenhahn <- function(log.dir, profile.files, diffusion, ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    eps <- do.call(rbind, parallelMap(function(pro.file) {
            pro.head <- read.table(file.path(log.dir, pro.file), header=1, 
                nrow=1, skip=1)
            age <- pro.head$star_age/10**9
            pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
            eps_nuc <- pro.DF$eps_nuc
            cutoff <- eps_nuc > 0.1
            eps_nuc <- log10(eps_nuc[cutoff])
            mass <- pro.DF$mass[cutoff]
            contours <- approx(mass, eps_nuc, m.vals, rule=2)$y
            data.frame(age, rbind(contours))
        }, profile.files))#[early]))
    eps <- eps[order(eps$age),]
    
    conv <- do.call(rbind, parallelMap(function(pro.file) {
            pro.head <- read.table(file.path(log.dir, pro.file), header=1, 
                nrow=1, skip=1)
            age <- pro.head$star_age/10**9
            pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
            grad.diff <- pro.DF$gradr - pro.DF$grada
            mass <- pro.DF$mass
            contours <- approx(mass, grad.diff, m.vals, rule=2)$y
            data.frame(age, rbind(contours))
        }, profile.files))#[early]))
    conv <- conv[order(conv$age),]
    
    filled.contour(eps$age, log10(m.vals), as.matrix(eps[,-1]),
        color=colorRampPalette(brewer.pal(11, "Spectral")),
        zlim=c(-1, 2.1),
        xlim=c(0, 10),
        xaxs='i', yaxs='i',
        key.axes={
            axis(4, tcl=0, line=0)
            mtext(expression("Nuclear Energy"~log[10](epsilon/erg/g/s)), 
                side=4, las=3, line=3)
        },
        plot.axes={
            contour(conv$age, log10(m.vals), as.matrix(conv[,-1]), levels=c(0),
                drawlabels=F, col='black', add=TRUE)
            abline(v=max(eps$age), lty=2)
            magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0),
               unlog='y')
            legend("topleft", bty='n', text.col='white', inset=c(-0.05, 0),
                legend=as.expression(bquote(D==.(diffusion))))
        },
        plot.title={
            title(xlab=expression("Star age"~tau/"Gyr"), line=2)
            title(ylab=expression(m/M["*"]), line=2)
        }
    )
}

just.contour <- function(log.dir, profile.files, diffusion, ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    filled.contour(eps$age, log10(m.vals), as.matrix(eps[,-1]),
        color=colorRampPalette(brewer.pal(11, "Spectral")),
        zlim=c(-1, 2.1),
        xlim=c(0, 10),
        xaxs='i', yaxs='i',
        key.axes={
            axis(4, tcl=0, line=0)
            mtext(expression("Nuclear Energy"~log[10](epsilon/erg/g/s)), 
                side=4, las=3, line=3)
        },
        plot.axes={
            contour(conv$age, log10(m.vals), as.matrix(conv[,-1]), levels=c(0),
                drawlabels=F, col='black', add=TRUE)
            abline(v=max(eps$age), lty=2)
            magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0),
               unlog='y')
            legend("topleft", bty='n', text.col='white', inset=c(-0.05, 0),
                legend=as.expression(bquote(D==.(diffusion))))
        },
        plot.title={
            title(xlab=expression("Star age"~tau/"Gyr"), line=2)
            title(ylab=expression(m/M["*"]), line=2)
        }
    )
}

hydro_Kippenhahn <- function(log.dir, profile.files, diffusion, ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    h1s <- do.call(rbind, parallelMap(function(pro.file) {
            pro.head <- read.table(file.path(log.dir, pro.file), header=1, 
                nrow=1, skip=1)
            age <- pro.head$star_age/10**9
            pro.DF <- read.table(file.path(log.dir, pro.file), header=1, skip=5)
            radius <- pro.DF$radius
            h1 <- pro.DF$h1
            contours <- approx(radius, h1, m.vals, rule=2)$y
            data.frame(age, rbind(contours))
        }, profile.files))
    h1s <- h1s[order(h1s$age),]
    
    contour(h1s$age, log10(m.vals), as.matrix(h1s[,-1]), levels=h1.vals, axes=F,
        xlim=c(0, 10), xaxs='i')
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0),
           unlog='y')
    abline(v=max(h1s$age), lty=2)
    title(xlab=expression("Star age"~tau/"Gyr"), line=2)
    title(ylab=expression(r/R["*"]), line=2)
    legend("topleft", legend=as.expression(bquote(D==.(diffusion))), bty='n',
        inset=c(-0.05, 0))
}

surface_composition <- function(history.file, diffusion, ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    age <- history.file$star_age / 10**9
    plot(NA, axes=F, xaxs='i', yaxs='i', log='y', 
         xlim=c(0, 10), 
         ylim=c(10^-6, 1),
         xlab=expression("Star age"~tau/"Gyr"),
         ylab=expression("Surface Abundance"))
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0), 
        unlog='y')
    pal <- brewer.pal(8, 'Spectral')
    lines(age, history.file$surface_h1, lwd=3, col=pal[1])
    lines(age, history.file$surface_he3, lwd=3, col=pal[2])
    lines(age, history.file$surface_he4, lwd=3, col=pal[3])
    lines(age, history.file$surface_c12, lwd=3, col=pal[4])
    lines(age, history.file$surface_n14, lwd=3, col=pal[5])
    lines(age, history.file$surface_o16, lwd=3, col=pal[6])
    lines(age, history.file$surface_ne20, lwd=3, col=pal[7])
    lines(age, history.file$surface_mg24, lwd=3, col=pal[8])
    abline(v=max(age), lty=2)
    abline(v=4.57, lty=3)
    legend("bottomright", col=pal, lwd=3, xpd=1, cex=0.75,
        #horiz=1, #pch=20, #bty='n', 
        #inset=c(-0.2, 0),
        inset=c(0.015, 0.01),
        #text.width=0.3, 
        #cex=0.75, 
        legend=expression(""^1*"H", ""^3*"He", ""^4*"He", ""^12*"C", ""^14*"N",
            ""^16*"O", ""^20*"Ne", ""^24*"Mg"))
    legend("topleft", legend=as.expression(bquote(D==.(diffusion))), bty='n',
        inset=c(-0.05, 0))
}

hrd <- function(history.file, diffusion, 
        L_range=c(0.4, 2.1), Teff_range=c(6000, 5000),
        L_true=1, Teff_true=5777, 
        L_unc=0, Teff_unc=0,
        age_true=4.57, age_unc=0,
        ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    age <- history.file$star_age / 10**9
    plot(NA, axes=F, xaxs='i', yaxs='i', #log='xy', 
         ylim=L_range, 
         xlim=Teff_range, 
         xlab=expression("Temperature"~T["eff"]/K),
         ylab=expression("Luminosity"~L/L["☉"]))
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0))
    pal <- brewer.pal(8, 'Spectral')
    
    if (Teff_unc == 0) {
        abline(h=L_true, lty=3)
        abline(v=Teff_true, lty=3)
        points(Teff_true, L_true, cex=0.1, pch=20)
        points(Teff_true, L_true)
    } else {
        polygon(c(10000, 10000, 0, 0),
                c(L_true-L_unc, L_true+L_unc,
                  L_true+L_unc, L_true-L_unc),
                col="#00000010", border=NA)
        polygon(c(Teff_true-Teff_unc, Teff_true-Teff_unc, 
                  Teff_true+Teff_unc, Teff_true+Teff_unc),
                c(-100, 100, 100, -100),
                col="#00000010", border=NA)
    }
    
    lines(10**history.file$log_Teff, 10**history.file$log_L, lwd=2)
    
    ages <- if (age_unc == 0) {
        age_true 
    } else {
        seq(age_true-age_unc, age_true+age_unc, 0.01)
    }
    T_age <- approx(age, 10**history.file$log_Teff, ages)$y
    L_age <- approx(age, 10**history.file$log_L, ages)$y
    lines(T_age, L_age, lwd=2.5, col='darkred')
    
    legend("topleft", legend=as.expression(bquote(D==.(diffusion))), bty='n',
        inset=c(-0.05, 0))
}

temp_feh <- function(history.file, diffusion, 
        FeH_range=c(0.4, 2.1), Teff_range=c(6000, 5000),
        FeH_true=0, Teff_true=5777, 
        FeH_unc=0, Teff_unc=0,
        age_true=4.57, age_unc=0,
        ...,
        text.cex=1, font=utils.font, mgp=utils.mgp) {
    
    age <- history.file$star_age / 10**9
    metals <- FeH(10**history.file$log_surf_cell_z, history.file$center_h1)
    
    plot(NA, axes=F, xaxs='i', yaxs='i', #log='xy', 
         ylim=FeH_range, 
         xlim=Teff_range, 
         xlab=expression("Temperature"~T["eff"]/K),
         ylab=expression("["*Fe/H*"]"))
    magaxis(side=1:4, family='Palatino', tcl=-0.25, labels=c(1,1,0,0))
    pal <- brewer.pal(8, 'Spectral')
    
    if (Teff_unc == 0) {
        abline(h=FeH_true, lty=3)
        abline(v=Teff_true, lty=3)
        points(Teff_true, FeH_true, cex=0.1, pch=20)
        points(Teff_true, FeH_true)
    } else {
        polygon(c(10000, 10000, -10000, -10000),
                c(L_true-L_unc, L_true+L_unc,
                  L_true+L_unc, L_true-L_unc),
                col="#00000010", border=NA)
        polygon(c(Teff_true-Teff_unc, Teff_true-Teff_unc, 
                  Teff_true+Teff_unc, Teff_true+Teff_unc),
                c(-10000, 10000, 10000, -10000),
                col="#00000010", border=NA)
    }
    
    lines(10**history.file$log_Teff, metals, lwd=2)
    
    ages <- if (age_unc == 0) {
        age_true 
    } else {
        seq(age_true-age_unc, age_true+age_unc, 0.01)
    }
    T_age <- approx(age, 10**history.file$log_Teff, ages)$y
    FeH_age <- approx(age, metals, ages)$y
    lines(T_age, FeH_age, lwd=2.5, col='darkred')
    
    legend("topleft", legend=as.expression(bquote(D==.(diffusion))), bty='n',
        inset=c(-0.05, 0))
}

parallelMap(function(sim.dir) {
        log.dir <- file.path(sim.dir, 'LOGS')
        log.files <- list.files(log.dir)
        profile.files <- log.files[grepl('profile\\d+\\.data$', log.files)]
        
        history.file <- read.table(file.path(log.dir, 'history.data'),
            header=1, skip=5)
        
        make_plots(just.contour, 
                paste0('nuclear-', get_diffusion(sim.dir)), 
                filepath=file.path('plots', 'diffusion', 'Kippenhahn'), 
                log.dir=log.dir, 
                profile.files=profile.files, diffusion=get_diffusion(sim.dir), 
                thin=F, paper=F, mar=c(4, 5, 1, 7))
        
        make_plots(surface_composition, 
                paste0('surface-', get_diffusion(sim.dir)), 
                filepath=file.path('plots', 'diffusion', 'Kippenhahn'), 
                log.dir=log.dir, 
                history.file=history.file, diffusion=get_diffusion(sim.dir), 
                paper=F, mar=c(3, 4, 1, 1))
        
        
        
        
        make_plots(nuclear_Kippenhahn, 
                paste0('nuclear-', get_diffusion(sim.dir)), 
                filepath=file.path('plots', 'diffusion', 'Kippenhahn'), 
                log.dir=log.dir, 
                profile.files=profile.files, diffusion=get_diffusion(sim.dir), 
                thin=F, short=F, paper=F, mar=c(4, 5, 1, 7))

        make_plots(hydro_Kippenhahn, 
                paste0('hydro-', get_diffusion(sim.dir)), 
                filepath=file.path('plots', 'diffusion', 'Kippenhahn'), 
                log.dir=log.dir, 
                profile.files=profile.files, diffusion=get_diffusion(sim.dir),
                thin=F, short=F, paper=F, make_pdf=F, mar=c(4, 5, 1, 2))
    }, dirs[c(1, 579, 581, 655, 754, 676, 824)])#dirs)












Map(function(sim.dir, diffusion) {
        log.dir <- file.path(sim.dir, 'LOGS')
        log.files <- list.files(log.dir)
        profile.files <- log.files[grepl('profile\\d+\\.data$', log.files)]
        
        print(log.dir)
        
        history.file <- read.table(file.path(log.dir, 'history.data'),
            header=1, skip=5)
        
        #print(range(10**history.file$log_L))
        #print(range(10**history.file$log_Teff))
        #return(1)
        
        make_plots(surface_composition, 
                paste0('surface-', sprintf('%4.2f', diffusion)), 
                filepath=file.path('plots', 'diffusion3', 'Kippenhahn'), 
                log.dir=log.dir, 
                history.file=history.file, diffusion=diffusion, 
                paper=F, mar=c(3, 4, 1, 1))
            
        make_plots(hrd, 
                paste0('hrd-', sprintf('%4.2f', diffusion)), 
                filepath=file.path('plots', 'diffusion3', 'Kippenhahn'), 
                log.dir=log.dir, 
                history.file=history.file, diffusion=diffusion, 
                #L_true=2.35, Teff_true=5622, 
                #L_unc=0.12, Teff_unc=106, 
                #L_range=c(0.5, 2.75), Teff_range=c(6100, 5000), 
                L_true=0.344, Teff_true=5046, 
                L_unc=0.022, Teff_unc=74, 
                L_range=c(0, 0.8), Teff_range=c(4650, 5300), 
                age_true=10.3, age_unc=0.96,
                paper=F, mar=c(3, 4, 1, 1))
            
        make_plots(temp_feh, 
                paste0('temp_feh-', sprintf('%4.2f', diffusion)), 
                filepath=file.path('plots', 'diffusion3', 'Kippenhahn'), 
                log.dir=log.dir, 
                history.file=history.file, diffusion=diffusion, 
                #L_true=2.35, Teff_true=5622, 
                #L_unc=0.12, Teff_unc=106, 
                #L_range=c(0.5, 2.75), Teff_range=c(6100, 5000), 
                FeH_true=-0.37, Teff_true=5046, 
                FeH_unc=0.09, Teff_unc=74, 
                FeH_range=c(-.5, 3), Teff_range=c(4650, 5300), 
                age_true=10.3, age_unc=0.96,
                paper=F, mar=c(3, 4, 1, 1))
        
        #make_plots(nuclear_Kippenhahn, 
        #        paste0('nuclear-', sprintf('%4.2f', diffusion)), 
        #        filepath=file.path('plots', 'diffusion', 'Kippenhahn'), 
        #        log.dir=log.dir, 
        #        profile.files=profile.files, diffusion=diffusion, 
        #        thin=F, short=F, make_pdf=F, paper=F, mar=c(4, 5, 1, 7))
        
    }, list.dirs(file.path('simulations2-6278762'), recursive=F),
       diffusion=c(0, 0.05, 0.5, 1, 10, 100, 2, 25, 5, 5.23, 50))
       
       #dirs[c(1, 416, 540, 579, 581, 655, 753, 754, 676, 824)],
       #diffusion=c(0, 0.05, 0.5, 1, 10, 2, 50, 5, 25, 100))

