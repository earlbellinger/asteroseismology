#### Mesh and scatterplot analysis of evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

## Libraries
source(file.path('..', 'scripts', 'utils.R'))
library(magicaxis)
library(RColorBrewer)
library(akima)
library(parallel)
library(parallelMap)
library(data.table)
library(lattice)

## Load data
seis.DF <- data.table(read.table('grid.dat', header=1))
setkey(seis.DF, M, Y, Z, alpha)
keys <- key(seis.DF)

solar_vals <- read.table(
        file.path('..', 'inverse', 'perturb', 'Sun_perturb.dat'), 
    nrow=1, header=1)

col.pal <- colorRampPalette(brewer.pal(11, "Spectral"))(1000)

# Sort data
combos <- unique(seis.DF[,keys, with=0])
ages <- unlist(Map(function(i) max(merge(seis.DF, combos[i,])$age), 
    1:nrow(combos)))
combos <- combos[order(ages),]

## Labels
Z_names <- names(seis.DF)[1:9]

Z_labels <- list(
    M=expression(M/M["\u0298"]), 
    Y=expression(Y[0]), 
    Z=expression(Z[0]),
    alpha=expression(alpha["MLT"]), 
    age='Age/Gyr', 
    radius=expression(R/R["\u0298"]), 
    H=expression(X(H)), 
    He=expression(X(He)), 
    Hc=expression(H[c])
)

Z_levels <- list(
    M=seq(0.7, 1.3, 0.1),
    Y=seq(0.22, 0.34, 0.01),
    Z=log10(seq(10**1e-04, 10**0.04, length=10)),
    alpha=seq(1.5, 2.5, 0.1),
    age=0:14, 
    radius=seq(0.6, 2.1, 0.1),
    H=seq(0.54, 0.78, 0.02),
    He=seq(0.22, 0.45, 0.02),
    Hc=seq(0, 0.78, 0.05)
)

color_levels <- list(
    age=seq(0, 13.8, 0.5),
    M=seq(0.7, 1.3, 0.025),
    Y=seq(0.22, 0.34, 0.0025),
    Z=log10(seq(10**1e-04, 10**0.04, length=20)),
    alpha=seq(1.5, 2.5, 0.025),
    radius=seq(0.6, 2.1, 0.025),
    H=seq(0.54, 0.78, 0.01),
    He=seq(0.22, 0.45, 0.005),
    Hc=seq(0, 0.78, 0.015)
)

seis.labs <- expression(
    M, Y[0], Z[0], alpha["MLT"], 
    tau, R, "X(H)", "X(He)", H[c], 
    log~g, L, T["eff"], "Fe"/"H", 
    "<"*Delta*nu*">", "<"*d*Delta*nu/d*nu*">", 
    "<"*Delta*nu[0]*">", "<"*d*Delta*nu[0]/d*nu*">", 
    "<"*delta*nu[0*","*2]*">", "<"*d*delta*nu[0*","*2]/d*nu*">", 
    "<"*r[0*","*2]*">", "<"*d*r[0*","*2]/d*nu*">", 
    "<"*r[0*","*1]*">", "<"*d*r[0*","*1]/d*nu*">", 
    "<"*delta*nu[1*","*3]*">", "<"*d*delta*nu[1*","*3]/d*nu*">",  
    "<"*r[1*","*3]*">", "<"*d*r[1*","*3]/d*nu*">",
    "<"*r[1*","*0]*">", "<"*d*r[1*","*0]/d*nu*">"
)

seis.latex <- c(
    "M", "$Y_0$", "$Z_0$", "$\\alpha_{\\text{\"MLT\"}}$", 
    "$\\tau$", "mass", "R", "X(H)", "X(Hc)", "$H_c$", 
    "$\\log g$", "L", "$T_{\text{\"eff\"}}$", "Fe/H", 
    
    "$\\langle\\Delta\\nu\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu}{d\nu}\\rangle$", 
    
    "$\\langle\\Delta\\nu_0\\rangle$", 
    "$\\langle\\frac{d\\Delta\\nu_0}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{02}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{02}\\rangle$", 
    "$\\langle\\frac{dr_{02}}{d\nu}\\rangle$", 
    
    "$\\langle r_{01}\\rangle$", 
    "$\\langle\\frac{dr_{01}}{d\nu}\\rangle$", 
    
    "$\\langle\\delta\\nu_{13}\\rangle$", 
    "$\\langle\\frac{d\\delta\\nu_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{13}\\rangle$", 
    "$\\langle\\frac{dr_{13}}{d\nu}\\rangle$", 
    
    "$\\langle r_{10}\\rangle$", 
    "$\\langle\\frac{dr_{10}}{d\nu}\\rangle$"
)

## Plot histograms
#tmp <- data.frame(seis.DF[,1:8, with=0])
#colnames(tmp) <- seis.labs[-exclude][1:8]
#d <- melt(tmp)
#ggplot(d, aes(x = value)) +#, y=..density..)) +
#    geom_histogram(aes(y=..ncount..), fill="#c0392b", alpha=0.75) + 
#    fte_theme() + 
#    geom_density(aes(y = ..scaled..)) +#col=2) + 
#    scale_y_continuous(labels=comma) + 
#    geom_hline(yintercept=0, size=0.4, color="black") +
#    facet_wrap(~variable,scales = "free_x", nrow=4) +
#    ggtitle(seis.labs[-exclude][1:8])


# Make inputs diagram
inputs_plot <- function(...) {
    H <- 1-combos$Y-combos$Z
    varmax <- max(H)
    varmin <- min(H)
    splom(combos,
        cex=0.001, pch=3,
        col=col.pal[floor((H-varmin)/(varmax-varmin)*(length(col.pal)-1))+1],
        xlab=NULL, ylab=NULL, 
        axis.text.cex=0.25,
        axis.text.lineheight=0.0001,
        axis.line.tck=0.25,
        xaxs='n', yaxs='n',
        varname.cex=0.5, varnames=Z_labels[1:4])
}
make_plots(inputs_plot, "inputs", filepath=file.path("plots", "mesh"),
    mar=c(0,0,0,0), mgp=c(0,0,0), oma=c(0,0,0,0))
#png(file.path(plot_dir, 'inputs.png'), res=400, 
#    width=150*plot_width, height=150*plot_width, 
#    family=font)
#par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0), oma=c(0, 0, 0, 0))
#H <- 1-combos$Y-combos$Z
#varmax <- max(H)
#varmin <- min(H)
#cols <- col.pal[floor((H-varmin) / (varmax-varmin) * (length(col.pal)-1))+1]
#splom(combos, cex=0.001, pch=3,
#      col=cols,
#      #col.pal[floor(ages/max(ages)*length(col.pal))],
#      xlab=NULL, ylab=NULL, 
#      axis.text.cex=0.25,
#      axis.text.lineheight=0.0001,
#      axis.line.tck=0.25,
#      xaxs='n', yaxs='n',
#      varname.cex=0.5,
#      varnames=c(expression(M[0]), expression(Y[0]), expression(Z[0]), 
#               expression(alpha["MLT"])))
#dev.off()

#png(file.path(plot_dir, 'inputs-legend.png'), res=400, 
#    width=150*plot_width/8, height=150*plot_width, 
#    family=font)
#par(mar=c(0, 0, 0, 0), mgp=c(0, 0, 0), oma=c(0, 0, 0, 0))
#color.legend(par()$usr[2], par()$usr[1], par()$usr[4], par()$usr[3], 
#             signif(quantile(seq(varmin, varmax, length=1000), 
#                    c(0.05, 0.275, 0.5, 0.725, 0.95)), 3), 
#             col.pal[1:length(col.pal)], gradient='y', align='rb')
#mtext(expression(H_0), 4, line=4.5, cex=1.3)
#dev.off()


# HR scatter
Z_name <- 'M'
varmax <- round(max(seis.DF[[Z_name]]), 2)
varmin <- round(min(seis.DF[[Z_name]]), 2)
png(file.path(plot_dir, 'HR-M.png'), 
    family=font, res=400, width=plot_width*250, height=plot_height*250)
par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
for (simulation_i in 1:nrow(combos)) {
    DF <- merge(seis.DF, combos[simulation_i,])
    relation <- log10(DF$L) ~ DF$Teff
    color <- col.pal[
        floor((DF[[Z_name]]-varmin)/(varmax-varmin)*length(col.pal))+1]
    cex <- 0.01
    if (simulation_i == 1) {
        plot(relation, 
            pch=1, axes=FALSE,
            col=color, cex=cex, tcl=0,
            ylim=range(log10(seis.DF$L)),#, 0, 1),
            xlim=rev(range(seis.DF$Teff)), 
            xlab=expression(T["eff"]/K), 
            ylab=expression(L / L['\u0298']))
        abline(v=5777, lty=3, col='black')
        abline(h=0, lty=3, col='black')
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                majorn=c(4, 3, 4, 3),
                unlog='y', mgp=c(2, 0.25, 0))
    } else {
        points(relation, col=color, pch=20, cex=cex)
    }
}
points(5777, 0, pch=1, cex=1)
points(5777, 0, pch=20, cex=0.1)
var1range <- diff(par()$usr)[1]
color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
             par()$usr[2]+0.10*var1range, par()$usr[4], 
    signif(quantile(seq(varmin, varmax, length=1000), 
        c(0, 0.25, 0.5, 0.75, 1)), 3), 
    col.pal[1:length(col.pal)], gradient='y', align='rb')
mtext(expression(M/M['\u0298']), 4, line=4.5, cex=1.3)
dev.off()




# HR mesh
mesh <- with(seis.DF, interp(Teff, log10(L), M,
    xo=seq(min(Teff), max(Teff), length=40),
    yo=seq(log10(min(L)), log10(max(L)), length=40)))
    
spacing <- c(0.001, seq(2, 98, 2)/100, 0.999)
mesh <- interp(log10(seis.DF[['Teff']]), log10(seis.DF[['L']]), 
    seis.DF[['M']],
    xo=seq(quantile(log10(seis.DF[['Teff']]), 0.001),
           quantile(log10(seis.DF[['Teff']]), 0.999),
           length=50),
    yo=seq(quantile(log10(seis.DF[['L']]), 0.001),
           quantile(log10(seis.DF[['L']]), 0.999),
           length=50))

cairo_pdf(file.path(plot_dir, 'mesh-HR-M.pdf'), 
    width=plot_width, height=plot_height, family=font)
par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
filled.contour(mesh,
    ylim=log10(range(seis.DF$L)),
    xlim=rev(log10(range(seis.DF$Teff))), 
    xaxs='r', yaxs='r',
    levels=color_levels[['M']], 
    color=colorRampPalette(brewer.pal(11, "Spectral")),
    key.axes={
        axis(4, cex.axis=1.5, tcl=0, line=0)
        mtext(expression(M/M['\u0298']), side=4, las=3, line=4, cex=2)
    },
    plot.axes={
        contour(mesh, add=TRUE, labcex=1, levels=Z_levels[['M']],
           )#method="simple")
        points(log10(5777), 0, pch=1, cex=1)
        points(log10(5777), 0, pch=20, cex=0.1)
        abline(v=log10(5777), lty=3, col=adjustcolor('black', alpha.f=0.25))
        abline(h=0, lty=3, col=adjustcolor('black', alpha.f=0.25))
        magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                majorn=c(4, 3, 4, 3), cex.axis=1.5,
                unlog='xy')#, mgp=c(2, 0.25, 0))
    },
    plot.title={
        title(xlab=expression(T["eff"]/K), cex.lab=2, line=3)
        title(ylab=expression(L / L['\u0298']), cex.lab=2, line=3)
    })
dev.off()




## Color by age, mass, Y0, X(He), metallicity, mix length, and core hydrogen
scatter_mesh <- function(plotname, var1, var2, var3, label1, label2, label3,
    seis.DF) { 
    # scatter 
    varmax <- max(seis.DF[[var3]])
    varmin <- min(seis.DF[[var3]])
    png(file.path(hrcdr_dir, paste0(plotname, '-', var3, '.png')), 
        family=font, res=400, width=plot_width*250, height=plot_height*250)
    par(mar=c(3, 4, 1, 6), mgp=c(2, 0.25, 0), cex.lab=1.3)
    for (simulation_i in 1:nrow(combos)) {
        DF <- merge(seis.DF, combos[simulation_i,])
        relation <- DF[[var2]] ~ DF[[var1]]
        color <- col.pal[floor((DF[[var3]]-varmin)/(varmax-varmin)*
            length(col.pal))+1]
        cex <- 0.01
        if (simulation_i == 1) {
            plot(relation, 
                pch=20, axes=FALSE,
                col=color, cex=cex, tcl=0,
                ylim=quantile(seis.DF[[var2]], c(0.001, 0.999)), 
                xlim=quantile(seis.DF[[var1]], c(0.001, 0.999)),
                xlab=label1, ylab=label2)
            abline(v=solar_vals[[var1]], lty=3, col='black')
            abline(h=solar_vals[[var2]], lty=3, col='black')
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.25, 0))
        } else {
            points(relation, col=color, pch=20, cex=cex)
        }
    }
    points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
    points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
    var1max <- max(seis.DF[[var1]])
    var1range <- diff(par()$usr)[1]
    color.legend(par()$usr[2]+0.05*var1range, par()$usr[3], 
                 par()$usr[2]+0.10*var1range, par()$usr[4], 
        signif(quantile(seq(varmin, varmax, length=1000), 
            c(0, 0.25, 0.5, 0.75, 1)), 2), 
        col.pal[1:length(col.pal)], gradient='y', align='rb')
    mtext(label3, 4, line=4.5, cex=1.3)
    dev.off()
    
       
    # mesh
    spacing <- c(0.001, (1:99)/100, 0.999)
    spacing <- c(0.001, seq(2, 98, 2)/100, 0.999)
    mesh <- interp(seis.DF[[var1]], seis.DF[[var2]], seis.DF[[var3]],
        xo=quantile(seis.DF[[var1]], spacing),
        yo=quantile(seis.DF[[var2]], spacing))
        #xo=seq(min(seis.DF[[var1]]), max(seis.DF[[var1]]), length=100),
        #yo=seq(min(seis.DF[[var2]]), max(seis.DF[[var2]]), length=100))
    cairo_pdf(file.path(hrcdr_dir, 
            paste0('mesh-', plotname, '-', var3, '.pdf')), 
        width=plot_width, height=plot_height, family=font)
    par(mar=c(5, 6, 1, 0), mgp=c(2, 0.25, 0), cex.lab=1)
    filled.contour(mesh,
        levels=color_levels[[var3]], 
        color=colorRampPalette(brewer.pal(11, "Spectral")),
        key.axes={
            axis(4, cex.axis=1.5, tcl=0, line=0)
            mtext(label3, side=4, las=3, line=4, cex=2)
        },
        plot.axes={
            contour(mesh, add=TRUE, labcex=0.5, levels=Z_levels[[var3]])
            if (any(grepl(var1, names(solar_vals))) && 
                any(grepl(var2, names(solar_vals)))) {
                points(solar_vals[[var1]], solar_vals[[var2]], pch=1, cex=1)
                points(solar_vals[[var1]], solar_vals[[var2]], pch=20, cex=0.1)
            }
            magaxis(side=1:4, family=font, tcl=0.25, labels=c(1,1,0,0),
                    mgp=c(2, 0.5, 0), cex.axis=1.5)
        },
        plot.title={
            title(xlab=label1, cex.lab=2, line=3)
            title(ylab=label2, cex.lab=2, line=3)
        })
    dev.off()
}

for (Z_name in Z_names) {
    scatter_mesh('JCD', 'Dnu0_median', 'dnu02_median', Z_name, 
        expression(Delta*nu/mu*Hz), 
        expression(delta*nu[0*","*2]/mu*Hz), 
        Z_labels[which(Z_name==Z_names)], seis.DF)
    
    #scatter_mesh('HR', 'Teff', 'L', Z_name, 
    #    expression(T[eff]~"["*K*"]"), 
    #    expression(L/L['\u0298']), 
    #    Z_labels[which(Z_name==Z_names)])
}

pca <- prcomp(seis.DF[,-1:-8, with=0], center=TRUE, scale.=TRUE)
pcs <- pca$x[,cumsum(pca$sdev)/sum(pca$sdev)<0.99]

scatter_mesh('PC23', 'PC2', 'PC3', 'age', 
    expression(PC[2]), 
    expression(PC[3]), 
    Z_labels[which('age'==Z_names)], cbind(seis.DF, pcs))

scatter_mesh('PC12', 'PC1', 'PC2', 'Hc', 
    expression(PC[1]), 
    expression(PC[2]), 
    Z_labels[which('Hc'==Z_names)], cbind(seis.DF, pcs))

scatter_mesh('PC14', 'PC1', 'PC4', 'radius', 
    expression(PC[1]), 
    expression(PC[4]), 
    Z_labels[which('radius'==Z_names)], cbind(seis.DF, pcs))

