library(magicaxis)
library(Bolstad) # for using Simpson's rule 

Rsun <- 69598946770.0
col.pal <- colorRampPalette(c("#332288", "#88CCEE", "#44AA99", "#117733", 
    "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"))(4)
freqs <- read.table(file.path('modelS', 'modelS', 'ModelS-freqs.dat'), 
    col.names=c('l', 'n', 'nu'))

ikit_dir <- file.path('compare', 'ik')

kern_dir <- 'modelS_ker'
kernels <- list.files(kern_dir)
kernels <- kernels[grepl('.dat', kernels)]

for (pair in c('c2-rho', 'Gamma1-rho')) {
    pair_names <- strsplit(pair, '-')[[1]]
    f1 <- pair_names[1]
    f2 <- pair_names[2]
    
    diffs <- data.frame()
    
    kpair_fnames <- kernels[grepl(pair, kernels)]
    for (kpair_fname in kpair_fnames) {
        kern <- read.table(file.path(kern_dir, kpair_fname),
            col.names=c('x', 'Kf1', 'Kf2'))
        
        parts <- strsplit(kpair_fname, '_')
        ell <- strtoi(strsplit(parts[[1]][2], '=')[[1]][2])
        nnn <- strtoi(strsplit(strsplit(parts[[1]][3], '=')[[1]][2], 
            '.dat')[[1]][1])
        
        nu <- freqs[freqs$l == ell & freqs$n == nnn,]$nu
        
        ikit <- merge(
            read.table(file.path(ikit_dir, pair,
                                 paste0(ell, '_', nnn)),
                       col.names=c('x', 'Kf1')),
            read.table(file.path(ikit_dir, paste0(f2, '-', f1),
                                 paste0(ell, '_', nnn)),
                       col.names=c('x', 'Kf2')))
        
        
        # make a fake Gaussian perturbation and calculate the change in freq
        x <- kern$x
        
        df_f <- exp( -(x - 0.3)**2/0.3**2 )
        df_f <- (df_f - min(df_f))/10
        
        # calculate dw/w = int K_f * df_f dr
        dw_w_ikit = integrate( splinefun(x, ikit$Kf1 * df_f), min(x), max(x) )$value
        dw_w_kern = integrate( splinefun(x, kern$Kf1 * Rsun * df_f), min(x), max(x) )$value
        
        abs_diff_f1 <- abs(dw_w_ikit - dw_w_kern)
        rel_diff_f1 <- abs((dw_w_ikit - dw_w_kern) / dw_w_ikit * 100)
        
        dw_w_ikit = integrate( splinefun(x, ikit$Kf2 * df_f), min(x), max(x) )$value
        dw_w_kern = integrate( splinefun(x, kern$Kf2 * Rsun * df_f), min(x), max(x) )$value
        
        abs_diff_f2 <- abs(dw_w_ikit - dw_w_kern)
        rel_diff_f2 <- abs((dw_w_ikit - dw_w_kern) / dw_w_ikit * 100)
        
        diffs <- rbind(diffs, 
            data.frame(l=ell, n=nnn, nu=nu, 
                abs_diff_f1=abs_diff_f1, rel_diff_f1=rel_diff_f1,
                abs_diff_f2=abs_diff_f2, rel_diff_f2=rel_diff_f2,
                max_diff_f1=max(abs(ikit$Kf1 - kern$Kf1*Rsun)),
                max_diff_f2=max(abs(ikit$Kf2 - kern$Kf2*Rsun))))
    }
    
    
    
    f1_exp <- if (f1 == 'c2') bquote(c^2) else if (f1 == 'Gamma1') bquote(Gamma[1])
    f2_exp <- if (f2 == 'rho') bquote(rho)
    
    cairo_pdf(paste0('abs_diff-', pair, '.pdf'), width=5, height=5)
    par(mar=c(4,6,1,1))
    plot(NA, axes=F, tcl=0, log='y',
         ylab=as.expression(bquote(abs(
             integral(K[.(f1_exp)*","*.(f2_exp)]^"ikit" * 
                      frac(d*.(f1_exp), .(f1_exp)) * dr) - 
             integral(K[.(f1_exp)*","*.(f2_exp)]^"me" * 
                      frac(d*.(f1_exp), .(f1_exp)) * dr))%.%nu/mu*Hz)),
         xlab=expression("Frequency"~nu/mu*Hz),
         ylim=range(diffs$nu*diffs$abs_diff_f1),
         xlim=range(diffs$nu))
    for (ell in unique(diffs$l)) {
        ells <- diffs[diffs$l == ell,]
        ells <- ells[order(ells$nu),]
        with(ells, points(nu, nu*abs_diff_f1, pch=ell+1, col=col.pal[ell+1]))
        with(ells, lines(nu, nu*abs_diff_f1, lty=ell+1, col=col.pal[ell+1]))
    }
    legend('bottomright', legend=unique(diffs$l), 
        lty=unique(diffs$l)+1, 
        col=col.pal[unique(diffs$l)+1], 
        pch=unique(diffs$l)+1)
    magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25)
    dev.off()
    
    
    
    cairo_pdf(paste0('abs_diff-', f2, '-', f1, '.pdf'), width=5, height=5)
    par(mar=c(4,6,1,1))
    plot(NA, axes=F, tcl=0, log='y',
         ylab=as.expression(bquote(abs(
             integral(K[.(f2_exp)*","*.(f1_exp)]^"ikit" * 
                      frac(d*.(f2_exp), .(f2_exp)) * dr) - 
             integral(K[.(f2_exp)*","*.(f1_exp)]^"me" * 
                      frac(d*.(f2_exp), .(f2_exp)) * dr))%.%nu/mu*Hz)),
         xlab=expression("Frequency"~nu/mu*Hz),
         ylim=range(diffs$nu*diffs$abs_diff_f2),
         xlim=range(diffs$nu))
    for (ell in unique(diffs$l)) {
        ells <- diffs[diffs$l == ell,]
        ells <- ells[order(ells$nu),]
        with(ells, points(nu, nu*abs_diff_f2, pch=ell+1, col=col.pal[ell+1]))
        with(ells, lines(nu, nu*abs_diff_f2, lty=ell+1, col=col.pal[ell+1]))
    }
    legend('bottomright', legend=unique(diffs$l), 
        lty=unique(diffs$l)+1, 
        col=col.pal[unique(diffs$l)+1], 
        pch=unique(diffs$l)+1)
    magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25)
    dev.off()
    
    
    
    
}









cairo_pdf(paste0('frac_diff.pdf'), width=5, height=5)
plot(ikit$x, ikit$Kf1 - kern$Kf1*Rsun, 
    ylim=c(-0.007, 0.004), 
    axes=F, tcl=0,
    xlab=expression("Fractional radius x"),
    ylab=expression("Difference" ~ K["ikit"] - K["me"]),
    type='l')
magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25)
lines(ikit$x, ikit$Kf2 - kern$Kf2*Rsun, col='red', lty=2)
legend('bottomleft', lty=c(1,2), col=c('black', 'red'),
    legend=c(expression(K[c^2*","*rho]^"ikit" - K[c^2*","*rho]^"me"), 
             expression(K[rho*","*c^2]^"ikit" - K[rho*","*c^2]^"me")))
dev.off()





















Rsun <- 69598946770.0 # Model S
      # 69598946770.0 # Centre
      # 69599062580.0 # 5b

basu <- merge(
    read.table('basu/kercrho/eig_l=0_n=21.dat',
        col.names=c("x", "Kc2rho")), 
    read.table('basu/kerrhoc/eig_l=0_n=21.dat',
        col.names=c("x", "Krhoc2")))

kern <- read.table('modelS_ker/c2-rho_l=1_n=21.dat',
    col.names=c("x", "Kc2rho", "Krhoc2"))

kern <- read.table('modelS_ker/Gamma1-rho_l=1_n=21.dat',
    col.names=c("x", "Kgamma1rho", "Krhogamma1"))

ikit <- merge(merge(
    read.table('compare/ik-c2rho_1.dat',
        col.names=c("x", "Kc2rho")), 
    read.table('compare/ik-rhoc2_1.dat',
        col.names=c("x", "Krhoc2"))),
    read.table('compare/ik-rhogamma1_1.dat',
        col.names=c("x", "Krhogamma1")))

plot(ikit$x, ikit$Krhoc2 - kern$Krhoc2*Rsun, 
    ylim=c(-0.002, 0.002), 
    axes=F, tcl=0,
    xlab=expression("Fractional radius x"),
    ylab=expression("Difference" ~ K["ikit"] - K["me"]),
    type='l')
magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25)
lines(ikit$x, ikit$Kc2rho - kern$Kc2rho*Rsun, col='red', lty=2)
#legend('bottomleft', lty=c(1,2), col=c('black', 'red'),
#    legend=c("K_c2_rho", "K_rho_c2"))
dev.off()



basu_Kc2rho <- with(basu, splinefun(x, Kc2rho))
kern_Kc2rho <- with(kern, splinefun(x, Kc2rho * Rsun))

basu_Krhoc2 <- with(basu, splinefun(x, Krhoc2))
kern_Krhoc2 <- with(kern, splinefun(x, Krhoc2 * Rsun))

x <- basu$x

plot(basu$x, basu$Kc2rho - kern_Kc2rho(basu$x), 
    ylim=c(-0.01, 0.03),
    axes=F, tcl=0,
    xlab=expression("Fractional radius x"),
    ylab=expression("Difference" ~ K["Basu"] - K["Me"]),
    type='l')
magaxis(1:4, labels=c(1,1,0,0), tcl=-0.25)
lines(basu$x, basu$Krhoc2 - kern_Krhoc2(basu$x), col='red', lty=2)
#legend('bottomleft', lty=c(1,2), col=c('black', 'red'),
#    legend=c("K_c2_rho", "K_rho_c2"))
dev.off()



# test c^2 kernel 
x <- ikit$x

dc_c <- exp( -(x - 0.3)**2/0.3**2 )
dc_c <- (dc_c - min(dc_c))/10

# calculate dw/w = int K_c2_rho * dc_c dr
dw_w_ikit = integrate( splinefun(x, ikit$Kc2rho * dc_c), min(x), max(x) )$value
dw_w_kern = integrate( splinefun(x, kern$Kc2rho * Rsun * dc_c), min(x), max(x) )$value
print(abs(dw_w_ikit - dw_w_kern))
print(paste0(abs((dw_w_ikit - dw_w_kern) / dw_w_ikit * 100), '%'))

dw_w_ikit = integrate( splinefun(x, ikit$Krhoc2 * dc_c), min(x), max(x) )$value
dw_w_kern = integrate( splinefun(x, kern$Krhoc2 * Rsun * dc_c), min(x), max(x) )$value
print(abs(dw_w_ikit - dw_w_kern))
print(paste0(abs((dw_w_ikit - dw_w_kern) / dw_w_ikit * 100), '%'))



# make a small gaussian 
dc_c <- exp( -(basu$x - 0.3)**2 )
dc_c <- (dc_c - min(dc_c))/10
# calculate dw/w = int K_c2_rho * dc_c dr
dw_w_basu = integrate( splinefun(x, basu$Kc2rho * dc_c), min(x), max(x) )$value

dc_c <- exp( -(kern$x - 0.3)**2 )
dc_c <- (dc_c - min(dc_c))/10
dw_w_kern = integrate( splinefun(kern$x * Rsun, kern$Kc2rho * dc_c), min(x), Rsun )$value


# test rho kernel
dc_c <- exp( -(basu$x - 0.5)**2 )
dc_c <- (dc_c - min(dc_c))/10
# calculate dw/w = int K_c2_rho * dc_c dr
dw_w_basu = integrate( splinefun(x, basu$Krhoc2 * dc_c), min(x), max(x) )$value


dc_c <- exp( -(kern$x - 0.5)**2 )
dc_c <- (dc_c - min(dc_c))/10
dw_w_kern = integrate( splinefun(kern$x, Rsun * kern$Krhoc2 * dc_c), min(x), max(x) )$value




