library(magicaxis)
library(Bolstad) # for using Simpson's rule 

Rsun <- 69598946770.0 # Model S
      # 69598946770.0 # Centre
      # 
      # 69599062580.0 # 5b

basu <- merge(
    read.table('basu/kercrho/eig_l=0_n=21.dat',
        col.names=c("x", "Kc2rho")), 
    read.table('basu/kerrhoc/eig_l=0_n=21.dat',
        col.names=c("x", "Krhoc2")))
kern <- read.table('modelS_ker/c2-rho_l=0_n=21.dat',
    col.names=c("x", "Kc2rho", "Krhoc2"))

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
legend('bottomleft', lty=c(1,2), col=c('black', 'red'),
    legend=c("K_c2_rho", "K_rho_c2"))
dev.off()



# test c^2 kernel 
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

