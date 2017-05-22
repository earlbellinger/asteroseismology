#### Obtain model properties from evolutionary tracks 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

# parse reference model 
source('../scripts/utils.R') 
source('../scripts/seismology.R') 

wd <- getwd()
setwd('../inversion')
source('kernels.R')
source('models.R')
source('frequencies.R')

target.name <- 'modmix'
ref.name <- 'no_diffusion'
mode.set <- 'CygA'
k.pair <- u_Y
freqs <- get_freqs(target.name=target.name, mode.set=mode.set, perturb=T) 
m1 <- get_model(freqs=freqs, model.name=ref.name, target.name=target.name, 
                k.pair=k.pair, square.Ks=F) 

setwd(wd)
parallelStartMulticore(16)

n.l <- freqs[,1:2]
nl.names <- paste0('l', n.l$l, '_n', n.l$n)
xs <- seq(0.01, 1, 0.01) #seq(0.05, 0.35, 0.01)
m1.f1 <- m1$f1.spl(xs)
m1.f2 <- m1$f2.spl(xs)

# asdf <- (m1.f1 - asdf2)/m1.f1
#d.m1.f1 <- splinefun(m1$r, m1$d.f1.true)(xs) 
#plot(m1$r, m1$d.f1.true, xlim=range(xs), ylim=range(m1$d.f1.true, asdf+std, asdf-std), type='l')
#points(xs, asdf); points(xs, asdf+std); points(xs, asdf-std); dev.off()

calc_seps <- function(nus) {
    c(get_separations('r_sep', nus, 0)$separations,
      get_separations('r_sep', nus, 1)$separations,
      get_separations('r_avg', nus, 0)$separations,
      get_separations('r_avg', nus, 1)$separations)#,
      #get_separations('Dnu', nus, 0)$separations,
      #get_separations('Dnu', nus, 1)$separations,
      #get_separations('Dnu', nus, 2)$separations,
      #get_separations('Dnu', nus, 3)$separations,
      #get_separations('dnu', nus, 0)$separations,
      #get_separations('dnu', nus, 1)$separations)
}

m1.freqs <- m1$nus[,1:3]
names(m1.freqs)[3] <- 'nu'
m1.r <- calc_seps(m1.freqs) 

rem.surf.term <- function(nu.x, nu.y, E, nu_ac, dnu=1) { 
    
    #nu <- nu.x / nu_ac 
    #Xpinv <- ginv( matrix(c(nu**-2/E, nu**2/E, nu) / dnu, ncol=3) ) 
    #a.surf <- Xpinv %*% ( (nu.x-nu.y)/(nu.y*dnu) ) 
    #F_surf <- ( a.surf[[1]]*nu**-2/E + 
    #            a.surf[[2]]*nu**2/E + 
    #            a.surf[[3]]*nu ) 
    #(nu.x-nu.y)/nu.x - F_surf 
    
    nu <- nu.x / nu_ac 
    Xpinv <- ginv( matrix(c(nu**-1/E, nu**3/E, nu) / dnu, ncol=3) ) 
    a.surf <- Xpinv %*% ( (nu.x-nu.y)/dnu ) 
    F_surf <- ( a.surf[[1]]*nu**-1/E + 
                a.surf[[2]]*nu**3/E  ) 
    
    diffs <- (nu.x-nu.y-F_surf)/nu.x
    
    diffs #- weighted.mean(diffs, 1/dnu)
    
} 

test.mat <- do.call(plyr:::rbind.fill, parallelMap(function(ii) {
    nu.y <- rnorm(nrow(freqs), freqs$nu, freqs$dnu)
    freqs2 <- freqs
    freqs2$nu <- nu.y
    
    ratios <- calc_seps(freqs2)
    result <- ratios #(m1.r - ratios) / m1.r
    names(result) <- c(paste0('r_', 1:length(result)))
    
    #result <- rem.surf.term(m1$nus$nu.x, nu.y, m1$nus$E, m1$nu_ac, freqs$dnu)
    #names(result) <- paste0('r_', 1:length(result)) #nl.names
    data.frame(t(result))
}, ii=1:1000))

write.table(test.mat, 'test_mat_abs.dat', quote=F, row.names=F, sep='\t')
#write.table(test.mat, 'test_mat.dat', quote=F, row.names=F, sep='\t')



# parse target frequencies 
#ptm <- proc.time()
simulations.dir <- 'simulations'
training.mat <- do.call(plyr:::rbind.fill, parallelMap(function(simulation) {
    do.call(plyr:::rbind.fill, Map(function(logs.dir) {
        directory <- file.path(simulations.dir, simulation, logs.dir)
        print(directory)
        log.files <- list.files(directory)
        freq.fnames <- log.files[grepl('-freqs.dat', log.files)]
        do.call(plyr:::rbind.fill, parallelMap(function(freq.fname) {
            
            #print(freq.fname)
            
            nus <- tryCatch(
                read.table(file.path(directory, freq.fname), 
                    col.names=c('l', 'n', 'nu', 'E')),
                error=function(e) NULL)
            if (is.null(nus)) return(NULL)
            
            nus <- merge(nus, n.l, by=c('n', 'l')) 
            if (nrow(nus) != nrow(n.l)) return(NULL) 
            
            fgong.fname <- sub('-freqs.dat', '.data.FGONG.dat', freq.fname) 
            if (!file.exists(file.path(directory, fgong.fname))) return(NULL)
            fgong <- tryCatch(
                read.table(file.path(directory, fgong.fname), header=1),
                error=function(e) NULL)
            if (is.null(fgong)) return(NULL)
            
            m2.f1 <- splinefun(fgong$x, fgong[[k.pair$f1]])(xs) 
            m2.f2 <- splinefun(fgong$x, fgong[[k.pair$f2]])(xs) 
            
            df1.fx <- m2.f1 #(m1.f1 - m2.f1) / m1.f1 #
            df2.fx <- m2.f2 #(m1.f2 - m2.f2) / m1.f2 #
            
            #d.nu <- rem.surf.term(m1$nus$nu.x, nus$nu, m1$nus$E, m1$nu_ac, 
            #    freqs$dnu)
                        
            ratios <- calc_seps(nus) 
            if (length(ratios) != length(m1.r)) return(NULL)
            d.nu <- ratios #(m1.r - ratios) / m1.r
            
            result <- c(d.nu, df1.fx, df2.fx)
            names(result) <- c(paste0('r_', 1:length(d.nu)),
            #names(result) <- c(nl.names, 
                               paste0('df1_', sub('\\.', '', paste0(xs))),
                               paste0('df2_', sub('\\.', '', paste0(xs))))
            
            data.frame(t(result))
        }, freq.fname=freq.fnames))
    }, logs.dir=c('LOGS_MS')))
}, simulation=list.files(simulations.dir)))
#proc.time() - ptm

write.table(training.mat, 'training_mat_abs.dat', quote=F, row.names=F, sep='\t')
#write.table(training.mat, 'training_mat.dat', quote=F, row.names=F, sep='\t')


