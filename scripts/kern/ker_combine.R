#### Collate extracted kernel ascii files into one big file 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

args <- commandArgs(TRUE)
kpair <- if (length(args)>0) as.numeric(args[1]) else 5 # default to u-Y kernels

ker.type <- ifelse(kpair < 10, "E_K", "L_K")
ker.names <- list(
  '1'  =  c("c2", "rho"),
  '2'  =  c("Gamma1", "rho"),
  '3'  =  c("c2", "Gamma1"),
  '4'  =  c("u", "Gamma1"),
  '5'  =  c("u", "Y"),
  '6'  =  c("rho", "Y"),
  '7'  =  c("c", "Gamma1_over_c"),
  '11' =  c("c2", "rho"),
  '12' =  c("Gamma1", "rho"),
  '13' =  c("c2", "Gamma1"),
  '14' =  c("u", "Gamma1"),
  '17' =  c("c", "Gamma1_over_c")
)

ker.name <- ker.names[kpair][[1]]

for (pair in 1:2) {
    ker <- NULL
    directory <- paste0("ichkr-", kpair, "_", pair)
    files <- list.files(directory)
    files <- files[grep('.dat', files)]
    for (file in files) {
        DF <- read.table(file.path(directory, file), 
            col.names=c('x', substr(file, 1, nchar(file)-4)))
        ker <- if (is.null(ker)) DF else merge(ker, DF, by='x')
    }
    
    # sort by l and n
    cols <- names(ker)[-1]
    ell <- as.numeric(Map(function(x) x[[3]], strsplit(cols, '[l|n|.|_]')))
    n <- as.numeric(Map(function(x) x[[6]], strsplit(cols, '[l|n|.|_]')))
    i <- data.frame(ell, n)
    ker <- ker[,c(1, 1+order(i[,1], i[,2]))]
    
    write.table(ker, quote=F, sep='\t', row.names=F,
        file=paste0(ker.type, "_", ker.name[pair], '-', 
                    ker.name[ ifelse(pair==1, 2, 1) ],
                    '.dat'))
}



