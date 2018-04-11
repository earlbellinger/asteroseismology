#### Makes models corresponding to the mode of the RF posterior distribution 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

options(scipen=99999999)

phys <- c('inv', 'inv-ov', 'inv-D', 'inv-Dov')
ov <- c(0, 0.2, 0, 0.2)
D <- c(0, 0, 1, 1)

for (phy_i in 1:length(phys)) {
    phy <- phys[phy_i]
    
    #dir.create(file.path('models', phy), showWarnings=F)
    
    DF <- read.table(file.path('..', 'regression', 
            paste0('learn-', phy), paste0('tables-', phy), 
            'inversions_modes.dat'), 
        header=1)
    
    DF <- DF[unique(DF$Name),]
    DF <- DF[!is.na(strtoi(DF$Name)),]
    DF <- DF[order(strtoi(DF$Name)),]
    
    command <- paste0("maybe_sub.sh -e -p 1 ./dispatch.sh -f",
        " -d ref_mods/", phy, 
        " -n ", DF$Name, 
        " -M ", DF$M, 
        " -Y ", DF$Y, 
        " -Z ", DF$Z, 
        " -a ", DF$alpha, 
        " -c ", DF$age, 'd9', 
        " -D ", D[phy_i], 
        " -g ", D[phy_i], 
        " -o ", ov[phy_i], 
        ifelse(D[phy_i], " -t", ""))
    
    print(command)
    for (ii in command) system(ii)
}
