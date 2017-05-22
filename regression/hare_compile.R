#### Make a table of the fundamental and chemical properties for SpaceInn H-n-H
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar predictions & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

fundamental <- read.table(header=0, 
    col.names=c("Model", "M", "radius", "rho", "log_g", "age", 
                "r_bcz", "tau", "tau_bcz", "alpha", "overshoot", "diffusion"), 
    textConnection("
Coco 0.78 0.815 2.029 4.508 9.616 0.746 2999 1293 1.6708 0 0
Aardvark 1.00 0.959 1.596 4.474 3.058 0.731 3371 1405 1.6708 0 0
Elvis 1.00 1.087 1.097 4.365 6.841 0.727 4045 1656 1.6708 0 0
Henry 1.10 1.138 1.051 4.367 2.055 0.839 4199 2280 1.8045 0 1
Izzy 1.10 1.141 1.041 4.364 2.113 0.849 4219 2353 1.6708 0 0
Blofeld 1.22 1.359 0.684 4.257 2.595 0.838 5107 2811 1.6738 0.10 1
Diva 1.22 1.353 0.693 4.261 4.622 0.767 5082 2288 1.6708 0.10 0
Felix 1.33 1.719 0.369 4.091 2.921 0.842 7014 3830 1.6708 0.05 0
George 1.33 1.697 0.383 4.102 2.944 0.875 6930 4156 1.6708 0.25 0
Jam 1.33 1.468 0.592 4.228 1.681 0.905 5666 3718 1.6708 0.05 0"))

chemical <- read.table(header=0, 
    col.names=c("Model", "X", "Z", "X_surf", "Z_surf", 
                "z_div_x_surf", "z_div_x_surf_solar"), 
    textConnection("
Aardvark 0.71550 0.01755 0.71543 0.01755 0.02453 0.0245
Blofeld 0.71400 0.02000 0.78280 0.01579 0.02018 0.0165
Coco 0.74140 0.00360 0.74132 0.00360 0.00486 0.0245
Diva 0.72600 0.02600 0.72593 0.02600 0.03582 0.0245
Elvis 0.71550 0.01755 0.71543 0.01755 0.02453 0.0245
Felix 0.71550 0.01755 0.71543 0.01755 0.02453 0.0245
George 0.71550 0.01755 0.71543 0.01755 0.02453 0.0245
Henry 0.72600 0.01000 0.78010 0.00825 0.01058 0.0245
Izzy 0.72600 0.01000 0.72593 0.01000 0.01378 0.0245
Jam 0.71550 0.01755 0.71543 0.01755 0.02453 0.0245"))

DF <- merge(fundamental, chemical)

output <- DF[c('Model', 'M', 'Z', 'age', 'R', 'log_g', 'alpha', 'ov', 
    'diffusion')]
output <- cbind(output, data.frame(
    Y = 1-DF$X-DF$Z,
    Y_surf = 1 - DF$X_surf - DF$Z_surf))

write.table(output, file.path('data', 'hares.dat'), row.names=F, quote=F,
    sep='\t')

