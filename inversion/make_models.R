#### Takes a covariance file and makes 9 models fitted to span their M and R 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)
options(scipen=99999)
library(randomForest)

### LIBRARIES 
args <- commandArgs(TRUE)
star <- if (length(args)>0) args[1] else '8006161'
covs.directory <- if (length(args)>1) args[2] else file.path('learn-inversions', 
    'covs-simulations', 'inversions')
obs.directory <- if (length(args)>2) args[3] else 'inversions'

DF <- read.table(file.path('..', 'regression', covs.directory, 
    paste0(star, '.dat')), header=1)

params <- read.table(file.path('..', 'regression', 'data', obs.directory, 
    paste0(star, '-obs.dat')), header=1)

M.mean <- mean(DF$M)
M.sd <- sd(DF$M)

if ('radius' %in% params$name) {
    R.mean <- params$value[params$name == 'radius']
    R.sd <- params$uncertainty[params$name == 'radius']
    print('radius in params')
    #DF <- cbind(DF, read.table('../regression/perturb/
} else {
    R.mean <- mean(DF$radius)
    R.sd <- sd(DF$radius)
    
    X <- data.frame(M=DF$M, R=DF$radius)
    mu <- colMeans(X)
    
    # uncorrelate with PCA 
    pca <- prcomp(X, center=T)#, scale=T)
    
    # calculate 1-sigma distance 
    D2 <- mahalanobis(pca$x, rep(0, times=length(mu)), cov=cov(pca$x))
    one.sigma <- pca$x[D2<=1,]
    
    # calculate oval size 
    r1 <- max(abs(range(one.sigma[,1])))
    r2 <- max(abs(range(one.sigma[,2])))
    
    # go around the oval picking out points and projecting back to orig space 
    new.points <- data.frame(t(sapply(0:7, function(ii)
        scale(c(r1*cos(pi*2*ii/8), r2*sin(pi*2*ii/8)) %*% t(pca$rotation), 
            center=-mu, scale=F))))
    
    # order the data
    new.points <- new.points[order(new.points[,1]),]
    new.points[1:3,] <- new.points[1:3,][order(new.points[1:3,2]),]
    new.points[4:5,] <- new.points[4:5,][order(new.points[4:5,2]),]
    new.points[6:8,] <- new.points[6:8,][order(new.points[6:8,2]),]
    new.points <- rbind(new.points[1:4,], c(M.mean, R.mean), new.points[5:8,])
}

age <- floor(mean(DF$age)*10**9)
luminosity <- if ('luminosity' %in% params$names) {
    params$value[params$name == 'luminosity']
} else mean(DF$L)
FeH <- params$value[params$name == 'Fe/H']

rf.Y     <- randomForest(Y ~ .,     as.data.frame(cbind(X, data.frame(Y=DF$Y,         age=DF$age, L=DF$L))))
rf.alpha <- randomForest(alpha ~ ., as.data.frame(cbind(X, data.frame(alpha=DF$alpha, age=DF$age, L=DF$L))))

lowhigh <- c('low', 'mean', 'high') 
for (M in 1:3) { 
    for (R in 1:3) { 
        mass    <- c(M.mean - M.sd, M.mean, M.mean + M.sd)[M] 
        radius  <- c(R.mean - R.sd, R.mean, R.mean + R.sd)[R] 
        Y <- mean(DF$Y)
        alpha <- mean(DF$alpha)
        if (exists('new.points')) {
            mass   <- new.points[(M-1)*3+R, 1]
            radius <- new.points[(M-1)*3+R, 2]
            Y     <- predict(rf.Y,     data.frame(M=mass, R=radius, age=age, L=luminosity))
            alpha <- predict(rf.alpha, data.frame(M=mass, R=radius, age=age, L=luminosity))
        }
        command <- paste("maybe_sub.sh -p 1 Rscript calibrate.R", 
            star, 
            paste0(lowhigh[R], 'R', lowhigh[M], 'M'), 
            mass, 
            log10(radius), 
            age, 
            log10(luminosity), 
            FeH, 
            Y, 
            mean(DF$Z), 
            alpha)
        print(command) 
        system(command) 
    } 
} 


if (F) {
X <- data.frame(M=DF$M, R=DF$radius)
mu <- colMeans(X)
pca <- prcomp(X, center=T)#, scale=T)
D2 <- mahalanobis(pca$x, c(0, 0), cov=cov(pca$x))
plot(X, #pca$x[,1], pca$x[,2], 
    pch=20, cex=0.5, 
    col=c('red', 'black', 'blue')[ifelse(D2<1, 1, ifelse(D2<2, 2, 3))])
one.sigma <- pca$x[D2<=1,]
r1 <- max(abs(range(one.sigma[,1])))
r2 <- max(abs(range(one.sigma[,2])))

new.points <- t(sapply(0:7, function(ii)
    scale(c(r1*cos(pi*2*ii/8), r2*sin(pi*2*ii/8)) %*% t(pca$rotation), 
        center=-mu, scale=F)))
new.points <- new.points[order(new.points[,1], new.points[,2]),]


for (ii in 0:7) {
    theta <- pi*2*ii/8
    X2 <- scale(c(r1*cos(theta), r2*sin(theta)) %*% t(pca$rotation), center=-mu, scale=F)
    #new.M <- r1*cos(theta)
    #new.R <- r2*sin(theta)
    points(X2, pch=20, cex=1, col='orange')
}
dev.off()
#X2 <- scale(pca$x %*% t(pca$rotation), center=-mu, scale=F)


#DF2 <- DF2[order(DF2$M),]
D2 <- mahalanobis(DF2, c(mean(DF2$M), mean(DF2$R)), cov=cov(DF2))
with(DF2, 
    plot(M, R, pch=20, cex=0.5, 
        col=c('red', 'black', 'blue')[ifelse(D2<1, 1, ifelse(D2<2, 2, 3))])
)

pca <- prcomp(DF2, center=T, scale=T)
D2 <- mahalanobis(pca$x, c(0, 0), cov=cov(pca$x))
plot(pca$x[,1], pca$x[,2], pch=20, cex=0.5, 
    col=c('red', 'black', 'blue')[ifelse(D2<1, 1, ifelse(D2<2, 2, 3))])
}


