library(randomForest)
library(magicaxis)

get_data <- function(perturb, cov) {
    print(c(perturb, cov))
    a <- read.table(perturb, header=1)[1:100,]
    b <- read.table(cov, header=1)[1:100,]
    print(c(nrow(a), nrow(b)))
    if (nrow(a) != nrow(b)) return(NULL)
    cbind(a, b)
}

perturbs <- c(file.path('perturb', 
    grep('.dat', list.files('perturb'), value=T)), 
    file.path('perturb', 'kages', 
                     list.files(file.path('perturb', 'kages'))))
perturbs <- perturbs[order(basename(perturbs))]
covs <- file.path('learn_covs-simulations', 'kages', 
            list.files(file.path('learn_covs-simulations', 'kages')))
covs <- covs[order(basename(covs))]

DF <- do.call(plyr:::rbind.fill, 
    Map(function(perturb, cov) get_data(perturb, cov), 
    perturb=perturbs, cov=covs))
#data <- data[,c(-26:-30)]
#data <- data[complete.cases(data),]



data <- read.table('predict.dat', header=1)

DF <- data[,!grepl("^d", names(data))]
DF.unc <- cbind(data.frame(Name=data$Name), data[,grepl("^d", names(data))])

new_vals <- NULL
new_uncs <- NULL

for (filename in grep('-obs.dat$', 
        c(file.path('data', list.files(file.path('data'))), 
          file.path('data', 'kages', list.files(file.path('data', 'kages')))), 
        value=1)) {
    Name <- sub('-.+', '', basename(filename))
    print(Name)
    if (!any(DF$Name == Name)) next
    obs <- read.table(filename, header=1)
    vals <- data.frame(rbind(obs$value))
    uncs <- data.frame(rbind(obs$uncertainty))
    names(vals) <- sub("Fe/H", "FeH", obs$name)
    names(uncs) <- sub("Fe/H", "dFeH", obs$name)
    names(uncs) <- sub("Teff", "dTeff", obs$name)
    names(uncs) <- sub("nu_max", "dTeff", obs$name)
    vals <- cbind(Name, vals)
    uncs <- cbind(Name, uncs)
    
    new_vals <- plyr:::rbind.fill(new_vals, vals)
    new_uncs <- plyr:::rbind.fill(new_uncs, uncs)
    #DF <- merge(DF, vals, all=TRUE)
    #DF.unc <- merge(DF.unc, uncs, all=TRUE)
    #print(DF)
    #pause()
}

DF <- merge(DF, new_vals)
DF.unc <- merge(DF.unc, new_uncs)

selection <- !apply(is.na(DF), 2, any)
DF <- DF[,selection]
DF.unc <- DF.unc[,selection]

#library(rgp)
#result1 <- symbolicRegression(D~., data=DF[,-1], 
#    populationSize=1000, individualSizeLimit=128,
#    breedingTries=1000,
#    extinctionPrevention=1,
#    penalizeGenotypeConstantIndividuals=1)

attach(DF)
attach(DF.unc)

plot(NA, xlab="Mass M", ylab="Diffusion D",
    xlim=c(min(M-dM), max(M+dM)), 
    ylim=c(min(D-dD), max(D+dD)))
arrows(M, D-dD, M, D+dD, length=0, angle=90, code=3, col="gray")
arrows(M-dM, D, M+dM, D, length=0, angle=90, code=3, col="gray")
points(M, D, pch=19)

plot(NA, axes=FALSE,
    xlab="Metallicity Z", ylab="Diffusion D",
    xlim=c(min(Z-dZ), max(Z+dZ)), 
    ylim=c(min(D-dD), max(D+dD)),
    xaxs='i', yaxs='i')
magaxis(side=1:4, family="Palatino", tcl=0.25, labels=c(1,1,0,0))
arrows(Z, D-dD, Z, D+dD, length=0, angle=90, code=3, col="gray")
arrows(Z-dZ, D, Z+dZ, D, length=0, angle=90, code=3, col="gray")
points(Z, D, pch=19)
relation <- lm(D~I(Z**3), data=DF)
new_Z <- seq(0, 0.04, 0.0001)
lines(new_Z, predict(relation, newdata=data.frame(Z=new_Z)))
legend("topright", bty='n',
   legend=as.expression(bquote(
       D == .(coef(relation)[2]) %*% Z + .(coef(relation)[1]) )))

plot(NA, axes=FALSE,
    xlab="Metallicity Y", ylab="Diffusion D",
    xlim=c(min(Y-dY), max(Y+dY)), 
    ylim=c(min(D-dD), max(D+dD)),
    xaxs='i', yaxs='i')
magaxis(side=1:4, family="Palatino", tcl=0.25, labels=c(1,1,0,0))
arrows(Y, D-dD, Y, D+dD, length=0, angle=90, code=3, col="gray")
arrows(Y-dY, D, Y+dY, D, length=0, angle=90, code=3, col="gray")
points(Y, D, pch=19)



