library(lpSolve)

Hc <- read.table('M=1.000000_Y=0.280000_Z=0.020508_alpha=2.250000.dat', 
                 header=1)$Hc
N <- 50
x <- sample(Hc, N)
n <- 10
ideal <- seq(max(x), min(x), length=n)

png('line.png', height=175, width=780, family="Helvetica")
par(xaxs='i',yaxs='i',mar=c(5,1,1,1))
plot(NA, xlim=c(0, 0.7), ylim=c(0, 1), axes=F, ann=F)
axis(1)
title(xlab=expression("Core Hydrogen Abundance"~H[c]))
points(x, rep(0.35, length(x)), pch=4, cex=1, xpd=NA, col='darkblue')
points(ideal, rep(0.1, n), pch=20, cex=1, xpd=NA)
legend("top", legend=c("Data", "Ideal Spacing"), pch=c(4, 20), 
       col=c('darkblue', 'black', 'darkred'), horiz=1)
dev.off()

cost.mat  <- outer(ideal, x, function(x, y) abs(x-y))
row.signs <- rep("==", n)
row.rhs   <- rep(1, n)
col.signs <- rep("<=", N)
col.rhs   <- rep(1, N)
sol <- lp.transport(cost.mat, "min", row.signs, row.rhs,
                    col.signs, col.rhs)$solution
y <- x[apply(sol, 1, which.max)]

png('line2.png', height=175, width=780, family="Helvetica")
par(xaxs='i',yaxs='i',mar=c(5,1,1,1))
plot(NA, xlim=c(0, 0.7), ylim=c(0, 1), axes=F, ann=F)
axis(1)
title(xlab=expression("Core Hydrogen Abundance"~H[c]))
points(x, rep(0.35, length(x)), pch=4, cex=1, col='darkblue', xpd=NA)
points(y, rep(0.225, length(y)), pch=3, cex=1, col='darkred', xpd=NA)
points(ideal, rep(0.1, 10), pch=20, cex=1, xpd=NA)
legend("top", legend=c("Data", "Ideal Spacing", "Nearly-Equidistant Data"), 
       pch=c(4, 20, 3), col=c('darkblue', 'black', 'darkred'), horiz=1)
dev.off()
