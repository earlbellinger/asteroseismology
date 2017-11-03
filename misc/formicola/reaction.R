# 14N(p,g)15O reaction rate given by Formicola et al. 2004

z1 <- 7
z2 <- 8
m1 <- 14.0067
m2 <- 15.999

SE0 <- 1.7 # KeV 

rate <- 7.83 * 10**9 ((z1 * z2)/(A*T9))**(1/3) * SE0 * MeV_barn *
    exp(-4.2487 * ((z1**2*z2**2*A)/T9)**(1/3))
