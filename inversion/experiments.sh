maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON CygAwball mod_Gauss 128 perturbed.CygA.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygAwball mod_Gauss 128 perturbed.CygA.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygA mod_Gauss 128 perturbed.CygA.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON modmix mod_Gauss 128 perturbed.model.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA modmix mod_Gauss 128 perturbed.model.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA BiSON mod_Gauss 128 perturbed.model.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON CygB mod_Gauss 128 perturbed.CygB.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygBwball mod_Gauss 128 perturbed.CygB.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygB mod_Gauss 128 perturbed.CygB.names



maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA BiSON modmix mod_Gauss 128 perturbed.model.names 1 1 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA BiSON CygAwball mod_Gauss 128 perturbed.CygA.names 1.08 1.22 

maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA CygA modmix mod_Gauss 128 perturbed.model.names 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA CygA CygAwball mod_Gauss 128 perturbed.CygA.names 1.08 1.22 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygB CygB CygBwball mod_Gauss 128 perturbed.CygB.names 1.03 1.12 0.015 0.02 


mode.set        <- if (length(args)>0)   args[1] else 'CygA'
error.set       <- if (length(args)>1)   args[2] else 'CygA'
target.name     <- if (length(args)>2)   args[3] else 'CygAwball'
ref.mod.name    <- if (length(args)>3)   args[4] else 'diffusion'
star.name       <- if (length(args)>4)   args[5] else 'CygA'
targ.kern.type  <- if (length(args)>5)   args[6] else 'mod_Gauss'
initial.M       <- if (length(args)>6)   as.numeric(args[7])  else 1.08
initial.R       <- if (length(args)>7)   as.numeric(args[8])  else 1.22
sigma.M         <- if (length(args)>8)   as.numeric(args[9])  else 0.016
sigma.R         <- if (length(args)>9)   as.numeric(args[10]) else 0.02
n_trials        <- if (length(args)>10)  as.numeric(args[11]) else 128
perturb         <- if (length(args)>11)  as.logical(args[12]) else T
cov.MR          <- if (length(args)>12)  as.numeric(args[13]) else 0
cov.RM          <- if (length(args)>13)  as.numeric(args[14]) else 0
subtract.mean   <- if (length(args)>14)  as.logical(args[15]) else F

## TESTS ON MODELS 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB modmix diffusion Sun mod_Gauss 1 1 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygAwball diffusion CygA mod_Gauss 1.08 1.22 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygBwball diffusion CygB mod_Gauss 1.03 1.12 0.015 0.02 


## TESTS ON THE SUN 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB BiSON diffusion Sun mod_Gauss 1 1 0.016 0.02 128 F 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R BiSON BiSON BiSON diffusion Sun mod_Gauss 1 1 0.016 0.02 128 F 
#maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygB diffusion Sun mod_Gauss 1 1 0.016 0.02 128 F 


## REAL STARS
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygA CygAwball CygA mod_Gauss 1.08 1.22 0.016 0.02 128 F 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygB CygBwball CygB mod_Gauss 1.03 1.12 0.015 0.02 128 F 


#maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R KIC_8006161 KIC_8006161 KIC_8006161 8006161_meanRmeanM 8006161 mod_Gauss 1.004 0.952 0.026 0.021 128 F 
#maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R KIC_6106415 KIC_6106415 KIC_6106415 6106415_meanRmeanM 6106415 mod_Gauss 1.1 1.289 0.043 0.037 128 F 

# test mean subtraction

## TESTS ON MODELS 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB modmix diffusion Sun mod_Gauss 1 1 0.016 0.02 128 T 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygAwball diffusion CygA mod_Gauss 1.08 1.22 0.016 0.02 128 T 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygBwball diffusion CygB mod_Gauss 1.03 1.12 0.015 0.02 128 T 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB BiSON diffusion Sun mod_Gauss 1 1 0.016 0.02 128 F 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R BiSON BiSON BiSON diffusion Sun mod_Gauss 1 1 0.016 0.02 128 F 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygA CygAwball CygA mod_Gauss 1.08 1.22 0.016 0.02 128 F 0 0 T
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygB CygBwball CygB mod_Gauss 1.03 1.12 0.015 0.02 128 F 0 0 T 


