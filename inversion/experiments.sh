maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON CygAwball mod_Gauss 128 perturbed.CygA.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygAwball mod_Gauss 128 perturbed.CygA.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygA mod_Gauss 128 perturbed.CygA.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON modmix mod_Gauss 128 perturbed.model.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA modmix mod_Gauss 128 perturbed.model.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA BiSON mod_Gauss 128 perturbed.model.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON CygB mod_Gauss 128 perturbed.CygB.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygBwball mod_Gauss 128 perturbed.CygB.names
maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygB mod_Gauss 128 perturbed.CygB.names



maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 BiSON BiSON modmix mod_Gauss 128 perturbed.model.names 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 BiSON BiSON CygAwball mod_Gauss 128 perturbed.CygA.names 1.08 1.22 

maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA CygA modmix mod_Gauss 128 perturbed.model.names 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygA CygA CygAwball mod_Gauss 128 perturbed.CygA.names 1.08 1.22 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R 0 CygB CygB CygBwball mod_Gauss 128 perturbed.CygB.names 1.03 1.12 0.015 0.02 



## TESTS ON MODELS 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA modmix 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygAwball diffusion perturbed.CygA.names mod_Gauss 1.08 1.22 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygBwball diffusion perturbed.CygB.names mod_Gauss 1.03 1.12 0.015 0.02 


## TESTS ON THE SUN 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R BiSON BiSON BiSON diffusion perturbed.model.names mod_Gauss 1 1 0.016 0.02 128 F 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygA diffusion perturbed.model.names mod_Gauss 1 1 0.016 0.02 128 F 


## REAL STARS
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygA CygA CygA CygAwball perturbed.CygA.names mod_Gauss 1.08 1.22 0.016 0.02 128 F 
maybe_sub.sh -m 6000000 -p 9 Rscript OLA_optimize_MR2.R CygB CygB CygB CygBwball perturbed.CygB.names mod_Gauss 1.03 1.12 0.015 0.02 128 F 


