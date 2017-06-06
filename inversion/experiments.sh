maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA modmix mod_Gauss 32 perturbed.model.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON modmix mod_Gauss 32 perturbed.model.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 BiSON BiSON CygAwball mod_Gauss 32 perturbed.CygA.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygAwball mod_Gauss 32 perturbed.CygA.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygA CygA CygA mod_Gauss 32 perturbed.CygA.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygB mod_Gauss 32 perturbed.CygB.names

maybe_sub.sh -p 9 Rscript OLA_optimize.R 0 CygB CygB CygBwball mod_Gauss 32 perturbed.CygB.names

