# Asteroseismology
Earl Bellinger's Ph.D. repository on forward and inverse problems in asteroseismology. 

---

## Workflow 

#### Forward 

`forward/`
  1. `python3 sobol_dispatcher.py` -- generate simulations with inputs varied in a quasi-random fashion 
    * `Rscript summarize.R` -- summarize each evolutionary track into a matrix table (this is called by `sobol_dispatcher`)
  2. `Rscript collate.R` -- collect nearly-evenly-spaced points from each summarized simulation into one big data file `simulations.dat`; this facilitates the inverse problem 
  3. `Rscript mesh.R` -- create meshes and scatterplot visualizations (e.g. H-R and C-D diagrams) on the collated tracks 
  4. `Rscript rankcorr.R` -- perform rank correlation tests and principal components analysis on the collated tracks 

#### Inverse

`inverse/`
  1. `Rscript tagesstern.R` -- degrade BiSON solar frequencies to the level of what is observable from the 16 Cyg stars for the sake of fair evaluation & comparison 
  2. `Rscript perturb.R` -- make Monte-Carlo perturbations of solar, Tagesstern, 16 Cyg, kages, and hares data to account for uncertainties in observed data 
  3. `python3 learn.py` -- learn what relates observable data to model properties from `../forward/simulations.dat` and predict the properties of the stars in `perturb/`
  4. `python3 constrain.py` -- constrain stellar systems like 16 Cyg by assuming that they have the same age and initial composition 

---

## Utilities 

#### Scripts

`scripts/`
- `maybe_sub.sh` -- shell script for submitting jobs to the condor queuing system 
- `dispatch.sh` -- shell script for generating a MESA evolutionary track 
- `mesa_template/` -- directory containing default instructions for a MESA evolutionary track (copied automatically by `dispatch`)
- `fgong2freqs.sh` -- shell script for redistributing a MESA model mesh and calculating adiabatic pulsation frequencies via ADIPLS 
- `seismology.R` -- R script for making seismological calculations from a frequencies data file 
- `utils.R` -- R utility script for plotting, constants, etc 
- `sobol_lib.py` -- python library for generating Sobol (quasi-random) numbers 

#### Miscellaneous

`misc/`
- `python3 plot_sph_harm.py` -- make spherical harmonics plots to visualize the pulsation frequencies of solar-like oscillators 
- `python3 plot_grids.py` -- make plots of linear, random, and quasi-random (Sobol) grids to justify the use of the latter 
- `python3 plot_classification.py` -- make plots of linear and non-linear (in this case, XOR) classification problems to illustrate the limitations and usefulness of basic and advanced ML routines 
