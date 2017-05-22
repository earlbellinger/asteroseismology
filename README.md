# Asteroseismology
[Earl Bellinger](http://earlbellinger.com)'s Ph.D. repository on forward and inverse problems in asteroseismology. 

If any of these programs are useful to you, please consider citing one or more of the following:

  * Bellinger, E. P., Angelou, G. C., Hekker, S., Basu, S., Ball, W., Guggenberger, E. (2016). [Fundamental Parameters of Main-Sequence Stars in an Instant with Machine Learning](http://adsabs.harvard.edu/abs/2016ApJ...830...31B). *The Astrophysical Journal*, 830 (1), 20.

  * Angelou, G. C., Bellinger, E. P., Hekker, S., Basu, S. (2017). [On the Statistical Properties of the Lower Main Sequence](http://adsabs.harvard.edu/abs/2017ApJ...839..116A). *The Astrophysical Journal*, 839 (2), 116.

  * Bellinger, E. P., Angelou, G., Hekker, S., Basu, S., Ball, W., Guggenberger, E. (2017). [Fundamental Parameters in an Instant with Machine Learning: Application to Kepler LEGACY Targets](https://arxiv.org/abs/1705.06759). Seismology of the Sun and Distant Stars, submitted.

See [Releases](https://github.com/earlbellinger/asteroseismology/releases) for the versions of this repository corresponding to those papers. 

---

## Workflow 

#### Forward 

`forward/`
  1. `mesa_template/` -- directory containing default instructions for a MESA evolutionary track (copied automatically by `dispatch`)
  
  
  2. `python3 sobol_dispatcher.py` -- generate initial conditions varied in a quasi-random fashion, calls the following files:
    * `./dispatch.sh` -- shell script for generating a parameterized MESA evolutionary track (e.g. -M 1 for a 1 solar mass model) 
    * `Rscript discontinuity.R` -- detect discontinuities in the simulated evolution 
    * `Rscript summarize.R` -- summarize an evolutionary track into a matrix 
    
    
  3. `Rscript collate.R` -- collect nearly-evenly-spaced points from each summarized simulation into one big data file `simulations.dat`; this facilitates the inverse problem
  
  
  4. `./amp.sh` -- Run `sobol_dispatcher.py` with settings that facilitate comparison with the asteroseismic modeling portal (AMP)



#### Analyze Models 

`analyze/`

  5. `Rscript diffusion.R` -- plots the initial and final simulation metallicities as a function of mass and diffusion
  
  
  6. `Rscript inputs.R` -- creates a diagram showing the initial conditions of the grid based off of `../forward/initial_conditions.dat` which is generated by `sobol_dispatcher.py`


#### Random Forest Regression 

`regression/`

  7. `Rscript tagesstern.R` -- degrade BiSON solar frequencies to the level of what is observable from the 16 Cyg stars for the sake of fair evaluation & comparison 

  8. `Rscript hare_compile.R` -- turn the Hare & Hound data into a format I can parse 
  
  9. `Rscript perturb.R` -- make Monte-Carlo perturbations of solar, Tagesstern, 16 Cyg, kages, and hares data to account for uncertainties in observed data 
  
  10. `python3 subsetter.py` -- determine the number of evolutionary tracks, models per evolutionary track, and trees in the forest that are needed for `learn.py` to work well 
    * `Rscript forest_evaluate.R` -- visualize the output of `subsetter.py`
  
  11. `python3 learn.py` -- learn what relates observable data to model properties from `../forward/simulations.dat` and predict the properties of the stars in `perturb/`
    * `Rscript importances.R` -- plots the feature importances of the random forests obtained in `learn.py`
    * `Rscript cyg.R` -- plots the predicted quantities of 16 Cyg from `learn.py` against literature values 
    * `Rscript us-vs-them.R`-- plots the predicted quantities of the KAGES stars and the Hare-and-Hound exercise against the literature values; also creates the diffusion plot for the KAGES stars 
  
  12. `Rscript legacy.R` -- plots the cumulative distribution functions for estimate uncertainties for the LEGACY targets

#### Asteroseismic Inversions 

`inversion/`

  * Coming soon...

---

## Utilities 

#### Scripts

`scripts/`
- `maybe_sub.sh` -- shell script for submitting jobs to the condor queuing system 
- `fgong2freqs.sh` -- shell script for redistributing a MESA model mesh and calculating adiabatic pulsation frequencies via ADIPLS 
- `seismology.R` -- R script for making seismological calculations from a frequencies data file 
- `utils.R` -- R utility script for plotting, constants, etc 
- `sobol_lib.py` -- python library for generating Sobol (quasi-random) numbers 
- `kerexact.sh` -- generates kernel functions from stellar models 

#### Miscellaneous

`misc/`
- `jcd-kasc/`
  * `Rscript CD_diagram.R` -- plot an asteroseismic H-R diagram from a grid of MESA/GYRE models and overplot LEGACY data points on it 
- `python3 plot_sph_harm.py` -- make spherical harmonics plots to visualize the pulsation frequencies of solar-like oscillators
- `python3 plot_grids.py` -- make plots of linear, random, and quasi-random (Sobol) grids to justify the use of the latter 
- `python3 plot_classification.py` -- make plots of linear and non-linear (in this case, XOR) classification problems to illustrate the limitations and usefulness of basic and advanced ML routines 
- `Rscript 16CygB.R` -- make an annotated power spectrum of 16 Cyg B
- `matlab animate_sph_harm.m` -- create animations of spherical harmonics 
- `Rscript interp_vs_reg.R` -- plot the difference between linear interpolation and regression
- `Rscript plot_nearly-even-spacing.R` -- show the result of the linear transport problem on finding nearly-evenly spaced points


