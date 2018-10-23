options(scipen=5)
args <- commandArgs(TRUE)
filename <- if (length(args)>0) args[1] else file.path('LOGS', 'history.data')
if (dir.exists(filename)) {
    dir_files <- list.files(filename)
    DF <- data.frame()
    for (log_dir in file.path(filename, dir_files[grep('LOGS_', dir_files)])) {
        DF <- plyr:::rbind.fill(DF, 
            read.table(file.path(log_dir, 'history.data'), 
                header=TRUE, skip=5))
    }
    DF <- DF[order(DF$star_age),]
} else if (!file.exists(filename)) {
    exit(paste(filename), 'not found') 
} else {
    DF <- read.table(filename, header=1, skip=5)
}

# clip PMS
if (!('log_L' %in% names(DF))) {
    DF$log_L <- log10(DF$luminosity)
}
decreasing_L <- which(diff(DF$log_L) < 0 & DF$center_h1[-1] > 0.6)
if (any(decreasing_L)) {
    pms <- max(decreasing_L)
    print(paste("Clipping", pms, "points"))
    DF <- DF[-1:-pms,]
}


### Obtain observable properties from models 
summarize <- function(pro_file, freqs_file, ev.DF, dname) {
    print(pro_file)
    
    pro_body <- read.table(pro_file, header=TRUE, nrows=1, skip=5)
    hstry <- ev.DF[ev.DF$model_number==pro_header$model_number,]
    if (nrow(hstry) == 0) {#|| hstry$mass_conv_core > 0) 
        print(paste("Model", pro_file, "failed: no history information"))
        return(NULL)
    }
    
    obs.DF <- NULL
    
    obs.DF["age"] <- hstry$star_age / 10**9
    obs.DF["M_current"] <- hstry$star_mass
    
    obs.DF["mass_cc"] <- hstry$mass_conv_core/hstry$star_mass
    obs.DF["mass_he_core"] <- hstry$he_core_mass
    
    obs.DF["h_exh_core_mass"] <- hstry$h_exh_core_mass
    obs.DF["h_exh_core_radius"] <- hstry$h_exh_core_radius
    
    obs.DF["radius_cc"] <- if (obs.DF['mass_cc'] <= 0) 0 else {
        if (hstry$conv_mx1_bot_r <= 0.1) {
            hstry$conv_mx1_top_r / 10**(hstry$log_R)
        } else {
            hstry$conv_mx2_top_r / 10**(hstry$log_R)
        }
    }
    
    obs.DF["log_LH"] <- hstry$log_LH
    obs.DF["log_LHe"] <- hstry$log_LHe
    obs.DF["log_center_T"] <- hstry$log_center_T
    obs.DF["log_center_Rho"] <- hstry$log_center_Rho
    obs.DF["log_center_P"] <- hstry$log_center_P
    obs.DF["center_mu"] <- hstry$center_mu
    obs.DF["center_degeneracy"] <- hstry$center_degeneracy
    
    obs.DF["X_c"]  <- hstry$center_h1 + hstry$center_h2
    obs.DF["Y_c"]  <- hstry$center_he3 + hstry$center_he4 
    obs.DF["Li_c"] <- hstry$center_li7
    obs.DF["Be_c"] <- hstry$center_be7
    obs.DF["B_c"]  <- hstry$center_b8
    obs.DF["C_c"]  <- hstry$center_c12 + hstry$center_c13
    obs.DF["N_c"]  <- hstry$center_n13 + hstry$center_n14 + hstry$center_n15
    obs.DF["O_c"]  <- hstry$center_o14 + hstry$center_o15 + hstry$center_o16 +
        hstry$center_o17 + hstry$center_o18
    obs.DF["F_c"]  <- hstry$center_f17 + hstry$center_f18 + hstry$center_f19
    obs.DF["Ne_c"] <- hstry$center_ne18 + hstry$center_ne19 + 
        hstry$center_ne20 + hstry$center_ne22
    obs.DF["Mg_c"] <- hstry$center_mg22 + hstry$center_mg24
    
    obs.DF["H1_c"]   <- hstry$center_h1 
    obs.DF["H2_c"]   <- hstry$center_h2 
    obs.DF["He3_c"]  <- hstry$center_he3 
    obs.DF["He4_c"]  <- hstry$center_he4 
    obs.DF["Li7_c"]  <- hstry$center_li7 
    obs.DF["Be7_c"]  <- hstry$center_be7 
    obs.DF["B8_c"]   <- hstry$center_b8 
    obs.DF["C12_c"]  <- hstry$center_c12 
    obs.DF["C13_c"]  <- hstry$center_c13 
    obs.DF["N13_c"]  <- hstry$center_n13 
    obs.DF["N14_c"]  <- hstry$center_n14 
    obs.DF["N15_c"]  <- hstry$center_n15 
    obs.DF["O14_c"]  <- hstry$center_o14 
    obs.DF["O15_c"]  <- hstry$center_o15 
    obs.DF["O16_c"]  <- hstry$center_o16 
    obs.DF["O17_c"]  <- hstry$center_o17 
    obs.DF["O18_c"]  <- hstry$center_o18 
    obs.DF["F17_c"]  <- hstry$center_f17 
    obs.DF["F18_c"]  <- hstry$center_f18 
    obs.DF["F19_c"]  <- hstry$center_f19 
    obs.DF["Ne18_c"] <- hstry$center_ne18 
    obs.DF["Ne19_c"] <- hstry$center_ne19 
    obs.DF["Ne20_c"] <- hstry$center_ne20 
    obs.DF["Ne22_c"] <- hstry$center_ne22 
    obs.DF["Mg22_c"] <- hstry$center_mg22 
    obs.DF["Mg24_c"] <- hstry$center_mg24 
    
    obs.DF["X_surf"]  <- hstry$surface_h1 + hstry$surface_h2
    obs.DF["Y_surf"]  <- hstry$surface_he3 + hstry$surface_he4 
    obs.DF["Li_surf"] <- hstry$surface_li7
    obs.DF["Be_surf"] <- hstry$surface_be7
    obs.DF["B_surf"]  <- hstry$surface_b8
    obs.DF["C_surf"]  <- hstry$surface_c12 + hstry$surface_c13
    obs.DF["N_surf"]  <- hstry$surface_n13 + hstry$surface_n14 + 
        hstry$surface_n15
    obs.DF["O_surf"]  <- hstry$surface_o14 + hstry$surface_o15 + 
        hstry$surface_o16 + hstry$surface_o17 + hstry$surface_o18
    obs.DF["F_surf"]  <- hstry$surface_f17 + hstry$surface_f18 + 
        hstry$surface_f19
    obs.DF["Ne_surf"] <- hstry$surface_ne18 + hstry$surface_ne19 + 
        hstry$surface_ne20 + hstry$surface_ne22
    obs.DF["Mg_surf"] <- hstry$surface_mg22 + hstry$surface_mg24
    
    obs.DF["H1_surf"]   <- hstry$surface_h1 
    obs.DF["H2_surf"]   <- hstry$surface_h2 
    obs.DF["He3_surf"]  <- hstry$surface_he3 
    obs.DF["He4_surf"]  <- hstry$surface_he4 
    obs.DF["Li7_surf"]  <- hstry$surface_li7 
    obs.DF["Be7_surf"]  <- hstry$surface_be7 
    obs.DF["B8_surf"]   <- hstry$surface_b8 
    obs.DF["C12_surf"]  <- hstry$surface_c12 
    obs.DF["C13_surf"]  <- hstry$surface_c13 
    obs.DF["N13_surf"]  <- hstry$surface_n13 
    obs.DF["N14_surf"]  <- hstry$surface_n14 
    obs.DF["N15_surf"]  <- hstry$surface_n15 
    obs.DF["O14_surf"]  <- hstry$surface_o14 
    obs.DF["O15_surf"]  <- hstry$surface_o15 
    obs.DF["O16_surf"]  <- hstry$surface_o16 
    obs.DF["O17_surf"]  <- hstry$surface_o17 
    obs.DF["O18_surf"]  <- hstry$surface_o18 
    obs.DF["F17_surf"]  <- hstry$surface_f17 
    obs.DF["F18_surf"]  <- hstry$surface_f18 
    obs.DF["F19_surf"]  <- hstry$surface_f19 
    obs.DF["Ne18_surf"] <- hstry$surface_ne18 
    obs.DF["Ne19_surf"] <- hstry$surface_ne19 
    obs.DF["Ne20_surf"] <- hstry$surface_ne20 
    obs.DF["Ne22_surf"] <- hstry$surface_ne22 
    obs.DF["Mg22_surf"] <- hstry$surface_mg22 
    obs.DF["Mg24_surf"] <- hstry$surface_mg24 
    
    obs.DF["log_Teff"] <- hstry$log_Teff
    obs.DF["log_L"] <- hstry$log_L
    obs.DF["log_R"] <- hstry$log_R
    obs.DF["log_g"] <- hstry$log_g
    
    obs.DF["Fe_H"] <- log10(10**hstry$log_surf_cell_z / 
            hstry$surface_h1 / Z_div_X_solar)
    
    obs.DF["cz_base"] <- hstry$cz_bot_radius
    obs.DF["acoustic_cutoff"] <- hstry$acoustic_cutoff
    obs.DF["acoustic_radius"] <- hstry$acoustic_radius
    
    obs.DF["surface_mu"] <- pro_body$mu
    obs.DF["delta_Pg_asym"] <- hstry$delta_Pg
    obs.DF["nu_max_classic"] <- scaling_nu_max(R=obs.DF[["radius"]], 
        M=hstry[["star_mass"]], Teff=obs.DF[["Teff"]])
    obs.DF["nu_max"] <- scaling_nu_max_Viani(R=obs.DF[["radius"]], 
        M=hstry[["star_mass"]], Teff=obs.DF[["Teff"]], 
        mu=obs.DF[["surface_mu"]]) 
    obs.DF["delta_nu_asym"] <- hstry$delta_nu
    
    freqs <- parse_freqs(freqs_file, gyre=T)
    seis <- seismology(freqs, nu_max=obs.DF[["nu_max_classic"]],
        min_points=3, check_nu_max=T)
    if ("Dnu0" %in% names(seis)) obs.DF["Dnu0_classic"] <- seis[["Dnu0"]]
    seis.DF <- seismology(freqs, nu_max=obs.DF[["nu_max"]],
        min_points=3, check_nu_max=T)
    
    #as.data.frame(cbind(obs.DF, seis.DF))
    merge(rbind(obs.DF), rbind(seis.DF))
}

