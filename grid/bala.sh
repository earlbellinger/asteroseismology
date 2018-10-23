maybe_sub.sh -p 2 ./dispatch.sh -d sun -n sunM2 -M 2 -D 0 -g 0
maybe_sub.sh -p 2 ./dispatch.sh -d sun -n sunM3 -M 3 -D 0 -g 0

maybe_sub.sh -p 1 ./dispatch.sh -d bala -n mp10 -M 1 -D 0 -g 0 -Y 0.258387 -Z 0.00842595
maybe_sub.sh -p 1 ./dispatch.sh -d bala -n mp15 -M 1.5 -D 0 -g 0 -Y 0.258387 -Z 0.00842595
maybe_sub.sh -p 1 ./dispatch.sh -d bala -n mp20 -M 2 -D 0 -g 0 -Y 0.258387 -Z 0.00842595

DF <- read.table('LOGS_BALA/history.data', header=1, skip=5)
DF2 <- with(DF, data.frame(star_age=star_age, star_mass=star_mass, cz_bot_radius=cz_bot_radius, log_L=log_L, log_Teff=log_Teff, log_R=log_R, center_h1=center_h1, center_h2=center_h2, center_he3=center_he3, center_he4=center_he4, center_li7=center_li7, center_be7=center_be7, center_b8=center_b8, center_c12=center_c12, center_c13=center_c13, center_n13=center_n13, center_n14=center_n14, center_n15=center_n15, center_o14=center_o14, center_o15=center_o15, center_o16=center_o16, center_o17=center_o17, center_o18=center_o18, center_f17=center_f17, center_f18=center_f18, center_f19=center_f19, center_ne18=center_ne18, center_ne19=center_ne19, center_ne20=center_ne20, center_ne22=center_ne22, center_mg22=center_mg22, center_mg24=center_mg24, surface_h1=surface_h1, surface_h2=surface_h2, surface_he3=surface_he3, surface_he4=surface_he4, surface_li7=surface_li7, surface_be7=surface_be7, surface_b8=surface_b8, surface_c12=surface_c12, surface_c13=surface_c13, surface_n13=surface_n13, surface_n14=surface_n14, surface_n15=surface_n15, surface_o14=surface_o14, surface_o15=surface_o15, surface_o16=surface_o16, surface_o17=surface_o17, surface_o18=surface_o18, surface_f17=surface_f17, surface_f18=surface_f18, surface_f19=surface_f19, surface_ne18=surface_ne18, surface_ne19=surface_ne19, surface_ne20=surface_ne20, surface_ne22=surface_ne22, surface_mg22=surface_mg22, surface_mg24=surface_mg24))
DF2$star_age <- with(DF2, (star_age-min(star_age))/10**9)
write.table(DF2, 'track-M2.dat', quote=F, row.names=F, sep='\t')
