#### Helio- and astero-seismic inversions
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

set.seed(0)

### LIBRARIES 
source('../scripts/utils.R') 
source('models.R')
source('frequencies.R')
source('kernels.R')
source('OLA_invert.R')
source('OLA_plots.R')
#num_procs <- max(1, as.integer(Sys.getenv()[['OMP_NUM_THREADS']]))
#parallelStartMulticore(num_procs)

### CONSTANTS 
k.pair  = u_Y
rs      = seq(0.05, 0.3, 0.05) 
sampler = c(T)
models  = get_model_list() 
kern.interp.xs = seq(0, 1, 0.001)
targ.kern.type <- 'mod_Gauss'


targ.mode <- '-p_CygAwball-m_CygA-e_CygA-r_diffusion-mod_Gauss'

#freqs <- get_freqs(target.name='CygAwball', mode.set='CygA', 
#    error.set='CygA', perturb=F) 

#ref.mod <- get_model(freqs=NULL, model.name='diffusion', 
#    target.name='CygAwball', k.pair=k.pair, square.Ks=F)


load(file.path('save', paste0('ref.mod',         targ.mode)))
load(file.path('save', paste0('inv.lists',       targ.mode)))
load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
load(file.path('save', paste0('cross.lists',     targ.mode)))
load(file.path('save', paste0('MRs',             targ.mode)))
load(file.path('save', paste0('inv.params',      targ.mode)))

#model.names <- get(model.list.name) 
#model.list <- parallelMap(function(model.name) 
#        get_model(freqs=freqs, model.name=model.name, 
#            target.name=target.name, k.pair=k.pair, square.Ks=F), 
#    model.name=model.names)
#names(model.list) <- model.names

#make_plots(plot_ref_mods, paste0('ref_mods', targ.mode),
#    model.list=model.list, inversion=inversion, 
#    xlim=c(0, 0.35), ylim=c(0.7, 1.7)) 

inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)

make_plots_inversion_all(ref.mod, inversion, kern.interp.xs=kern.interp.xs,
    k.str=targ.mode, cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype",
    make_xlabs=c(F,F,T), 
    make_ylabs=c(T,T,T), 
    inversion_ylim=c(-0.03, 0.15),
    #inversion_ylim=c(-0.2, 0.1),
    #col.pal="#F46D43", sampler=c(F,F,T,F,F,F),
    #inversion_ylim=c(0, 0.15),
    cross_kern_ylim=c(-0.5, 0.3))







targ.mode <- '-p_CygBwball-m_CygB-e_CygB-r_diffusion-mod_Gauss'

#freqs <- get_freqs(target.name='CygAwball', mode.set='CygA', 
#    error.set='CygA', perturb=F) 

#ref.mod <- get_model(freqs=NULL, model.name='diffusion', 
#    target.name='CygAwball', k.pair=k.pair, square.Ks=F)


load(file.path('save', paste0('ref.mod',         targ.mode)))
load(file.path('save', paste0('inv.lists',       targ.mode)))
load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
load(file.path('save', paste0('cross.lists',     targ.mode)))
load(file.path('save', paste0('MRs',             targ.mode)))
load(file.path('save', paste0('inv.params',      targ.mode)))

#model.names <- get(model.list.name) 
#model.list <- parallelMap(function(model.name) 
#        get_model(freqs=freqs, model.name=model.name, 
#            target.name=target.name, k.pair=k.pair, square.Ks=F), 
#    model.name=model.names)
#names(model.list) <- model.names

#make_plots(plot_ref_mods, paste0('ref_mods', targ.mode),
#    model.list=model.list, inversion=inversion, 
#    xlim=c(0, 0.35), ylim=c(0.7, 1.7)) 

inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)

make_plots_inversion_all(ref.mod, inversion, kern.interp.xs=kern.interp.xs,
    k.str=targ.mode, cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype",
    make_xlabs=c(F,F,T), 
    make_ylabs=c(F,F,F), 
    inversion_ylim=c(-0.03, 0.15),
    #inversion_ylim=c(-0.2, 0.1),
    #col.pal="#F46D43", sampler=c(F,F,T,F,F,F),
    #inversion_ylim=c(0, 0.15),
    cross_kern_ylim=c(-0.4, 0.3))

#ref.mod <- get_model(freqs=NULL, model.name='CygAyoung', 
#    target.name='CygA', k.pair=k.pair, square.Ks=F)











target.name <- 'CygA'
mode.set <- 'CygA'
error.set <- 'CygA'
ref.mod.name <- 'CygAwball'
model.list.name <- 'perturbed.CygA.names'
targ.mode <- paste0(
    '-p_', target.name, 
    '-m_', mode.set, 
    '-e_', error.set, 
    '-r_', ref.mod.name, 
    paste0("-", targ.kern.type))
model.names <- get(model.list.name) 
model.list <- parallelMap(function(model.name) 
        get_model(freqs=NULL, model.name=model.name, 
            target.name=target.name, k.pair=k.pair, square.Ks=F), 
    model.name=model.names)
names(model.list) <- model.names

#targ.mode <- '-p_CygA-m_CygA-e_CygA-r_CygAwball-mod_Gauss'
load(file.path('save', paste0('ref.mod',         targ.mode)))
load(file.path('save', paste0('inv.lists',       targ.mode)))
load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
load(file.path('save', paste0('cross.lists',     targ.mode)))
load(file.path('save', paste0('MRs',             targ.mode)))
load(file.path('save', paste0('inv.params',      targ.mode)))
inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)
inversion$k.pair$f1.units <- bquote(cm^2~s^-2)
print_latex_table(inversion)
make_plots_inversion_all(ref.mod, inversion, kern.interp.xs=kern.interp.xs,
    k.str=targ.mode, cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype",
    make_xlabs=c(T,T,T), 
    make_ylabs=c(T,T,T), 
    caption="16 Cyg A",
    #inversion_ylim=c(-0.02, 0.2),
    inversion_ylim=c(-0.2, 0.1),
    #col.pal="#F46D43", sampler=c(F,F,T,F,F,F),
    #inversion_ylim=c(0, 0.15),
    cross_kern_ylim=c(-0.8, 0.3))
make_plots(plot_ref_mods, paste0('ref_mods', targ.mode),
    ref.mod=model.list[['CygAhighRlowM']],
    model.list=model.list, inversion=inversion, caption="16 Cyg A", 
    xlim=c(0, 0.35), ylim=c(0.7, 1.7), use.cairo=T, font="Palatino Linotype") 


make_plots(plot_ref_rel_diffs, paste0('ref_mods_rel', targ.mode),
    model.list=list(#CyA=model.list[['CygADiffusion']],
                    CygAhighRlowM=model.list[['CygAhighRlowM']], 
                    CygAlowRhighM=model.list[['CygAlowRhighM']]), 
        xlim=c(0, 0.45), ylim=c(-0.2, 0.1), col.pal=c(1, red), 
        #c(1, blue, orange), 
        rs=rs, inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
        cross.lists=cross.lists, inv.params=inv.params, 
        kern.interp.xs=kern.interp.xs, 
        caption="16 Cyg A", wide=F, tall=F, make_png=F, slides=F, 
        #xlim=c(0, 0.35), ylim=c(0.7, 1.7), 
        use.cairo=T, font="Palatino Linotype") 

CygA.mod <- model.list[['CygADiffusion']]
CygAyoung.mod <- get_model(freqs=NULL, model.name='CygAyoung', 
     target.name='CygA', k.pair=k.pair, square.Ks=F)
CygAyounger.mod <- get_model(freqs=NULL, model.name='CygAyounger', 
     target.name='CygA', k.pair=k.pair, square.Ks=F)
make_plots(plot_ref_rel_diffs, paste0('ref_mods_rel-young', targ.mode),
    model.list=list(CygA=CygA.mod,
                    CygAyoung=CygAyoung.mod, 
                    CygAyounger=CygAyounger.mod), 
        xlim=c(0, 0.45), ylim=c(-0.2, 0.1), col.pal=c(orange, '#551A8B', blue), 
        rs=rs, inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
        cross.lists=cross.lists, inv.params=inv.params, 
        kern.interp.xs=kern.interp.xs, 
        caption="16 Cyg A", wide=F, tall=F, make_png=F, slides=F, 
        #xlim=c(0, 0.35), ylim=c(0.7, 1.7), 
        use.cairo=T, font="Palatino Linotype") 



ref.mod <- get_model(freqs=NULL, model.name='CygAyoung', 
     target.name='CygA', k.pair=k.pair, square.Ks=F)
inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)
make_plots(plot_inversion, filename="CygAyoung", 
    model=ref.mod, inversion=inversion, kern.interp.xs=kern.interp.xs,
    k.str='CygAyoung', cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype",
    caption=expression(tau == 6~"Gyr"),
    ylim=c(-0.2, 0.1), xlim=c(0, 0.35))
ref.mod <- get_model(freqs=NULL, model.name='CygAyounger', 
     target.name='CygA', k.pair=k.pair, square.Ks=F)
inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)
make_plots(plot_inversion, filename="CygAyounger", 
    model=ref.mod, inversion=inversion, kern.interp.xs=kern.interp.xs,
    k.str='CygAyounger', cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype", make_ylab=F, 
    caption=expression(tau == 5~"Gyr"),
    ylim=c(-0.2, 0.1), xlim=c(0, 0.35))





target.name <- 'CygB'
mode.set <- 'CygB'
error.set <- 'CygB'
ref.mod.name <- 'CygBwball'
model.list.name <- 'perturbed.CygB.names'
targ.mode <- paste0(
    '-p_', target.name, 
    '-m_', mode.set, 
    '-e_', error.set, 
    '-r_', ref.mod.name, 
    paste0("-", targ.kern.type))
model.names <- get(model.list.name) 
model.list <- parallelMap(function(model.name) 
        get_model(freqs=NULL, model.name=model.name, 
            target.name=target.name, k.pair=k.pair, square.Ks=F), 
    model.name=model.names)
names(model.list) <- model.names

load(file.path('save', paste0('ref.mod',         targ.mode)))
load(file.path('save', paste0('inv.lists',       targ.mode)))
load(file.path('save', paste0('avg.kerns.lists', targ.mode)))
load(file.path('save', paste0('cross.lists',     targ.mode)))
load(file.path('save', paste0('MRs',             targ.mode)))
load(file.path('save', paste0('inv.params',      targ.mode)))
inversion <- lists_to_inversion(model=ref.mod, rs=rs, 
    inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
    cross.lists=cross.lists, inv.params=inv.params,
    kern.interp.xs=kern.interp.xs)
inversion$k.pair$f1.units <- bquote(cm^2~s^-2)
print_latex_table(inversion)
make_plots_inversion_all(ref.mod, inversion, kern.interp.xs=kern.interp.xs,
    k.str=targ.mode, cross.inset="bottomright", use.cairo=T, 
    font="Palatino Linotype",
    make_xlabs=c(T,T,T), 
    make_ylabs=c(F,T,T), 
    caption="16 Cyg B",
    #inversion_ylim=c(-0.02, 0.2),
    inversion_ylim=c(-0.2, 0.1),
    #col.pal="#F46D43", sampler=c(F,F,T,F,F,F),
    #inversion_ylim=c(0, 0.15),
    cross_kern_ylim=c(-0.8, 0.3))
make_plots(plot_ref_mods, paste0('ref_mods', targ.mode),
    model.list=model.list, inversion=inversion, 
    make_ylab=F, caption="16 Cyg B", 
    xlim=c(0, 0.35), ylim=c(0.7, 1.7), use.cairo=T, font="Palatino Linotype") 

make_plots(plot_ref_rel_diffs, paste0('ref_mods_rel', targ.mode),
    model.list=list(#CygB=model.list[['CygBDiffusion']],
                    CygBhighRlowM=model.list[['CygBhighRlowM']],
                    CygBlowRhighM=model.list[['CygBlowRhighM']]), 
        xlim=c(0, 0.45), ylim=c(-0.2, 0.1), col.pal=c(1, red),
        #c(1, blue, orange), 
        rs=rs, inv.lists=inv.lists, avg.kerns.lists=avg.kerns.lists, 
        cross.lists=cross.lists, inv.params=inv.params, 
        kern.interp.xs=kern.interp.xs, make_ylab=F, 
        caption="16 Cyg B", wide=F, tall=F, make_png=F, slides=F, 
        #xlim=c(0, 0.35), ylim=c(0.7, 1.7), 
        use.cairo=T, font="Palatino Linotype") 


