#### Helio- and astero-seismic inversions with the Backus-Gilbert method 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 
#### Department of Astronomy, Yale University 

### LIBRARIES 
library(Bolstad) # sintegral
library(pracma) # cumtrapzlibrary(parallel)
library(parallelMap)
library(magicaxis)
library(RColorBrewer)

get_ln <- function(mode) {
    data.frame(l=strsplit(strsplit(mode, '_')[[1]][1], '\\.')[[1]][2],
               n=strsplit(strsplit(mode, '_')[[1]][2], '\\.')[[1]][2])
}

make_symmetric <- function(mat) {
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    mat
}

get_square_Ks <- function(modes, K, x0=NULL) {
    Ks <- do.call(plyr:::rbind.fill, Map(function(mode_i) {
        data.frame(do.call(cbind, Map(function(mode_j) {
            if (mode_i > mode_j) return(0)
            K_ij.integrand <- K[[modes[mode_i]]] * K[[modes[mode_j]]]
            MOLA <- if (!is.null(x0)) (K$x - x0)**2 else 1
            sintegral(K$x, K_ij.integrand * MOLA)$value
        }, mode_j=1:length(modes))))
    }, mode_i=1:length(modes)))
    names(Ks) <- modes
    rownames(Ks) <- modes
    make_symmetric(Ks)
}

get_K.ints <- function(modes, K) {
    do.call(c, Map(function(mode) {
        sintegral(K$x, K[[mode]])$value
    }, mode=modes))
}

get_knots <- function(x, num.knots, degree=4) {
    x <- x[x>0]
    n <- num.knots-degree+1#(2*degree)
    c(rep(min(x), degree), 
      seq(min(x), max(x), length=n),
      rep(max(x), degree))
}

get_F_surf <- function(nus, use.BG=T, num.knots=0, nu_ac=5000, nu.deg=4) {
    if (num.knots > 0) nu.knots <- get_knots(nus$nu.x, num.knots, degree=nu.deg)
    do.call(plyr:::rbind.fill, with(nus, Map(function(ell, nn, nu, E, Q_norm) {
        mode <- paste0('l.', ell, '_', 'n.', nn)
        
        if (use.BG) {
            BG_surfs <- data.frame(do.call(cbind, Map(function(pow) 
                ifelse(nu != 0, (nu/nu_ac)**pow / E, 0), # Q_norm, 0), #
                pow=c(-2, 2))))
            names(BG_surfs) <- paste0('a_', 1:length(BG_surfs)) 
        } else BG_surfs <- NULL 
        
        if (num.knots > 0) {
            #deg <- (length(nu.knots)-length(unique(nu.knots)))/2 
            F_surfs <- data.frame(do.call(cbind, Map(function(knot) 
                    if (nu==0) 0 else B(nu, knot, nu.deg, nu.knots) / Q_norm, 
                knot=1:(length(nu.knots)-(nu.deg+1)))))
            names(F_surfs) <- paste0('F_', 1:length(F_surfs))
        } else F_surfs <- NULL
        
        if (!is.null(BG_surfs) && !is.null(F_surfs)) {
            cbind(BG_surfs, F_surfs)
        } else if (!is.null(BG_surfs) && is.null(F_surfs)) {
            BG_surfs
        } else if (is.null(BG_surfs) && !is.null(F_surfs)) {
            F_surfs
        } else {
            data.frame()
        }
    
    }, ell=l, nn=n, nu=nu.x, E=E, Q_norm=Q_norm)))
}

get_A.mat <- function(cross.term, error.sup, modes, nus, F_surf, 
        K=NULL, C=NULL, K_ijs=NULL, C_ijs=NULL, K.ints=NULL, x0=NULL) {
    if (length(K.ints) <= 1) {
        cat("Calculating integrals\n")
        K.ints <- get_K.ints(modes, K)
    }
    if (length(K_ijs) <= 1 || length(C_ijs) <= 1) {
        cat("Calculating square matrices\n")
        K_ijs <- get_square_Ks(modes, K, x0)
        C_ijs <- get_square_Ks(modes, C, x0)
    }
    
    first.ij <- K_ijs + cross.term * C_ijs + error.sup * diag(nus$d.r.diff**2)
    mat <- cbind(rbind(cbind(first.ij, K.ints), c(K.ints, 0)), rbind(F_surf, 0))
    surf.zeros <- matrix(0L, nrow=ncol(F_surf), ncol=ncol(F_surf)+1)
    tsurf <- cbind(t(F_surf), surf.zeros)
    colnames(tsurf) <- names(mat)
    rbind(mat, tsurf)
}

mod_Gauss <- function(r, r_0, width, normalization_factor=1) {
    normalization_factor * r * exp( -( (r-r_0)/width + width/(2*r_0) )**2 )
} 

mod_sinc <- function(r, r_0, width, normalization_factor=1) {
    normalization_factor * ifelse(r==r_0, 
        (1-r_0)**4 * r_0 * width * 1000, 
        (1-r  )**4 * r   * sin(width*1000*(r-r_0)) / (r-r_0))
}

get_target_kernel <- function(r_0, r_f, width.r_f, K.x, f.spl, targ.kern.type) {
    targ.fun <- if (targ.kern.type == 'mod_Gauss') {
        mod_Gauss
    } else if (targ.kern.type == 'mod_sinc') {
        mod_sinc        
    } else stop(paste('Invalid target kernel type:', targ.kern.type))
    width.k <- width.r_f * f.spl(r_0) / f.spl(r_f)
    norm_fac <- 1/sintegral(K.x, targ.fun(r=K.x, r_0=r_0, width=width.k, 
        normalization_factor=1))$value
    if (is.finite(norm_fac)) {
        sapply(K.x, function(r) targ.fun(r=r, r_0=r_0, width=width.k, 
            normalization_factor=norm_fac))
    } else 0*K.x
} 

get_SOLA_v.vec <- function(r_0, r_f, width.r_f, K, modes, f.spl, 
        targ.kern.type) {
    target_kernel <- get_target_kernel(r_0=r_0, r_f=r_f, 
                                       width.r_f=width.r_f, 
                                       K.x=K$x, f.spl=f.spl,
                                       targ.kern.type=targ.kern.type)
    get_SOLA_v.vec.(target_kernel=target_kernel, K=K, modes=modes)
}

get_SOLA_v.vec. <- function(target_kernel, modes, K) {
    sapply(modes, function(mode) {
        sintegral(K$x, K[[mode]] * target_kernel)$value
    })
}

get_averaging_kernel <- function(coefs, modes, K) { 
    rowSums(do.call(cbind, Map(function(mode.idx) { 
        coefs[mode.idx] * K[modes[mode.idx]]
    }, mode.idx=1:length(modes))))
}

get_target_kernels <- function(rs, width, f.spl, K.x, r_f=0) {
    do.call(cbind, parallelMap(function(target_radius) {
        get_target_kernel(r_0=target_radius, r_f=r_f, 
            width.r_f=width, f.spl=f.spl, K.x=K.x,
            targ.kern.type=targ.kern.type)
    }, target_radius=rs))
}

invert.OLA <- function(model, rs, cross.term, error.sup, width=NULL,
        use.BG=T, num.knots=0, F_surf=NULL,
        MOLA.K_ijs=NULL, MOLA.C_ijs=NULL, 
        targ.kern.type='mod_Gauss', targ_kerns=NULL,
        get_cross_kerns=T) {
    K <- model$k1
    C <- model$k2
    nus <- model$nus
    modes <- model$modes
    nu_ac <- model$nu_ac
    K.ints <- model$K.ints
    K_ijs <- model$K_ijs
    C_ijs <- model$C_ijs
    f.spl <- model$cs.spl
    
    if (is.null(F_surf)) F_surf <- 
        get_F_surf(nus, num.knots=num.knots, use.BG=use.BG, nu_ac=nu_ac)
    
    
    if (is.null(width)) { # MOLA
        
        cat("MOLA inversion\n")
        parornot <- if (length(rs) == 1) Map else parallelMap
        c.vecs <- do.call(cbind, parornot(function(x0) {
            cat(paste("Inverting radius", x0, '\n'))
            MOLA.K_ij <- if (!is.null(MOLA.K_ijs) & 
                             paste0(x0) %in% names(MOLA.K_ijs)) {
                MOLA.K_ijs[[paste0(x0)]]
            } else MOLA.K_ijs
            MOLA.C_ij <- if (!is.null(MOLA.C_ijs) & 
                             paste0(x0) %in% names(MOLA.C_ijs)) {
                MOLA.C_ijs[[paste0(x0)]]
            } else MOLA.C_ijs
            A.mat <- get_A.mat(cross.term=cross.term, 
                               error.sup=error.sup, 
                               modes=modes, nus=nus, K=K, C=C, 
                               K_ijs=MOLA.K_ij, C_ijs=MOLA.C_ij, 
                               K.ints=K.ints, x0=x0, 
                               F_surf=F_surf)
            A.inv <- solve(A.mat, system="LDLt", tol=0)
            
            v.vec <- c(rep(0, length(modes)), 1, 
                       rep(0, nrow(A.mat)-length(modes)-1))
            A.inv %*% v.vec
        }, x0=rs))
        
    } else { # SOLA
        
        cat("SOLA inversion\n")
        A.mat <- get_A.mat(cross.term=cross.term, error.sup=error.sup, 
                           modes=modes, nus=nus, K=K, C=C, 
                           K_ijs=K_ijs, C_ijs=C_ijs, 
                           K.ints=K.ints, x0=NULL, 
                           F_surf=F_surf)
        A.inv <- solve(A.mat, system="LDLt", tol=0)
        
        targ_kerns <- if (is.null(targ_kerns)) {
            get_target_kernels(rs=rs, width=width, f.spl=f.spl, K.x=K$x, r_f=0)
        }
        
        c.vecs <- do.call(cbind, parallelMap(function(targ_kern_i) {
            targ_kern <- targ_kerns[,targ_kern_i]
            v.vec <- get_SOLA_v.vec.(target_kernel=targ_kern, modes=modes, K=K)
            v.vec <- c(v.vec, 1, rep(0, nrow(A.mat)-length(v.vec)-1))
            A.inv %*% v.vec
        }, targ_kern_i=1:ncol(targ_kerns)))
        
        #checksum <- sum(model$K.ints * c.vecs[1:length(modes)]) == 1
        
    }
    
    inv.coefs <- as.matrix(c.vecs[1:nrow(nus),]) 
    df_dr <- colSums(sweep(inv.coefs, MARGIN=1, nus$r.diff, `*`)) 
    err <- sqrt(colSums(sweep(inv.coefs**2, MARGIN=1, nus$d.r.diff**2, `*`))) 
    
    avg_kerns <- do.call(cbind, parallelMap(function(coefs_j) {
        get_averaging_kernel(coefs=inv.coefs[,coefs_j], modes=modes, K=K)
    }, coefs_j=1:ncol(inv.coefs)))
    cross_kerns <- if (get_cross_kerns)
        do.call(cbind, parallelMap(function(coefs_j) {
            get_averaging_kernel(coefs=inv.coefs[,coefs_j], modes=modes, K=C)
        }, coefs_j=1:ncol(inv.coefs))) else NULL
    
    #avg_kerns <- apply(inv.coefs, 2, function(coefs) {
    #    get_averaging_kernel(coefs=coefs, modes=modes, K=K)
    #})
    #cross_kerns <- apply(inv.coefs, 2, function(coefs) {
    #    get_averaging_kernel(coefs=coefs, modes=modes, K=C)
    #})
    
    errbars <- do.call(rbind, apply(avg_kerns, 2, function(avg_kern) {
        cumulative <- cumtrapz(K$x, avg_kern)
        hm <- max(avg_kern)/2 # half maximum
        fwhm.mid <- K$x[which.max(avg_kern)]
        
        to.left <- K$x < fwhm.mid
        negs.left <- which(avg_kern[to.left] < 0)
        if (length(negs.left)==0) negs.left <- 1
        left <- K$x[to.left][max(negs.left)]
        if (is.na(left)) left <- 0
        
        to.right <- K$x > fwhm.mid
        negs.right <- which(avg_kern[to.right] < 0)
        if (length(negs.right)==0) negs.right <- length(to.right)
        right <- K$x[to.right][min(negs.right)]
        if (is.na(right)) right <- 1
        
        inside <- K$x > left & K$x < right
        inside.left <- inside & to.left
        inside.right <- inside & to.right
        fwhm.left <- if (!any(inside.left)) 0 else
            K$x[inside.left][which.min((avg_kern[inside.left]-hm)**2)]
        fwhm.right <- if (!any(inside.right)) 1 else 
            K$x[inside.right][which.min((avg_kern[inside.right]-hm)**2)]
        data.frame(r.first_q=K$x[which.min(cumulative<=0.25)], 
                   r.median=K$x[which.min(cumulative<=0.5)], 
                   r.third_q=K$x[which.min(cumulative<=0.75)],
                   fwhm.left=fwhm.left, 
                   fwhm.mid=fwhm.mid, 
                   fwhm.right=fwhm.right)
    }))
    
    m.f1 <- model$f1.spl(errbars$fwhm.mid)
    f <- m.f1 - df_dr * m.f1
    f.err <- err * abs(m.f1)
    result <- cbind(rs, df_dr, err, errbars, f, f.err)
    
    if ('d.f1.spl' %in% names(model)) {
        result <- cbind(result, data.frame(true_df_dr=model$d.f1.spl(rs)))
    }
    
    params <- if (is.null(width)) {
        data.frame(cross.term=cross.term, 
           error.sup=error.sup,
           num.knots=num.knots)
    } else { 
        data.frame(cross.term=cross.term, 
           error.sup=error.sup, 
           width=width,
           num.knots=num.knots)
    }
    
    list(inv.coefs=inv.coefs,
         avg_kerns=avg_kerns,
         targ_kerns=if (!is.null(width)) targ_kerns else NULL, 
         cross_kerns=cross_kerns, 
         params=params,
         targ.kern.type=targ.kern.type, 
         c.vecs=c.vecs,
         result=result)
}

minimize_dist <- function(model, rs, use.BG=T, num.knots=0,
                          initial_params=c(100, 100, 0.01), max.width=0.2,
                          targ.kern.type='mod_Gauss', max.err=0.2) {
    if (length(initial_params) < 3) { # MOLA
        MOLA.K_ijs <- parallelMap(function(r) 
            get_square_Ks(model$modes, model$k1, r), r=rs)
        MOLA.C_ijs <- parallelMap(function(r) 
            get_square_Ks(model$modes, model$k2, r), r=rs)
        names(MOLA.K_ijs) <- rs
        names(MOLA.C_ijs) <- rs
    } else {
        MOLA.K_ijs <- NULL
        MOLA.C_ijs <- NULL
    }
    d.f1.spl <- splinefun(model$r, model$d.f1.true)
    optim_result <- optim(log10(initial_params), fn=function(params) {
        params <- 10**params
        cat(paste("Trying", 
            c("beta:", "mu:", "width:", "knots:")[1:length(params)], 
            params,'\n'))
        cross.term <- params[1]
        error.sup <- if (length(params) > 1) params[2] else 0
        width <- if (length(params) > 2) params[3] else NULL
        num.knots <- if (length(params) > 3) ceil(params[4]) else 0
        if (num.knots < 0) num.knots <- 0
        if (!is.null(width) && width > max.width) return(Inf) 
        inversion <- tryCatch({
            invert.OLA(model=model, rs=rs, 
                cross.term=cross.term, error.sup=error.sup, 
                width=width, use.BG=use.BG, num.knots=num.knots,
                targ.kern.type=targ.kern.type,
                MOLA.K_ijs=MOLA.K_ijs, MOLA.C_ijs=MOLA.C_ijs)
            }, error = function(e) 
                list(result=data.frame(fwhm.mid=0, df_dr=Inf)))
        
        #loc_avg <- apply(inversion$avg_kerns, 2, function(avg_kern) {
        #    sintegral(model$k1$x[-1], model$d.f1.true * avg_kern[-1])$value
        #})
        #chi2 <- sum(abs(d.f1.spl(rs) - loc_avg))
        
        #cross_terms <- apply(inversion$cross_kerns, 2, function(cross_kern) {
        #    sintegral(model$k1$x[-1], abs(cross_kern[-1]))$value
        #})
        
        result <- inversion$result 
        
        #chi2 <- with(result, sum( ((d.f1.spl(fwhm.mid) - df_dr)/err)**2 ) )
        chi2 <- with(result, 
            sum( ( (d.f1.spl(fwhm.mid) - df_dr)**2 ) ) #+ #* err ) #+
            #
            #sum( abs( (d.f1.spl(fwhm.mid) - df_dr) ) ) + #* err ) #+
            #sum( abs( (d.f1.spl(rs) - df_dr) ) ) + #* err ) #+
            #sum( abs( (fwhm.mid - rs) ) ) #+
            #
            #sum( (d.f1.spl(fwhm.mid) - df_dr)**2 ) +
            #sum( (fwhm.mid - rs) ) #+
            #
            #sum( ((m1$m2.f1.spl(fwhm.mid) - f) / f.err)**2 ) + 
            #sum( ((d.f1.spl(fwhm.mid) - df_dr) / err)**2 ) + 
            #
            #sum( ((fwhm.mid - rs) / (fwhm.right - fwhm.left)/2)**2 ) #+
            #sum( f.err**2 )
        )
        #chi2 <- sum(abs(d.f1.spl(result$fwhm.mid) - result$df_dr)**2) #+
            #with(result, sum((fwhm.mid-rs)**2))
        #**2 / result$err)
        #chi2 <- sum(abs(d.f1.spl(result$fwhm.mid)-result$df_dr)**2 / 
        #    result$err + result$err)
        
        #chi2 <- chi2 / (result$f / result$f.err)
        #if (any(result$err > max.err)) chi2 <- 1000*chi2
        
        cat(paste("Result:", chi2, '\n'))
        chi2
    }, control=list(abstol=0.001, reltol=0.00001))
    cat("Optimization complete\n")
    #, control=list(trace=999, parscale=c(1000, 0.1, 100)))
    best_params <- 10**optim_result$par
    inversion <- invert.OLA(model=model, rs=rs, 
               cross.term=best_params[1], 
               error.sup=best_params[2], 
               width=if (length(best_params) >= 3) 
                   best_params[3] else NULL, 
               use.BG=use.BG, 
               num.knots=if (length(best_params)>=4) 
                   best_params[4] else num.knots, 
               targ.kern.type=targ.kern.type,
               MOLA.K_ijs=MOLA.K_ijs, MOLA.C_ijs=MOLA.C_ijs) 
    inversion$score <- optim_result$value
    inversion
}

minimize_dist_individual <- function(model, rs, use.BG=T, num.knots=0,
        initial_params=c(100, 100, 0.01), max.width=0.2,
        targ.kern.type='mod_Gauss', max.err=0.1) {
    
    inversion.results <- parallelMap(function(r)
        minimize_dist(model=model, rs=r, use.BG=use.BG, num.knots=num.knots,
            initial_params=initial_params, max.width=max.width, 
            targ.kern.type=targ.kern.type, max.err=max.err),
        r=rs)
    
    Map(function(name) {
        combine <- if (name == "result" || name == "params") rbind else cbind
        do.call(combine, Map(function(result) result[[name]],
            result=inversion.results))
    }, name=names(inversion.results[[1]]))
}

