/* Inference in the sampling-covariance scale (no sqrt(N) rescaling).

   Focused inference uses a conservative outer envelope, NOT a finite grid
   masquerading as the sharp two-step extremum.  For candidate loading a,
   target t, nuisance mean mu, and Gaussian sampling error e, the error is

                  a*e + (a-t)*mu.

   The nonzero entries of a-t occur only in the precontrast block.  Over
   an ellipsoid centered at observed h, its support function is analytic.
   Taking extrema over ALL candidates also covers the selected candidate.
   For an additive output, summing component extrema is conservative even
   when the same nuisance ellipsoid applies to every component.  Each
   component nevertheless uses the SAME joint Gaussian draw, retaining
   all cross-cohort covariance.

   Monte Carlo error is handled with two DKW bounds and a union bound.
   Thus these are conservative Gaussian-limit marginal intervals.  They
   are not simultaneous intervals, and are not finite-sample coverage
   guarantees for an estimated covariance matrix. */

real rowvector jw_infer_normal(real scalar estimate, real scalar variance,
                               real scalar level)
{
    real scalar se, critical, pv
    se = sqrt(variance)
    critical = invnormal(1-(1-level/100)/2)
    if (se==0) pv = (estimate==0)
    else pv = min((1,2*normal(-abs(estimate/se))))
    return((estimate,se,estimate-critical*se,estimate+critical*se,pv))
}

/* Return S such that S*S' equals V, allowing genuine zero eigenvalues.
   Only negative eigenvalues at numerical roundoff scale are clipped. */
real matrix jw_infer_sqrt(real matrix V)
{
    real matrix E
    real rowvector lambda
    real scalar scale, tol, j

    symeigensystem((V+V')/2,E,lambda)
    if (hasmissing(E) | hasmissing(lambda)) {
        jw_core_error("eigendecomposition failed for Gaussian simulation")
    }
    scale = max((max(abs(lambda)),1e-300))
    tol = 1e-10*scale
    if (min(lambda)<-tol) {
        jw_core_error("covariance is not positive semidefinite")
    }
    for (j=1; j<=cols(lambda); j++) {
        if (lambda[j]<0) lambda[j] = 0
    }
    /* symeigensystem() puts eigenvectors in COLUMNS. */
    return(E*diag(sqrt(lambda)))
}

real matrix jw_infer(struct jw_model scalar M, struct jw_fit scalar F,
                     string scalar selector, real scalar level,
                     real scalar gamma, real scalar mcerr, real scalar reps)
{
    real matrix result, S, actions, targets, B, Q, E, X, meta
    real rowvector active, hcols, eig, shift, support, tsum
    real colvector idx, groups, gj, blockidx, lower, upper, lo, hi, gl, gu
    real colvector ls, us
    real scalar o, i, j, g, start, stop, nb, cstart, cstop, first
    real scalar nout, m, estimate, variance, alpha, alpha2, dkw, plo, phi
    real scalar ilo, ihi, radius, pv, nactive, scale, tol, initialized

    if (missing(level) | level<=0 | level>=100) {
        jw_core_error("level() must be strictly between 0 and 100")
    }
    nout = rows(F.A)
    m = rows(M.z)
    if (cols(F.A)!=m | rows(F.targets)!=nout | cols(F.targets)!=m) {
        jw_core_error("incompatible output dimensions in inference")
    }
    result = J(nout,5,.)
    for (o=1; o<=nout; o++) {
        estimate = F.A[o,.]*M.z
        variance = jw_fixed_variance(M.V,F.A[o,.])
        result[o,.] = jw_infer_normal(estimate,variance,level)
    }
    if (selector!="focused") return(result)

    if (rows(F.allA)!=rows(F.allT) | rows(F.allA)!=rows(F.allmeta) |
        cols(F.allA)!=m | cols(F.allT)!=m | cols(F.allmeta)!=2) {
        jw_core_error("focused inference requires the complete candidate envelope")
    }
    if (missing(gamma) | gamma<=0 | gamma>=1 | missing(mcerr) |
        mcerr<=0 | mcerr>=1 | missing(reps) | reps<1 | reps!=floor(reps)) {
        jw_core_error("invalid gamma(), mcerr(), or reps()")
    }
    alpha = 1-level/100
    alpha2 = alpha-gamma-mcerr
    if (alpha2<=0) {
        jw_core_error("gamma()+mcerr() must be smaller than 1-level()/100")
    }
    dkw = sqrt(ln(4/mcerr)/(2*reps))
    plo = alpha2/2-dkw
    phi = 1-alpha2/2+dkw
    ilo = floor(reps*plo)
    ihi = ceil(reps*phi)
    initialized = 0

    for (o=1; o<=nout; o++) {
        idx = selectindex(F.allmeta[,1]:==o)
        if (rows(idx)==0) {
            jw_core_error("an output has no candidates for focused inference")
        }
        actions = F.allA[idx,.]
        targets = F.allT[idx,.]
        meta = F.allmeta[idx,.]
        B = actions-targets
        tol = 1e-10*max((1,max(abs(actions)),max(abs(targets))))
        if (max(abs(B[,1..M.p]))>tol) {
            jw_core_error("candidate adjustments may only use preperiod contrasts")
        }

        /* All candidates of a component must target the same quantity;
           the components must sum to the output's declared target. */
        groups = uniqrows(sort(meta[,2],1))
        tsum = J(1,m,0)
        for (g=1; g<=rows(groups); g++) {
            gj = selectindex(meta[,2]:==groups[g])
            if (max(abs(targets[gj,.]-J(rows(gj),1,1)*targets[gj[1],.]))>tol) {
                jw_core_error("candidate targets differ within an inference component")
            }
            tsum = tsum+targets[gj[1],.]
        }
        if (max(abs(tsum-F.targets[o,.]))>tol) {
            jw_core_error("inference components do not sum to the output target")
        }
        if (M.k==0) continue
        hcols = (M.p+1)..m
        /* Do not threshold small, genuine adjustment coefficients. */
        active = selectindex(colsum(abs(B[,hcols])):>0)
        nactive = cols(active)
        if (nactive==0) continue

        if (plo<=0 | phi>=1 | ilo<1 | ihi>reps) {
            jw_core_error("reps() is too small for the requested level(), gamma(), and mcerr(); increase reps()")
        }
        hcols = hcols[active]
        Q = M.V[hcols,hcols]
        eig = symeigenvalues((Q+Q')/2)
        scale = max((max(abs(eig)),1e-300))
        if (hasmissing(eig) | min(eig)<=1e-12*scale) {
            jw_core_error("active precontrast covariance is singular or numerically ill conditioned")
        }
        radius = sqrt(invchi2(nactive,1-gamma))
        B = B[,hcols]
        shift = (B*M.z[hcols,1])'
        support = J(1,rows(B),0)
        for (i=1; i<=rows(B); i++) {
            support[i] = radius*sqrt(jw_fixed_variance(Q,B[i,.]))
        }
        if (!initialized) {
            S = jw_infer_sqrt(M.V)
            initialized = 1
        }

        lower = J(reps,1,.)
        upper = J(reps,1,.)
        /* Bounded temporary storage: 512 draws and 64 actions at a time.
           All components within an output share E, not independent draws. */
        for (start=1; start<=reps; start=start+512) {
            stop = min((reps,start+511))
            nb = stop-start+1
            E = rnormal(nb,m,0,1)*S'
            lo = J(nb,1,0)
            hi = J(nb,1,0)
            for (g=1; g<=rows(groups); g++) {
                gj = selectindex(meta[,2]:==groups[g])
                first = 1
                for (cstart=1; cstart<=rows(gj); cstart=cstart+64) {
                    cstop = min((rows(gj),cstart+63))
                    blockidx = gj[|cstart\cstop|]
                    X = E*actions[blockidx,.]'
                    ls = rowmin(X+J(nb,1,1)*(shift[1,blockidx']-support[1,blockidx']))
                    us = rowmax(X+J(nb,1,1)*(shift[1,blockidx']+support[1,blockidx']))
                    if (first) {
                        gl = ls
                        gu = us
                        first = 0
                    }
                    else {
                        gl = rowmin((gl,ls))
                        gu = rowmax((gu,us))
                    }
                }
                lo = lo+gl
                hi = hi+gu
            }
            lower[|start\stop|] = lo
            upper[|start\stop|] = hi
        }
        estimate = result[o,1]
        /* These equal-tail p-values invert the same conservative tails.
           The null tested here is zero for each reported output. */
        pv = gamma+mcerr+2*min((sum(lower:<=estimate)/reps+dkw,
                               sum(upper:>=estimate)/reps+dkw))
        lower = sort(lower,1)
        upper = sort(upper,1)
        result[o,3] = estimate-upper[ihi]
        result[o,4] = estimate-lower[ilo]
        result[o,5] = min((1,pv))
    }
    return(result)
}
