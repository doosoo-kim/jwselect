/* Numerical primitives. z=(post cells, optional pretrend contrasts).
   All covariance matrices are ordinary sampling covariance matrices.
   No Stata state is read or written by these routines. */

void jw_core_error(string scalar message)
{
    errprintf("jwselect: %s\n", message)
    _error(3498)
}

/* A candidate estimator is represented by one row A: estimate=A*z.
   target is the baseline target row. Only the chosen h columns are
   adjusted. This retains every covariance, including across cohorts. */
/* Several reporting targets sharing a mask use one covariance factorization.
   cholsolve accepts multiple right-hand sides; no inverse is formed. */
real matrix jw_loadings(real matrix V, real matrix targets,
                        real rowvector mask, real scalar k)
{
    real scalar p
    real rowvector j
    real matrix A, a

    p = cols(V)-k
    if (rows(V)!=cols(V) | cols(targets)!=cols(V) | cols(mask)!=k) {
        jw_core_error("incompatible dimensions in candidate loading")
    }
    A = targets
    if (k==0) return(A)
    if (sum(mask)==0) return(A)
    j = selectindex(mask:==1):+p
    a = cholsolve(V[j,j],V[j,.]*targets')
    if (hasmissing(a)) {
        jw_core_error("candidate restriction covariance is singular or not positive definite")
    }
    A[.,j] = A[.,j]-a'
    return(A)
}

real rowvector jw_loading(real matrix V, real rowvector target,
                          real rowvector mask, real scalar k)
{
    return(jw_loadings(V,target,mask,k))
}

real scalar jw_wald(real matrix V, real colvector z,
                    real rowvector mask, real scalar k)
{
    real scalar p, W, tol
    real rowvector j
    real colvector h, qh

    if (k==0) return(0)
    if (sum(mask)==0) return(0)
    p = rows(z)-k
    j = selectindex(mask:==1):+p
    h = z[j,1]
    qh = cholsolve(V[j,j],h)
    if (hasmissing(qh)) {
        jw_core_error("candidate restriction covariance is singular or not positive definite")
    }
    W = h'*qh
    tol = 1e-10*max((1,abs(W)))
    if (W < -tol | missing(W)) jw_core_error("invalid candidate Wald statistic")
    return(max((0,W)))
}

/* General fixed-action variance; unlike a v0-gain shortcut, this also
   works for weighted sums of cohort-specific adjustments. */
real scalar jw_fixed_variance(real matrix V, real rowvector A)
{
    real scalar v, scale
    v = A*V*A'
    scale = max((1e-300,sum(abs(A))^2*max(abs(diagonal(V)))))
    if (missing(v) | v < -1e-8*scale) {
        jw_core_error("negative or missing fixed-candidate variance")
    }
    return(max((0,v)))
}

/* Columns: r, Wald, WBIC, estimate, fixed-model SE, adjustment,
   covariance(baseline,adjustment), focused score. */
/* A is also returned by reference so selection and inference families can
   reuse the identical loading. The Wald solve and loading solve share one
   factorization of the candidate restriction covariance. */
real rowvector jw_candidate_stats_loading(struct jw_model scalar M,
                                 real rowvector target,
                                 real rowvector mask,
                                 real rowvector A)
{
    real rowvector R, j
    real colvector h, c
    real matrix solved
    real scalar r, W, wbic, estimate, variance, delta, gain, v0, fmsc, tol

    if (rows(M.V)!=cols(M.V) | cols(target)!=cols(M.V) | cols(mask)!=M.k) {
        jw_core_error("incompatible dimensions in candidate loading")
    }
    A = target
    r = sum(mask)
    W = 0
    if (r>0) {
        j = selectindex(mask:==1):+(rows(M.z)-M.k)
        h = M.z[j,1]
        c = M.V[j,.]*target'
        solved = cholsolve(M.V[j,j],(c,h))
        if (hasmissing(solved)) {
            jw_core_error("candidate restriction covariance is singular or not positive definite")
        }
        A[1,j] = A[1,j]-solved[.,1]'
        W = h'*solved[.,2]
        tol = 1e-10*max((1,abs(W)))
        if (W < -tol | missing(W)) jw_core_error("invalid candidate Wald statistic")
        W = max((0,W))
    }
    R = target-A
    wbic = .
    if (M.nunit>1 & M.nunit<.) wbic = W-r*ln(M.nunit)
    estimate = A*M.z
    variance = jw_fixed_variance(M.V,A)
    delta = R*M.z
    gain = target*M.V*R'
    v0 = jw_fixed_variance(M.V,target)
    fmsc = delta^2+v0-2*gain
    if (missing(estimate) | missing(delta) | missing(gain) | missing(fmsc)) {
        jw_core_error("missing candidate calculation")
    }
    return((r,W,wbic,estimate,sqrt(variance),delta,gain,fmsc))
}

/* Original public primitive retained for callers that only need scores. */
real rowvector jw_candidate_stats(struct jw_model scalar M,
                                 real rowvector target,
                                 real rowvector mask)
{
    real rowvector A
    return(jw_candidate_stats_loading(M,target,mask,A))
}

/* Reusable score for an arbitrary fixed linear adjustment. In particular,
   an aggregate of separate cohort adjustments need not equal joint MD. */
real scalar jw_focused_score(real colvector z, real matrix V,
                            real rowvector target, real rowvector A)
{
    real rowvector R
    real scalar delta
    R = target-A
    delta = R*z
    return(delta^2+target*V*target'-2*target*V*R')
}
