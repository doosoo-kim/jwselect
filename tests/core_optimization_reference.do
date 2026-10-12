*! Frozen jwselect0.2 numerical/candidate/selection routines for regression tests.
*! Test-only reference; not loaded by the command or shipped in runtime bundles.
*! Retains original append order and separate factorization paths.
version 16.1
capture mata: mata drop jw020_mask_exists()
capture mata: mata drop jw020_candidates()
capture mata: mata drop jw020_core_error()
capture mata: mata drop jw020_loading()
capture mata: mata drop jw020_wald()
capture mata: mata drop jw020_fixed_variance()
capture mata: mata drop jw020_candidate_stats()
capture mata: mata drop jw020_focused_score()
capture mata: mata drop jw020_unique_masks()
capture mata: mata drop jw020_joint_masks()
capture mata: mata drop jw020_run()
mata:
/* Local candidate families. Their masks always index the canonical H vector. */
real scalar jw020_mask_exists(real matrix masks, real rowvector mask)
{
    real scalar j
    if (cols(masks)!=cols(mask)) return(0)
    for (j=1; j<=rows(masks); j++) {
        if (sum(abs(masks[j,.]-mask))==0) return(1)
    }
    return(0)
}

struct jw_family scalar jw020_candidates(struct jw_model scalar M,
    string scalar selectmethod, string scalar candidates, real matrix prematrix,
    real matrix candmatrix, real scalar maxcand)
{
    struct jw_family scalar F
    real matrix localm
    real colvector hi, ci, ip, ord
    real rowvector mask
    real scalar G, i, j, s, g, h, kg, nc, bit, u, haspre, hascand

    if (selectmethod!="manual" & selectmethod!="wbic" & selectmethod!="focused")
        jw_design_error("select() must be manual, wbic, or focused")
    if (candidates!="nested" & candidates!="all")
        jw_design_error("candidates() must be nested or all")
    if (missing(maxcand) | maxcand<1 | maxcand!=floor(maxcand))
        jw_design_error("maxcand() must be a positive integer")
    haspre = (rows(prematrix)>0 | cols(prematrix)>0)
    hascand = (rows(candmatrix)>0 | cols(candmatrix)>0)
    if (haspre & selectmethod!="manual")
        jw_design_error("prematrix() is available only with select(manual)")
    if (hascand & selectmethod=="manual")
        jw_design_error("candmatrix() is not available with select(manual)")
    if (haspre) {
        if (cols(prematrix)!=2 | hasmissing(prematrix))
            jw_design_error("prematrix() must have nonmissing columns cohort time")
        for (i=1; i<=rows(prematrix); i++) {
            h = jw_findcell(M.hmap,prematrix[i,1],prematrix[i,2],M.gap)
            if (h==0)
                jw_design_error("prematrix() contains an ineligible date (do not list the baseline)")
            for (j=1; j<i; j++) {
                if (prematrix[i,1]==prematrix[j,1] &
                    abs(prematrix[i,2]-prematrix[j,2])<=1e-8*M.gap)
                    jw_design_error("prematrix() contains duplicate dates")
            }
        }
    }
    if (hascand) {
        if (cols(candmatrix)!=M.k+1 | hasmissing(candmatrix))
            jw_design_error("candmatrix() must contain cohort followed by one bit per eligible preperiod")
        for (j=1; j<=rows(candmatrix); j++) {
            if (!any(M.groups:==candmatrix[j,1]))
                jw_design_error("candmatrix() contains an unknown cohort")
            if (M.k) {
                mask = candmatrix[j,2..(M.k+1)]
                if (any((mask:!=0):&(mask:!=1)))
                    jw_design_error("candmatrix() masks must be binary")
                if (sum(mask:*(M.hmap[,1]':!=candmatrix[j,1]))>0)
                    jw_design_error("a local candidate may contain only its own cohort's preperiods")
            }
        }
    }
    F.masks = J(0,M.k,0)
    F.owners = J(0,1,.)
    F.ids = J(0,1,.)
    G = rows(M.groups)
    for (i=1; i<=G; i++) {
        g = M.groups[i]
        /* Scalar predicates can yield a 1x0 index; enforce column shape. */
        hi = vec(selectindex(M.hmap[,1]:==g))
        kg = rows(hi)
        mask = J(1,M.k,0)
        localm = J(1,M.k,0)
        if (selectmethod=="manual") {
            if (haspre) {
                ip = vec(selectindex(prematrix[,1]:==g))
                for (j=1; j<=rows(ip); j++) {
                    h = jw_findcell(M.hmap,g,prematrix[ip[j],2],M.gap)
                    mask[h] = 1
                }
            }
            localm = mask
        }
        else if (hascand) {
            ci = vec(selectindex(candmatrix[,1]:==g))
            /* Absence of rows for a cohort is an error, not an implicit restriction. */
            if (rows(ci)==0)
                jw_design_error(sprintf("candmatrix() has no rows for cohort %g",g))
            for (j=1; j<=rows(ci); j++) {
                mask = J(1,M.k,0)
                if (M.k) mask = candmatrix[ci[j],2..(M.k+1)]
                /* A supplied baseline is coalesced with the mandatory first row. */
                if (sum(mask)==0) continue
                if (jw020_mask_exists(localm,mask))
                    jw_design_error("candmatrix() contains duplicate nonbaseline candidates")
                localm = localm\mask
                if (rows(localm)>maxcand)
                    jw_design_error("local candidate count exceeds maxcand()")
            }
        }
        else if (candidates=="nested") {
            if (kg+1>maxcand)
                jw_design_error("local candidate count exceeds maxcand()")
            if (kg) {
                ord = order(-M.hmap[hi,2],1)
                hi = hi[ord]
                for (j=1; j<=kg; j++) {
                    mask[hi[j]] = 1
                    localm = localm\mask
                }
            }
        }
        else {
            if (kg>floor(ln(maxcand)/ln(2)+1e-10))
                jw_design_error("all-subsets candidate count exceeds maxcand()")
            nc = 2^kg
            if (nc>maxcand)
                jw_design_error("all-subsets candidate count exceeds maxcand()")
            if (kg) {
                ord = order(-M.hmap[hi,2],1)
                hi = hi[ord]
                for (s=1; s<nc; s++) {
                    mask = J(1,M.k,0)
                    for (j=1; j<=kg; j++) {
                        bit = mod(floor(s/2^(j-1)),2)
                        mask[hi[j]] = bit
                    }
                    localm = localm\mask
                }
            }
        }
        F.masks = F.masks\localm
        F.owners = F.owners\J(rows(localm),1,g)
        F.ids = F.ids\(1::rows(localm))
    }
    return(F)
}

/* Numerical primitives. z=(post cells, optional pretrend contrasts).
   All covariance matrices are ordinary sampling covariance matrices.
   No Stata state is read or written by these routines. */

void jw020_core_error(string scalar message)
{
    errprintf("jwselect: %s\n", message)
    _error(3498)
}

/* A candidate estimator is represented by one row A: estimate=A*z.
   target is the baseline target row. Only the chosen h columns are
   adjusted. This retains every covariance, including across cohorts. */
real rowvector jw020_loading(real matrix V, real rowvector target,
                          real rowvector mask, real scalar k)
{
    real scalar p
    real rowvector j, A
    real matrix Q
    real colvector c, a

    p = cols(V)-k
    if (rows(V)!=cols(V) | cols(target)!=cols(V) | cols(mask)!=k) {
        jw020_core_error("incompatible dimensions in candidate loading")
    }
    A = target
    if (k==0) return(A)
    if (sum(mask)==0) return(A)
    j = selectindex(mask:==1):+p
    Q = V[j,j]
    c = V[j,.]*target'
    a = cholsolve(Q,c)
    if (hasmissing(a)) {
        jw020_core_error("candidate restriction covariance is singular or not positive definite")
    }
    A[1,j] = A[1,j]-a'
    return(A)
}

real scalar jw020_wald(real matrix V, real colvector z,
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
        jw020_core_error("candidate restriction covariance is singular or not positive definite")
    }
    W = h'*qh
    tol = 1e-10*max((1,abs(W)))
    if (W < -tol | missing(W)) jw020_core_error("invalid candidate Wald statistic")
    return(max((0,W)))
}

/* General fixed-action variance; unlike a v0-gain shortcut, this also
   works for weighted sums of cohort-specific adjustments. */
real scalar jw020_fixed_variance(real matrix V, real rowvector A)
{
    real scalar v, scale
    v = A*V*A'
    scale = max((1e-300,sum(abs(A))^2*max(abs(diagonal(V)))))
    if (missing(v) | v < -1e-8*scale) {
        jw020_core_error("negative or missing fixed-candidate variance")
    }
    return(max((0,v)))
}

/* Columns: r, Wald, WBIC, estimate, fixed-model SE, adjustment,
   covariance(baseline,adjustment), focused score. */
real rowvector jw020_candidate_stats(struct jw_model scalar M,
                                 real rowvector target,
                                 real rowvector mask)
{
    real rowvector A, R
    real scalar r, W, wbic, estimate, variance, delta, gain, v0, fmsc

    A = jw020_loading(M.V,target,mask,M.k)
    R = target-A
    r = sum(mask)
    W = jw020_wald(M.V,M.z,mask,M.k)
    wbic = .
    if (M.nunit>1 & M.nunit<.) wbic = W-r*ln(M.nunit)
    estimate = A*M.z
    variance = jw020_fixed_variance(M.V,A)
    delta = R*M.z
    gain = target*M.V*R'
    v0 = jw020_fixed_variance(M.V,target)
    fmsc = delta^2+v0-2*gain
    if (missing(estimate) | missing(delta) | missing(gain) | missing(fmsc)) {
        jw020_core_error("missing candidate calculation")
    }
    return((r,W,wbic,estimate,sqrt(variance),delta,gain,fmsc))
}

/* Reusable score for an arbitrary fixed linear adjustment. In particular,
   an aggregate of separate cohort adjustments need not equal joint MD. */
real scalar jw020_focused_score(real colvector z, real matrix V,
                            real rowvector target, real rowvector A)
{
    real rowvector R
    real scalar delta
    R = target-A
    delta = R*z
    return(delta^2+target*V*target'-2*target*V*R')
}

/* Candidate selection and linear maps for every displayed estimate.
   Stable family row order supplies the tie break: the first exact minimum.
   No outcome-specific selection code is duplicated inside inference. */

real matrix jw020_unique_masks(real matrix C, real scalar k,
                            real scalar maxjoint, real scalar addempty)
{
    real matrix U
    real scalar i, j, seen, hasempty
    U = J(0,k,0)
    hasempty = 0
    for (i=1; i<=rows(C); i++) {
        seen = 0
        for (j=1; j<=rows(U); j++) {
            if (sum(C[i,.]:!=U[j,.])==0) {
                seen = 1
                break
            }
        }
        if (!seen) {
            if (rows(U)>=maxjoint) jw020_core_error("joint candidate family exceeds maxjoint")
            U = U \ C[i,.]
            if (sum(C[i,.])==0) hasempty = 1
        }
    }
    if (addempty & !hasempty) {
        if (rows(U)>=maxjoint) jw020_core_error("joint candidate family including baseline exceeds maxjoint")
        U = J(1,k,0) \ U
    }
    return(U)
}

real matrix jw020_joint_masks(struct jw_model scalar M,
                          struct jw_family scalar F,
                          real matrix supplied, real scalar maxjoint)
{
    real matrix C, next, localset
    real colvector ii
    real scalar g, i, j, nr, nnew, pos

    if (maxjoint<1 | maxjoint>=. | maxjoint!=floor(maxjoint)) {
        jw020_core_error("maxjoint must be a positive finite integer")
    }
    if (rows(supplied)>0) {
        if (cols(supplied)!=M.k) jw020_core_error("joint candidate matrix has the wrong number of columns")
        if (hasmissing(supplied) | sum((supplied:!=0):&(supplied:!=1))>0) {
            jw020_core_error("joint candidate matrix must be binary")
        }
        return(jw020_unique_masks(supplied,M.k,maxjoint,1))
    }
    C = J(1,M.k,0)
    for (g=1; g<=rows(M.groups); g++) {
        ii = selectindex(F.owners:==M.groups[g])
        if (rows(ii)==0) jw020_core_error("a cohort has no candidates")
        localset = jw020_unique_masks(F.masks[ii,.],M.k,maxjoint,0)
        nr = rows(localset)
        if (rows(C)>floor(maxjoint/nr)) {
            jw020_core_error("Cartesian joint candidate family exceeds maxjoint; supply a smaller prespecified family")
        }
        nnew = rows(C)*nr
        next = J(nnew,M.k,0)
        pos = 0
        for (i=1; i<=rows(C); i++) {
            for (j=1; j<=nr; j++) {
                pos++
                next[pos,.] = (C[i,.]+localset[j,.]):>0
            }
        }
        C = next
    }
    return(jw020_unique_masks(C,M.k,maxjoint,1))
}

struct jw_fit scalar jw020_run(struct jw_model scalar M,
                            struct jw_family scalar F,
                            real matrix jointmatrix,
                            string scalar selectstr,
                            string scalar focusstr,
                            real scalar maxjoint)
{
    struct jw_fit scalar out
    real scalar G, m, g, i, j, r, oi, bestrow, bestscore, score,
                isfocused, docohort, doaggregate, primaryjoint, aggindex
    real colvector ii, localchosen
    real rowvector stat, target, mask, chosenunion, A, U, unitrow,
                   aggregateA, aggregatetarget
    real matrix C, localmasks

    G = rows(M.groups)
    m = rows(M.z)
    if (G<1 | M.p<1 | M.p+M.k!=m) jw020_core_error("invalid model dimensions")
    if (rows(M.V)!=m | cols(M.V)!=m | rows(M.T)!=G | cols(M.T)!=m) {
        jw020_core_error("incompatible model covariance or target dimensions")
    }
    if (rows(M.postmap)!=M.p | cols(M.postmap)!=2 | rows(M.hmap)!=M.k | cols(M.hmap)!=2 | cols(M.w)!=M.p) {
        jw020_core_error("incompatible cell maps or aggregation weights")
    }
    if (hasmissing(M.z) | hasmissing(M.V) | hasmissing(M.T) | hasmissing(M.w)) {
        jw020_core_error("model estimates, covariance, targets, and weights must be finite")
    }
    for (r=1; r<=M.p; r++) {
        if (sum(M.groups:==M.postmap[r,1])!=1) jw020_core_error("post-treatment cell has no unique model cohort")
    }
    if (cols(F.masks)!=M.k | rows(F.masks)!=rows(F.owners) | rows(F.ids)!=rows(F.owners)) {
        jw020_core_error("incompatible candidate-family dimensions")
    }
    if (hasmissing(F.masks) | sum((F.masks:!=0):&(F.masks:!=1))>0) {
        jw020_core_error("candidate masks must be binary")
    }
    if (selectstr!="manual" & selectstr!="wbic" & selectstr!="focused") {
        jw020_core_error("unknown selection method")
    }
    if (focusstr!="cohort" & focusstr!="aggregate" & focusstr!="both") {
        jw020_core_error("unknown focused-selection scope")
    }
    if (selectstr=="wbic" & (M.nunit<=1 | M.nunit>=.)) {
        jw020_core_error("Wald-BIC requires a finite independent-unit count greater than one")
    }
    for (i=1; i<=rows(F.owners); i++) {
        if (sum(M.groups:==F.owners[i])!=1) jw020_core_error("candidate owner is not a unique model cohort")
        if (M.k>0) {
            if (sum(F.masks[i,.]:*(M.hmap[.,1]':!=F.owners[i]))>0) {
                jw020_core_error("a cohort candidate contains another cohort's restrictions")
            }
        }
    }

    out.A = J(0,m,.)
    out.targets = J(0,m,.)
    out.selected = J(0,M.k,.)
    out.scores = J(0,11,.)
    out.cellA = J(M.p,m,0)
    out.cell_selected = J(M.p,M.k,0)
    out.allA = J(0,m,.)
    out.allT = J(0,m,.)
    out.allmeta = J(0,2,.)
    out.selA = J(0,m,.)
    out.selT = J(0,m,.)
    out.outmap = J(0,2,.)
    localchosen = J(G,1,.)
    localmasks = J(G,M.k,0)
    chosenunion = J(1,M.k,0)
    aggregatetarget = (M.w,J(1,M.k,0))
    isfocused = selectstr=="focused"
    docohort = isfocused & (focusstr=="cohort" | focusstr=="both")
    doaggregate = isfocused & (focusstr=="aggregate" | focusstr=="both")
    primaryjoint = !isfocused | (focusstr=="aggregate")

    /* Local scores are useful diagnostics in all modes. In manual/WBIC,
       their candidate estimates use own-cohort restrictions; the final
       reported estimator subsequently uses the selected joint system. */
    for (g=1; g<=G; g++) {
        ii = selectindex(F.owners:==M.groups[g])
        if (rows(ii)==0) jw020_core_error("a cohort has no candidates")
        if (selectstr=="manual" & rows(ii)!=1) jw020_core_error("manual mode requires exactly one candidate per cohort")
        bestrow = ii[1]
        bestscore = .
        for (j=1; j<=rows(ii); j++) {
            i = ii[j]
            stat = jw020_candidate_stats(M,M.T[g,.],F.masks[i,.])
            out.scores = out.scores \ (1,M.groups[g],F.ids[i],stat)
            if (selectstr=="wbic") score = stat[3]
            else score = stat[8]
            if (score<bestscore) {
                bestrow = i
                bestscore = score
            }
        }
        localchosen[g] = bestrow
        localmasks[g,.] = F.masks[bestrow,.]
        chosenunion = (chosenunion+localmasks[g,.]):>0
    }

    /* Aggregate focus evaluates the aggregate target itself for each
       allowed JOINT restriction configuration. No greedy substitution. */
    C = J(0,M.k,.)
    if (doaggregate) {
        C = jw020_joint_masks(M,F,jointmatrix,maxjoint)
        bestrow = 1
        bestscore = .
        for (i=1; i<=rows(C); i++) {
            stat = jw020_candidate_stats(M,aggregatetarget,C[i,.])
            out.scores = out.scores \ (2,0,i,stat)
            if (stat[8]<bestscore) {
                bestrow = i
                bestscore = stat[8]
            }
        }
        mask = C[bestrow,.]
        if (focusstr=="aggregate") chosenunion = mask
    }

    /* Primary cohort estimates and the complete underlying cell vector. */
    for (g=1; g<=G; g++) {
        if (primaryjoint) U = chosenunion
        else U = localmasks[g,.]
        A = jw020_loading(M.V,M.T[g,.],U,M.k)
        out.A = out.A \ A
        out.targets = out.targets \ M.T[g,.]
        out.selected = out.selected \ U
        if (isfocused & focusstr=="aggregate") out.outmap = out.outmap \ (5,M.groups[g])
        else out.outmap = out.outmap \ (1,M.groups[g])
        for (r=1; r<=M.p; r++) {
            if (M.postmap[r,1]==M.groups[g]) {
                unitrow = J(1,m,0)
                unitrow[r] = 1
                out.cellA[r,.] = jw020_loading(M.V,unitrow,U,M.k)
                out.cell_selected[r,.] = U
            }
        }
        /* Families of possible loadings for selection-safe inference. */
        if (docohort) {
            ii = selectindex(F.owners:==M.groups[g])
            for (j=1; j<=rows(ii); j++) {
                i = ii[j]
                out.allA = out.allA \ jw020_loading(M.V,M.T[g,.],F.masks[i,.],M.k)
                out.allT = out.allT \ M.T[g,.]
                out.allmeta = out.allmeta \ (g,1)
                out.selA = out.selA \ jw020_loading(M.V,M.T[g,.],F.masks[i,.],M.k)
                out.selT = out.selT \ M.T[g,.]
            }
        }
        else if (doaggregate) {
            for (i=1; i<=rows(C); i++) {
                out.allA = out.allA \ jw020_loading(M.V,M.T[g,.],C[i,.],M.k)
                out.allT = out.allT \ M.T[g,.]
                out.allmeta = out.allmeta \ (g,1)
                out.selA = out.selA \ jw020_loading(M.V,aggregatetarget,C[i,.],M.k)
                out.selT = out.selT \ aggregatetarget
            }
        }
        else {
            out.allA = out.allA \ A
            out.allT = out.allT \ M.T[g,.]
            out.allmeta = out.allmeta \ (g,1)
            out.selA = out.selA \ A
            out.selT = out.selT \ M.T[g,.]
        }
    }

    /* This is exactly the weighted sum of primary cell estimates. With
       arbitrary cell weights it need not be a weighted sum of uniformly
       averaged cohort summaries; the displayed target metadata says so. */
    aggregateA = M.w*out.cellA
    aggindex = rows(out.A)+1
    out.A = out.A \ aggregateA
    out.targets = out.targets \ aggregatetarget
    out.selected = out.selected \ chosenunion
    if (docohort) out.outmap = out.outmap \ (2,0)
    else if (doaggregate) out.outmap = out.outmap \ (3,0)
    else out.outmap = out.outmap \ (4,0)

    if (docohort) {
        /* For the weighted aggregate of cohort-focused estimators, keep
           each cohort's reporting and selection families separately.
           Selection uses the cohort summary target even when custom cell
           weights make the reporting component a different target. */
        for (g=1; g<=G; g++) {
            target = J(1,m,0)
            for (r=1; r<=M.p; r++) {
                if (M.postmap[r,1]==M.groups[g]) target[r] = M.w[r]
            }
            ii = selectindex(F.owners:==M.groups[g])
            for (j=1; j<=rows(ii); j++) {
                i = ii[j]
                out.allA = out.allA \ jw020_loading(M.V,target,F.masks[i,.],M.k)
                out.allT = out.allT \ target
                out.allmeta = out.allmeta \ (aggindex,g)
                out.selA = out.selA \ jw020_loading(M.V,M.T[g,.],F.masks[i,.],M.k)
                out.selT = out.selT \ M.T[g,.]
            }
        }
    }
    else if (doaggregate) {
        for (i=1; i<=rows(C); i++) {
            out.allA = out.allA \ jw020_loading(M.V,aggregatetarget,C[i,.],M.k)
            out.allT = out.allT \ aggregatetarget
            out.allmeta = out.allmeta \ (aggindex,1)
            out.selA = out.selA \ jw020_loading(M.V,aggregatetarget,C[i,.],M.k)
            out.selT = out.selT \ aggregatetarget
        }
    }
    else {
        out.allA = out.allA \ aggregateA
        out.allT = out.allT \ aggregatetarget
        out.allmeta = out.allmeta \ (aggindex,1)
        out.selA = out.selA \ aggregateA
        out.selT = out.selT \ aggregatetarget
    }

    /* In both mode preserve the primary cohort-wise fit, and append the
       independently selected aggregate-focused estimate and its family. */
    if (isfocused & focusstr=="both") {
        oi = rows(out.A)+1
        out.A = out.A \ jw020_loading(M.V,aggregatetarget,mask,M.k)
        out.targets = out.targets \ aggregatetarget
        out.selected = out.selected \ mask
        out.outmap = out.outmap \ (3,0)
        for (i=1; i<=rows(C); i++) {
            out.allA = out.allA \ jw020_loading(M.V,aggregatetarget,C[i,.],M.k)
            out.allT = out.allT \ aggregatetarget
            out.allmeta = out.allmeta \ (oi,1)
            out.selA = out.selA \ jw020_loading(M.V,aggregatetarget,C[i,.],M.k)
            out.selT = out.selT \ aggregatetarget
        }
    }
    return(out)
}

end
