/* Reselection-based two-step focused inference, in sampling-variance units.

   At every trial nuisance mean mu and every common Gaussian draw E, select
   the first minimizer of the ORIGINAL focused criterion. Reporting uses
   allA/allT; selection uses selA/selT. These need not describe the same
   target. Additive components share the complete Gaussian draw.

   Numerical outer extrema are obtained by multistart projected coordinate
   pattern searches on a unit ball. This is NOT NOMAD/MADS and does not
   certify a global extremum. Independent validation controls Monte Carlo
   error at the FOUR located points only; it does not bound missed extrema.

   diagnostics columns:
   output objective nactive search_reps validation_reps evaluations
   starts_used converged_starts budget_exhausted search_value
   validated_value validation_gap radius_fraction
   Objectives: 1 lower quantile (min), 2 upper quantile (max),
               3 P(Y<=estimate) (max), 4 P(Y>=estimate) (max).
   extrema columns: output objective then all M.k nuisance means.
*/

struct jw_reselect_result {
    real matrix table
    real matrix diagnostics
    real matrix extrema
}

struct jw_rs_context {
    real matrix reportB
    real matrix selectB
    real colvector constant
    real colvector component
    real colvector groups
    real matrix EH
    real colvector baseline
    real matrix projection
}

struct jw_rs_search_result {
    real colvector x
    real scalar value
    real scalar evaluations
    real scalar starts
    real scalar converged
    real scalar exhausted
}

/* cachemb is a ceiling for the reusable search projection only, in decimal
   MB. Other working arrays require additional RAM. Validation passes zero. */
struct jw_rs_context scalar jw_rs_prepare_opts(struct jw_model scalar M,
                    struct jw_fit scalar F, real scalar o,
                    real matrix E, real rowvector hcols, real scalar cachemb)
{
    struct jw_rs_context scalar C
    real colvector ix, ci
    real matrix DA
    real rowvector target, twiceTV
    real scalar g, variance, ncache

    ix = selectindex(F.allmeta[,1]:==o)
    C.reportB = F.allA[ix,hcols]-F.allT[ix,hcols]
    C.selectB = F.selA[ix,hcols]-F.selT[ix,hcols]
    C.component = F.allmeta[ix,2]
    C.groups = uniqrows(sort(C.component,1))
    C.constant = J(rows(ix),1,.)
    DA = F.selT[ix,.]-F.selA[ix,.]
    /* Within a component the target is shared by every candidate. Reuse
       its covariance product instead of recomputing it for each subset. */
    for (g=1; g<=rows(C.groups); g++) {
        ci = selectindex(C.component:==C.groups[g])
        target = F.selT[ix[ci[1]],.]
        variance = jw_fixed_variance(M.V,target)
        twiceTV = (2*target)*M.V
        C.constant[ci,1] = J(rows(ci),1,variance)-DA[ci,.]*twiceTV'
    }
    C.EH = E[,hcols]
    C.baseline = E*F.targets[o,.]'
    C.projection = J(0,0,.)
    if (cachemb>0 & rows(E)>0) {
        ncache = min((rows(ix),floor(cachemb*1000000/(8*rows(E)))))
        if (ncache>0) C.projection = C.EH*C.selectB[1..ncache,.]'
    }
    return(C)
}

/* Source-compatible helper retained for tests and external Mata callers. */
struct jw_rs_context scalar jw_rs_prepare(struct jw_model scalar M,
                    struct jw_fit scalar F, real scalar o,
                    real matrix E, real rowvector hcols)
{
    return(jw_rs_prepare_opts(M,F,o,E,hcols,64))
}

void jw_rs_cache(struct jw_rs_context scalar C, real scalar cachemb)
{
    real scalar ncache
    C.projection = J(0,0,.)
    if (cachemb>0 & rows(C.EH)>0) {
        ncache = min((rows(C.selectB),floor(cachemb*1000000/(8*rows(C.EH)))))
        if (ncache>0) C.projection = C.EH*C.selectB[1..ncache,.]'
    }
}

real colvector jw_rs_draw_errors(struct jw_rs_context scalar C,
                                real colvector mu)
{
    real matrix W, Q, P
    real colvector y, ix, ci, changed, best, chosen, cached, uncached
    real rowvector offset
    real scalar g, start, stop, j, n, ncache

    n = rows(C.EH)
    ncache = cols(C.projection)
    W = C.EH+J(n,1,1)*mu'
    y = C.baseline
    for (g=1; g<=rows(C.groups); g++) {
        ix = selectindex(C.component:==C.groups[g])
        best = J(n,1,.)
        chosen = J(n,1,ix[1])
        for (start=1; start<=rows(ix); start=start+64) {
            stop = min((rows(ix),start+63))
            ci = ix[|start\stop|]
            if (max(ci)<=ncache) {
                offset = (C.selectB[ci,.]*mu)'
                P = C.projection[,ci']+J(n,1,1)*offset
            }
            else if (min(ci)>ncache) P = W*C.selectB[ci,.]'
            else {
                /* One block can straddle the cached prefix. Column order
                   stays unchanged, preserving first-candidate ties. */
                cached = selectindex(ci:<=ncache)
                uncached = selectindex(ci:>ncache)
                P = J(n,rows(ci),.)
                offset = (C.selectB[ci[cached],.]*mu)'
                P[,cached'] = C.projection[,ci[cached]']+J(n,1,1)*offset
                P[,uncached'] = W*C.selectB[ci[uncached],.]'
            }
            Q = P:^2+J(n,1,1)*C.constant[ci,1]'
            for (j=1; j<=rows(ci); j++) {
                /* Strict improvement preserves original family order. */
                changed = vec(selectindex(Q[,j]:<best))
                if (rows(changed)>0) {
                    best[changed,1] = Q[changed,j]
                    chosen[changed,1] = J(rows(changed),1,ci[j])
                }
            }
        }
        if (hasmissing(best)) jw_core_error("nonfinite simulated focused scores")
        y = y+rowsum(W:*C.reportB[chosen,.])
    }
    if (hasmissing(y)) jw_core_error("nonfinite simulated selected errors")
    return(y)
}

/* Deterministic, public verification helper. mu contains every optional
   precontrast mean, in model order; E is supplied, not generated here. */
real colvector jw_rs_errors(struct jw_model scalar M,
                            struct jw_fit scalar F, real scalar o,
                            real matrix E, real colvector mu)
{
    struct jw_rs_context scalar C
    real rowvector hcols
    if (cols(E)!=rows(M.z) | rows(mu)!=M.k | cols(mu)!=1) {
        jw_core_error("incompatible supplied draws or nuisance means")
    }
    if (M.k==0) return(E*F.targets[o,.]')
    hcols = (M.p+1)..rows(M.z)
    C = jw_rs_prepare(M,F,o,E,hcols)
    return(jw_rs_draw_errors(C,mu))
}

/* Empirical inverse-CDF convention, with no interpolated artificial draw. */
real rowvector jw_rs_statistics(real colvector y, real scalar tail,
                               real scalar estimate)
{
    real colvector ys
    real scalar n, il, iu
    n = rows(y)
    il = max((1,ceil(n*tail)))
    iu = min((n,ceil(n*(1-tail))))
    ys = sort(y,1)
    return((ys[il],ys[iu],sum(y:<=estimate)/n,sum(y:>=estimate)/n))
}

real scalar jw_rs_objective(struct jw_rs_context scalar C,
                            real colvector h, real matrix B,
                            real colvector x, real scalar tail,
                            real scalar estimate, real scalar objective)
{
    real colvector y, ys
    real scalar n, rank
    y = jw_rs_draw_errors(C,h+B*x)
    n = rows(y)
    /* Tail objectives require counts only, not an O(n log n) sort. */
    if (objective==3) return(sum(y:<=estimate)/n)
    if (objective==4) return(sum(y:>=estimate)/n)
    if (objective==1) rank = max((1,ceil(n*tail)))
    else rank = min((n,ceil(n*(1-tail))))
    ys = sort(y,1)
    return(ys[rank])
}

real colvector jw_rs_project(real colvector x)
{
    real scalar length
    length = sqrt(sum(x:^2))
    if (length>1) return(x/length)
    return(x)
}

/* The starts include the center, a boundary point toward the null mean,
   and alternating random boundary/interior points. They are generated
   once and reused for the four objective searches. */
real matrix jw_rs_starts(real colvector h, real matrix B,
                         real scalar nstarts)
{
    real matrix X
    real colvector v
    real scalar s, d, length
    d = rows(h)
    X = J(d,nstarts,0)
    for (s=2; s<=nstarts; s++) {
        if (s==2) v = qrsolve(B,-h)
        else v = rnormal(d,1,0,1)
        length = sqrt(sum(v:^2))
        if (hasmissing(v) | length==0) {
            v = J(d,1,0)
            v[1] = 1
            length = 1
        }
        v = v/length
        if (s>2 & mod(s,2)==0) v = .5*v
        X[,s] = v
    }
    return(X)
}

/* One independent starting-point path. Numerical order is the same as
   the original serial routine, including strict first-in-family ties.
   Progress counts objective evaluations, not candidate/draw comparisons. */
struct jw_rs_search_result scalar jw_rs_search_start(
                    struct jw_rs_context scalar C,
                    real colvector h, real matrix B, real colvector start,
                    real scalar tail, real scalar estimate,
                    real scalar objective, real scalar budget,
                    real scalar tol, real scalar progress)
{
    struct jw_rs_search_result scalar R
    real colvector x, trial, incumbent, direction
    real scalar j, sg, used, d, value, trialvalue,
                bestvalue, step, improved, factor, complete, length

    d = rows(h)
    factor = 1
    if (objective==1) factor = -1
    x = start
    value = jw_rs_objective(C,h,B,x,tail,estimate,objective)
    used = 1
    R.value = value
    R.x = x
    R.starts = 1
    R.converged = 0
    R.exhausted = 0
    step = .5
    while (step>tol & used<budget) {
        incumbent = x
        bestvalue = value
        improved = 0
        complete = 1
        for (j=1; j<=d+1; j++) {
            direction = J(d,1,0)
            if (j<=d) direction[j] = 1
            else {
                length = sqrt(sum(x:^2))
                if (length==0) continue
                direction = x/length
            }
            for (sg=-1; sg<=1; sg=sg+2) {
                if (used>=budget) {
                    complete = 0
                    break
                }
                trial = jw_rs_project(x+sg*step*direction)
                if (max(abs(trial-x))<1e-14) continue
                trialvalue = jw_rs_objective(C,h,B,trial,tail,estimate,objective)
                used++
                if (progress & mod(used,50)==0) {
                    printf("    evaluations %g/%g, mesh %g\n",used,budget,step)
                    displayflush()
                }
                if (factor*trialvalue>factor*bestvalue) {
                    incumbent = trial
                    bestvalue = trialvalue
                    improved = 1
                }
            }
            if (!complete) break
        }
        if (improved) {
            x = incumbent
            value = bestvalue
            if (factor*value>factor*R.value) {
                R.value = value
                R.x = x
            }
        }
        else if (complete) step = step/2
        if (!complete) break
    }
    R.evaluations = used
    if (step<=tol) R.converged = 1
    else R.exhausted = 1
    return(R)
}

real scalar jw_rs_start_budget(real scalar maxeval, real scalar nstarts,
                              real scalar s)
{
    return(floor(maxeval/nstarts)+(s<=mod(maxeval,nstarts)))
}

/* Merge in original start order; only strict improvement replaces a tie. */
void jw_rs_merge(struct jw_rs_search_result scalar R,
                 struct jw_rs_search_result scalar S, real scalar objective)
{
    real scalar factor
    factor = 1
    if (objective==1) factor = -1
    if (missing(R.value)) {
        R.value = S.value
        R.x = S.x
    }
    else if (factor*S.value>factor*R.value) {
        R.value = S.value
        R.x = S.x
    }
    R.evaluations = R.evaluations+S.evaluations
    R.starts = R.starts+S.starts
    R.converged = R.converged+S.converged
    R.exhausted = max((R.exhausted,S.exhausted))
}

struct jw_rs_search_result scalar jw_rs_search_opts(
                    struct jw_rs_context scalar C,
                    real colvector h, real matrix B, real matrix starts,
                    real scalar tail, real scalar estimate,
                    real scalar objective, real scalar maxeval,
                    real scalar tol, real scalar progress)
{
    struct jw_rs_search_result scalar R, S
    real scalar s, budget
    R.value = .
    R.x = starts[,1]
    R.evaluations = 0
    R.starts = 0
    R.converged = 0
    R.exhausted = 0
    for (s=1; s<=cols(starts); s++) {
        budget = jw_rs_start_budget(maxeval,cols(starts),s)
        if (progress) {
            printf("  objective %g, start %g/%g, budget %g\n",
                   objective,s,cols(starts),budget)
            displayflush()
        }
        S = jw_rs_search_start(C,h,B,starts[,s],tail,estimate,
                                objective,budget,tol,progress)
        jw_rs_merge(R,S,objective)
    }
    return(R)
}

struct jw_rs_search_result scalar jw_rs_search(
                    struct jw_rs_context scalar C,
                    real colvector h, real matrix B, real matrix starts,
                    real scalar tail, real scalar estimate,
                    real scalar objective, real scalar maxeval,
                    real scalar tol)
{
    return(jw_rs_search_opts(C,h,B,starts,tail,estimate,
                            objective,maxeval,tol,0))
}

/* Exact one-sided order-statistic bounds for a population p quantile.
   Lower rank k satisfies P(Bin(n,p)<k)<=delta; upper rank satisfies
   P(Bin(n,p)>=k)<=delta. Return 0/n+1 when finite bounds do not exist. */
real rowvector jw_rs_quantile_ranks(real scalar n, real scalar p,
                                   real scalar delta)
{
    real scalar lo, hi, mid, lower, upper
    lower = 0
    lo = 1
    hi = n
    while (lo<=hi) {
        mid = floor((lo+hi)/2)
        if (binomial(n,mid-1,p)<=delta) {
            lower = mid
            lo = mid+1
        }
        else hi = mid-1
    }
    upper = n+1
    lo = 1
    hi = n
    while (lo<=hi) {
        mid = floor((lo+hi)/2)
        if (binomialtail(n,mid,p)<=delta) {
            upper = mid
            hi = mid-1
        }
        else lo = mid+1
    }
    return((lower,upper))
}

/* Exact one-sided Clopper-Pearson upper probability bound. */
real scalar jw_rs_probability_upper(real scalar successes, real scalar n,
                                    real scalar delta)
{
    if (successes>=n) return(1)
    return(invibeta(successes+1,n-successes,1-delta))
}

/* Shared preparation/validation makes serial and distributed execution
   use identical numerical objects and RNG consumption order. */
struct jw_rs_output {
    real scalar output
    real scalar estimate
    real rowvector active
    real rowvector hcols
    real colvector fullmu
    real colvector h
    real matrix B
    real matrix starts
    struct jw_rs_context scalar C
    struct jw_rs_context scalar CV
}

real rowvector jw_rs_validate(struct jw_model scalar M,
                    struct jw_fit scalar F,
                    real scalar level, real scalar gamma,
                    real scalar mcerr, real scalar reps,
                    real scalar optstarts, real scalar optmaxeval,
                    real scalar opttol, real scalar optreps,
                    real rowvector outputs, real scalar cachemb)
{
    real scalar nout, m, alpha, tail, delta, il, iu
    real rowvector ranks
    nout = rows(F.A)
    m = rows(M.z)
    if (missing(level) | level<=0 | level>=100) {
        jw_core_error("level() must be strictly between 0 and 100")
    }
    if (cols(F.A)!=m | rows(F.targets)!=nout | cols(F.targets)!=m) {
        jw_core_error("incompatible output dimensions in reselection inference")
    }
    if (rows(F.allA)!=rows(F.allT) | rows(F.allA)!=rows(F.allmeta) |
        rows(F.selA)!=rows(F.allA) | rows(F.selT)!=rows(F.allA) |
        cols(F.allA)!=m | cols(F.allT)!=m | cols(F.selA)!=m |
        cols(F.selT)!=m | cols(F.allmeta)!=2 | hasmissing(F.allA) |
        hasmissing(F.allT) | hasmissing(F.selA) | hasmissing(F.selT)) {
        jw_core_error("reselection requires compatible reporting and selection candidate maps")
    }
    if (missing(gamma) | gamma<=0 | gamma>=1 | missing(mcerr) |
        mcerr<=0 | mcerr>=1 | missing(reps) | reps<1 |
        reps!=floor(reps) | missing(optreps) | optreps<100 |
        optreps!=floor(optreps)) {
        jw_core_error("invalid gamma(), mcerr(), reps(), or optreps()")
    }
    if (missing(optstarts) | optstarts<1 | optstarts!=floor(optstarts) |
        missing(optmaxeval) | optmaxeval<optstarts |
        optmaxeval!=floor(optmaxeval) | missing(opttol) |
        opttol<=0 | opttol>=.5) {
        jw_core_error("invalid optstarts(), optmaxeval(), or opttol()")
    }
    alpha = 1-level/100
    if (alpha-gamma-mcerr<=0) {
        jw_core_error("gamma()+mcerr() must be smaller than 1-level()/100")
    }
    tail = (alpha-gamma-mcerr)/2
    /* Four located points, four one-sided validation statements each.
       This is marginal per reported output, not simultaneous across outputs. */
    delta = mcerr/16
    ranks = jw_rs_quantile_ranks(reps,tail,delta)
    il = ranks[1]
    ranks = jw_rs_quantile_ranks(reps,1-tail,delta)
    iu = ranks[2]
    if (missing(cachemb) | cachemb<0) jw_core_error("cachemb() must be nonnegative")
    if (cols(outputs)>0) {
        if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout |
            any(outputs:!=floor(outputs))) jw_core_error("invalid inference output indices")
    }
    return((tail,delta,il,iu))
}

struct jw_reselect_result scalar jw_rs_result_init(
                    struct jw_model scalar M, struct jw_fit scalar F,
                    real scalar level)
{
    struct jw_reselect_result scalar R
    real scalar nout, o, j, r
    nout = rows(F.A)
    R.table = J(nout,5,.)
    R.diagnostics = J(4*nout,13,.)
    R.extrema = J(4*nout,2+M.k,.)
    for (o=1; o<=nout; o++) {
        R.table[o,.] = jw_infer_normal(F.A[o,.]*M.z,
                              jw_fixed_variance(M.V,F.A[o,.]),level)
        /* Only point estimate and fixed-candidate SE for unrequested rows. */
        R.table[o,3..5] = J(1,3,.)
        for (j=1; j<=4; j++) {
            r = 4*(o-1)+j
            R.diagnostics[r,1..2] = (o,j)
            R.extrema[r,1..2] = (o,j)
        }
    }
    return(R)
}

struct jw_rs_output scalar jw_rs_output_plan(
                    struct jw_model scalar M, struct jw_fit scalar F,
                    real scalar o, real scalar gamma)
{
    struct jw_rs_output scalar P
    real matrix Bfull, Dfull, Q
    real colvector ix, gi, groups
    real rowvector hc, tsum, eig
    real scalar m, d, g, tolerance, scale, radius
    m = rows(M.z)
    P.output = o
    P.estimate = F.A[o,.]*M.z
    P.active = J(1,0,.)
    P.hcols = J(1,0,.)
    P.fullmu = J(0,1,.)
    P.h = J(0,1,.)
    P.B = J(0,0,.)
    P.starts = J(0,0,.)
    hc = J(1,0,.)
    if (M.k>0) {
        hc = (M.p+1)..m
        P.fullmu = M.z[hc,1]
    }
    ix = selectindex(F.allmeta[,1]:==o)
    if (rows(ix)==0) jw_core_error("an output has no reselection candidates")
    Bfull = F.allA[ix,.]-F.allT[ix,.]
    Dfull = F.selA[ix,.]-F.selT[ix,.]
    tolerance = 1e-10*max((1,max(abs(F.allA[ix,.])),max(abs(F.selA[ix,.]))))
    if (max(abs(Bfull[,1..M.p]))>tolerance |
        max(abs(Dfull[,1..M.p]))>tolerance) {
        jw_core_error("focused candidate adjustments may only use precontrasts")
    }
    groups = uniqrows(sort(F.allmeta[ix,2],1))
    tsum = J(1,m,0)
    for (g=1; g<=rows(groups); g++) {
        gi = selectindex(F.allmeta[ix,2]:==groups[g])
        if (max(abs(F.allT[ix[gi],.]-J(rows(gi),1,1)*F.allT[ix[gi[1]],.]))>tolerance |
            max(abs(F.selT[ix[gi],.]-J(rows(gi),1,1)*F.selT[ix[gi[1]],.]))>tolerance) {
            jw_core_error("candidate targets differ within a reselection component")
        }
        tsum = tsum+F.allT[ix[gi[1]],.]
    }
    if (max(abs(tsum-F.targets[o,.]))>tolerance) {
        jw_core_error("reselection component targets do not sum to output target")
    }
    if (M.k>0) {
        if (max(abs(Bfull[,hc]))>0) {
            P.active = selectindex(colsum(abs(Bfull[,hc])+abs(Dfull[,hc])):>0)
        }
    }
    d = cols(P.active)
    if (d==0) return(P)
    P.hcols = hc[P.active]
    P.h = M.z[P.hcols,1]
    Q = M.V[P.hcols,P.hcols]
    eig = symeigenvalues((Q+Q')/2)
    scale = max((max(abs(eig)),1e-300))
    if (hasmissing(eig) | min(eig)<=1e-12*scale) {
        jw_core_error("active reselection precontrast covariance is singular or ill conditioned")
    }
    radius = sqrt(invchi2(d,1-gamma))
    P.B = radius*jw_infer_sqrt(Q)
    return(P)
}

void jw_rs_zero_output(struct jw_reselect_result scalar R,
                       struct jw_rs_output scalar P,
                       real scalar variance, real scalar level)
{
    real scalar j, r, k
    k = rows(P.fullmu)
    R.table[P.output,.] = jw_infer_normal(P.estimate,variance,level)
    for (j=1; j<=4; j++) {
        r = 4*(P.output-1)+j
        R.diagnostics[r,.] = (P.output,j,0,0,0,0,0,0,0,.,.,.,0)
        if (k>0) R.extrema[r,3..(2+k)] = P.fullmu'
    }
}

void jw_rs_record_search(struct jw_reselect_result scalar R,
                         struct jw_rs_output scalar P,
                         struct jw_rs_search_result scalar SR,
                         real scalar objective, real scalar optreps,
                         real scalar reps)
{
    real scalar r
    real colvector fullmu
    r = 4*(P.output-1)+objective
    R.diagnostics[r,.] = (P.output,objective,cols(P.active),optreps,reps,
            SR.evaluations,SR.starts,SR.converged,SR.exhausted,SR.value,
            .,.,sqrt(sum(SR.x:^2)))
    fullmu = P.fullmu
    fullmu[P.active,1] = P.h+P.B*SR.x
    R.extrema[r,.] = (P.output,objective,fullmu')
}

/* Every one of the four located points contributes all four bounds.
   Validation uncertainty is allocated mcerr/16 as in the original code. */
void jw_rs_finish_output(struct jw_reselect_result scalar R,
                         struct jw_rs_output scalar P,
                         real matrix located, real rowvector searchvalues,
                         real rowvector numerical, real scalar gamma,
                         real scalar mcerr, real scalar progress)
{
    real matrix bounds
    real colvector y, ys
    real rowvector v
    real scalar j, r, n, il, iu, tail, delta, lo, hi, pv
    tail = numerical[1]
    delta = numerical[2]
    il = numerical[3]
    iu = numerical[4]
    n = rows(P.CV.EH)
    bounds = J(4,4,.)
    for (j=1; j<=4; j++) {
        if (progress) {
            printf("  output %g: validating located point %g/4\n",P.output,j)
            displayflush()
        }
        y = jw_rs_draw_errors(P.CV,located[,j])
        ys = sort(y,1)
        bounds[j,.] = (ys[il],ys[iu],
            jw_rs_probability_upper(sum(y:<=P.estimate),n,delta),
            jw_rs_probability_upper(sum(y:>=P.estimate),n,delta))
        v = jw_rs_statistics(y,tail,P.estimate)
        r = 4*(P.output-1)+j
        R.diagnostics[r,11] = v[j]
        R.diagnostics[r,12] = abs(v[j]-searchvalues[j])
    }
    lo = min(bounds[,1])
    hi = max(bounds[,2])
    pv = min((1,gamma+mcerr+2*min((max(bounds[,3]),max(bounds[,4])))))
    R.table[P.output,3] = P.estimate-hi
    R.table[P.output,4] = P.estimate-lo
    R.table[P.output,5] = pv
}

struct jw_reselect_result scalar jw_infer_reselect_opts(
                    struct jw_model scalar M, struct jw_fit scalar F,
                    real scalar level, real scalar gamma,
                    real scalar mcerr, real scalar reps,
                    real scalar optstarts, real scalar optmaxeval,
                    real scalar opttol, real scalar optreps,
                    real rowvector outputs, real scalar cachemb,
                    real scalar progress)
{
    struct jw_reselect_result scalar R
    struct jw_rs_output scalar P
    struct jw_rs_search_result scalar SR
    real matrix S, E, EV, located
    real rowvector numerical, searchvalues
    real scalar o, j, d, initialized
    numerical = jw_rs_validate(M,F,level,gamma,mcerr,reps,optstarts,
                  optmaxeval,opttol,optreps,outputs,cachemb)
    R = jw_rs_result_init(M,F,level)
    initialized = 0
    S = J(0,0,.)
    E = J(0,0,.)
    for (o=1; o<=rows(F.A); o++) {
        if (cols(outputs)>0) {
            if (!any(outputs:==o)) continue
        }
        P = jw_rs_output_plan(M,F,o,gamma)
        d = cols(P.active)
        if (d==0) {
            jw_rs_zero_output(R,P,jw_fixed_variance(M.V,F.A[o,.]),level)
            continue
        }
        if (numerical[3]<1 | numerical[4]>reps) {
            jw_core_error("reps() too small for finite Monte Carlo quantile bounds; increase reps()")
        }
        if (!initialized) {
            S = jw_infer_sqrt(M.V)
            E = rnormal(optreps,rows(M.z),0,1)*S'
            initialized = 1
        }
        P.C = jw_rs_prepare_opts(M,F,o,E,P.hcols,cachemb)
        P.starts = jw_rs_starts(P.h,P.B,optstarts)
        /* No search routine consumes RNG. Generate validation here so
           frozen preparation and serial execution have identical draws. */
        EV = rnormal(reps,rows(M.z),0,1)*S'
        P.CV = jw_rs_prepare_opts(M,F,o,EV,P.hcols,0)
        located = J(d,4,.)
        searchvalues = J(1,4,.)
        if (progress) {
            printf("Focused inference: output %g, active dimension %g\n",o,d)
            displayflush()
        }
        for (j=1; j<=4; j++) {
            SR = jw_rs_search_opts(P.C,P.h,P.B,P.starts,numerical[1],
                        P.estimate,j,optmaxeval,opttol,progress)
            located[,j] = P.h+P.B*SR.x
            searchvalues[j] = SR.value
            jw_rs_record_search(R,P,SR,j,optreps,reps)
        }
        jw_rs_finish_output(R,P,located,searchvalues,numerical,gamma,mcerr,progress)
    }
    return(R)
}

/* Existing source-level entry point retains all outputs and 64 MB cache. */
struct jw_reselect_result scalar jw_infer_reselect(
                    struct jw_model scalar M, struct jw_fit scalar F,
                    real scalar level, real scalar gamma,
                    real scalar mcerr, real scalar reps,
                    real scalar optstarts, real scalar optmaxeval,
                    real scalar opttol, real scalar optreps)
{
    return(jw_infer_reselect_opts(M,F,level,gamma,mcerr,reps,optstarts,
                  optmaxeval,opttol,optreps,J(1,0,.),64,0))
}
