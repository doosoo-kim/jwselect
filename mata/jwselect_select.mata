/* Candidate selection and linear maps for every displayed estimate.
   Stable family row order supplies the tie break: the first exact minimum.
   No outcome-specific selection code is duplicated inside inference. */

real matrix jw_unique_masks(real matrix C, real scalar k,
                            real scalar maxjoint, real scalar addempty)
{
    real matrix U
    real colvector ord, keep, ii
    real scalar i, first, hasempty, n

    n = rows(C)
    if (n==0) {
        if (addempty) return(J(1,k,0))
        return(J(0,k,0))
    }
    /* Deduplicate in sorted order, then restore the ORIGINAL first-row
       order. Candidate order is part of the exact-tie convention. */
    keep = J(n,1,0)
    if (k==0) keep[1] = 1
    else {
        ord = order(C,1..k)
        first = ord[1]
        for (i=2; i<=n; i++) {
            if (sum(C[ord[i],.]:!=C[ord[i-1],.])==0) {
                first = min((first,ord[i]))
            }
            else {
                keep[first] = 1
                first = ord[i]
            }
        }
        keep[first] = 1
    }
    ii = selectindex(keep)
    if (rows(ii)>maxjoint) jw_core_error("joint candidate family exceeds maxjoint")
    U = C[ii,.]
    hasempty = 0
    for (i=1; i<=rows(U); i++) {
        if (sum(U[i,.])==0) {
            hasempty = 1
            break
        }
    }
    if (addempty & !hasempty) {
        if (rows(U)>=maxjoint) jw_core_error("joint candidate family including baseline exceeds maxjoint")
        U = J(1,k,0) \ U
    }
    return(U)
}

real matrix jw_joint_masks(struct jw_model scalar M,
                          struct jw_family scalar F,
                          real matrix supplied, real scalar maxjoint)
{
    real matrix C, next, localset
    real colvector ii
    real scalar g, i, j, nr, nnew, pos

    if (maxjoint<1 | maxjoint>=. | maxjoint!=floor(maxjoint)) {
        jw_core_error("maxjoint must be a positive finite integer")
    }
    if (rows(supplied)>0) {
        if (cols(supplied)!=M.k) jw_core_error("joint candidate matrix has the wrong number of columns")
        if (hasmissing(supplied) | sum((supplied:!=0):&(supplied:!=1))>0) {
            jw_core_error("joint candidate matrix must be binary")
        }
        return(jw_unique_masks(supplied,M.k,maxjoint,1))
    }
    C = J(1,M.k,0)
    for (g=1; g<=rows(M.groups); g++) {
        ii = selectindex(F.owners:==M.groups[g])
        if (rows(ii)==0) jw_core_error("a cohort has no candidates")
        localset = jw_unique_masks(F.masks[ii,.],M.k,maxjoint,0)
        nr = rows(localset)
        if (rows(C)>floor(maxjoint/nr)) {
            jw_core_error("Cartesian joint candidate family exceeds maxjoint; supply a smaller prespecified family")
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
    return(jw_unique_masks(C,M.k,maxjoint,1))
}

struct jw_fit scalar jw_run_scope(struct jw_model scalar M,
                            struct jw_family scalar F,
                            real matrix jointmatrix,
                            string scalar selectstr,
                            string scalar focusstr,
                            real scalar maxjoint,
                            string scalar wbicscope)
{
    struct jw_fit scalar out
    real scalar G, m, g, i, j, r, oi, bestrow, bestscore, score,
                isfocused, docohort, doaggregate, jointwbic, searchjoint,
                primaryjoint, aggindex,
                nout, nfamily, pos, ng, offset
    real colvector ii, localchosen, rr
    real rowvector stat, target, mask, chosenunion, A, U, unitrow,
                   aggregateA, aggregatetarget
    real matrix C, localmasks, localA, jointA, unitrows, blockA

    G = rows(M.groups)
    m = rows(M.z)
    if (G<1 | M.p<1 | M.p+M.k!=m) jw_core_error("invalid model dimensions")
    if (rows(M.V)!=m | cols(M.V)!=m | rows(M.T)!=G | cols(M.T)!=m) {
        jw_core_error("incompatible model covariance or target dimensions")
    }
    if (rows(M.postmap)!=M.p | cols(M.postmap)!=2 | rows(M.hmap)!=M.k | cols(M.hmap)!=2 | cols(M.w)!=M.p) {
        jw_core_error("incompatible cell maps or aggregation weights")
    }
    if (hasmissing(M.z) | hasmissing(M.V) | hasmissing(M.T) | hasmissing(M.w)) {
        jw_core_error("model estimates, covariance, targets, and weights must be finite")
    }
    for (r=1; r<=M.p; r++) {
        if (sum(M.groups:==M.postmap[r,1])!=1) jw_core_error("post-treatment cell has no unique model cohort")
    }
    if (cols(F.masks)!=M.k | rows(F.masks)!=rows(F.owners) | rows(F.ids)!=rows(F.owners)) {
        jw_core_error("incompatible candidate-family dimensions")
    }
    if (hasmissing(F.masks) | sum((F.masks:!=0):&(F.masks:!=1))>0) {
        jw_core_error("candidate masks must be binary")
    }
    if (selectstr!="manual" & selectstr!="wbic" & selectstr!="focused") {
        jw_core_error("unknown selection method")
    }
    if (focusstr!="cohort" & focusstr!="aggregate" & focusstr!="both") {
        jw_core_error("unknown focused-selection scope")
    }
    if (wbicscope!="cohort" & wbicscope!="joint") {
        jw_core_error("unknown Wald-BIC selection scope")
    }
    if (wbicscope=="joint" & selectstr!="wbic") {
        jw_core_error("joint Wald-BIC scope requires Wald-BIC selection")
    }
    if (selectstr=="wbic" & (M.nunit<=1 | M.nunit>=.)) {
        jw_core_error("Wald-BIC requires a finite independent-unit count greater than one")
    }
    for (i=1; i<=rows(F.owners); i++) {
        if (sum(M.groups:==F.owners[i])!=1) jw_core_error("candidate owner is not a unique model cohort")
        if (M.k>0) {
            if (sum(F.masks[i,.]:*(M.hmap[.,1]':!=F.owners[i]))>0) {
                jw_core_error("a cohort candidate contains another cohort's restrictions")
            }
        }
    }

    localchosen = J(G,1,.)
    localmasks = J(G,M.k,0)
    chosenunion = J(1,M.k,0)
    aggregatetarget = (M.w,J(1,M.k,0))
    isfocused = selectstr=="focused"
    docohort = isfocused & (focusstr=="cohort" | focusstr=="both")
    doaggregate = isfocused & (focusstr=="aggregate" | focusstr=="both")
    jointwbic = selectstr=="wbic" & wbicscope=="joint"
    searchjoint = doaggregate | jointwbic
    primaryjoint = !isfocused | (focusstr=="aggregate")
    C = J(0,M.k,.)
    if (searchjoint) C = jw_joint_masks(M,F,jointmatrix,maxjoint)

    /* Every row count is known before the loops. Preserve the original
       output/family order without repeatedly copying growing matrices. */
    nout = G+1+(isfocused & focusstr=="both")
    if (docohort) nfamily = 2*rows(F.masks)+rows(C)
    else if (doaggregate) nfamily = (G+1)*rows(C)
    else nfamily = G+1
    out.A = J(nout,m,.)
    out.targets = J(nout,m,.)
    out.selected = J(nout,M.k,.)
    out.scores = J(rows(F.masks)+rows(C),11,.)
    out.cellA = J(M.p,m,0)
    out.cell_selected = J(M.p,M.k,0)
    out.allA = J(nfamily,m,.)
    out.allT = J(nfamily,m,.)
    out.allmeta = J(nfamily,2,.)
    out.selA = J(nfamily,m,.)
    out.selT = J(nfamily,m,.)
    out.outmap = J(nout,2,.)
    localA = J(rows(F.masks),m,.)
    jointA = J(rows(C),m,.)
    pos = 0

    /* Local scores are useful diagnostics in all modes. In manual/WBIC,
       their candidate estimates use own-cohort restrictions; the final
       reported estimator subsequently uses the selected joint system. */
    for (g=1; g<=G; g++) {
        ii = selectindex(F.owners:==M.groups[g])
        if (rows(ii)==0) jw_core_error("a cohort has no candidates")
        if (selectstr=="manual" & rows(ii)!=1) jw_core_error("manual mode requires exactly one candidate per cohort")
        bestrow = ii[1]
        bestscore = .
        for (j=1; j<=rows(ii); j++) {
            i = ii[j]
            stat = jw_candidate_stats_loading(M,M.T[g,.],F.masks[i,.],A)
            localA[i,.] = A
            pos++
            out.scores[pos,.] = (1,M.groups[g],F.ids[i],stat)
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

    /* Joint WBIC evaluates the full restriction covariance for each joint
       configuration, not a sum of cohort scores. Aggregate focus uses the
       same candidate statistics but minimizes FMSC. The aggregate target
       supplies reporting diagnostics only under joint WBIC; its Wald/BIC
       criterion does not depend on aggregation weights. */
    if (searchjoint) {
        bestrow = 1
        bestscore = .
        for (i=1; i<=rows(C); i++) {
            stat = jw_candidate_stats_loading(M,aggregatetarget,C[i,.],A)
            jointA[i,.] = A
            pos++
            out.scores[pos,.] = (2,0,i,stat)
            if (jointwbic) score = stat[3]
            else score = stat[8]
            if (score<bestscore) {
                bestrow = i
                bestscore = score
            }
        }
        mask = C[bestrow,.]
        oi = bestrow
        if (jointwbic | focusstr=="aggregate") chosenunion = mask
    }

    /* A common selected joint system is factored once for every reporting
       target and underlying post-treatment cell, using multiple RHS. */
    if (primaryjoint) {
        unitrows = (I(M.p),J(M.p,M.k,0))
        blockA = jw_loadings(M.V,(M.T \ unitrows),chosenunion,M.k)
        out.A[1..G,.] = blockA[1..G,.]
        out.cellA = blockA[(G+1)..(G+M.p),.]
        for (r=1; r<=M.p; r++) out.cell_selected[r,.] = chosenunion
    }
    pos = 0
    for (g=1; g<=G; g++) {
        if (primaryjoint) U = chosenunion
        else {
            U = localmasks[g,.]
            out.A[g,.] = localA[localchosen[g],.]
            rr = selectindex(M.postmap[.,1]:==M.groups[g])
            ng = rows(rr)
            if (ng>0) {
                unitrows = J(ng,m,0)
                for (r=1; r<=ng; r++) unitrows[r,rr[r]] = 1
                out.cellA[rr,.] = jw_loadings(M.V,unitrows,U,M.k)
                for (r=1; r<=ng; r++) out.cell_selected[rr[r],.] = U
            }
        }
        A = out.A[g,.]
        out.targets[g,.] = M.T[g,.]
        out.selected[g,.] = U
        if (isfocused & focusstr=="aggregate") out.outmap[g,.] = (5,M.groups[g])
        else out.outmap[g,.] = (1,M.groups[g])

        /* Retain distinct reporting and selection maps, even when they
           coincide. Reuse each computed loading rather than solve twice. */
        if (docohort) {
            ii = selectindex(F.owners:==M.groups[g])
            for (j=1; j<=rows(ii); j++) {
                i = ii[j]
                pos++
                out.allA[pos,.] = localA[i,.]
                out.allT[pos,.] = M.T[g,.]
                out.allmeta[pos,.] = (g,1)
                out.selA[pos,.] = localA[i,.]
                out.selT[pos,.] = M.T[g,.]
            }
        }
        else if (doaggregate) {
            for (i=1; i<=rows(C); i++) {
                pos++
                out.allA[pos,.] = jw_loading(M.V,M.T[g,.],C[i,.],M.k)
                out.allT[pos,.] = M.T[g,.]
                out.allmeta[pos,.] = (g,1)
                out.selA[pos,.] = jointA[i,.]
                out.selT[pos,.] = aggregatetarget
            }
        }
        else {
            pos++
            out.allA[pos,.] = A
            out.allT[pos,.] = M.T[g,.]
            out.allmeta[pos,.] = (g,1)
            out.selA[pos,.] = A
            out.selT[pos,.] = M.T[g,.]
        }
    }

    /* This is exactly the weighted sum of primary cell estimates. With
       arbitrary cell weights it need not be a weighted sum of uniformly
       averaged cohort summaries; the displayed target metadata says so. */
    aggregateA = M.w*out.cellA
    aggindex = G+1
    out.A[aggindex,.] = aggregateA
    out.targets[aggindex,.] = aggregatetarget
    out.selected[aggindex,.] = chosenunion
    if (docohort) out.outmap[aggindex,.] = (2,0)
    else if (doaggregate) out.outmap[aggindex,.] = (3,0)
    else out.outmap[aggindex,.] = (4,0)

    if (docohort) {
        /* Selection uses the cohort summary target even when custom cell
           weights make the reporting component a different target. */
        for (g=1; g<=G; g++) {
            target = J(1,m,0)
            for (r=1; r<=M.p; r++) {
                if (M.postmap[r,1]==M.groups[g]) target[r] = M.w[r]
            }
            ii = selectindex(F.owners:==M.groups[g])
            for (j=1; j<=rows(ii); j++) {
                i = ii[j]
                pos++
                out.allA[pos,.] = jw_loading(M.V,target,F.masks[i,.],M.k)
                out.allT[pos,.] = target
                out.allmeta[pos,.] = (aggindex,g)
                out.selA[pos,.] = localA[i,.]
                out.selT[pos,.] = M.T[g,.]
            }
        }
    }
    else if (doaggregate) {
        for (i=1; i<=rows(C); i++) {
            pos++
            out.allA[pos,.] = jointA[i,.]
            out.allT[pos,.] = aggregatetarget
            out.allmeta[pos,.] = (aggindex,1)
            out.selA[pos,.] = jointA[i,.]
            out.selT[pos,.] = aggregatetarget
        }
    }
    else {
        pos++
        out.allA[pos,.] = aggregateA
        out.allT[pos,.] = aggregatetarget
        out.allmeta[pos,.] = (aggindex,1)
        out.selA[pos,.] = aggregateA
        out.selT[pos,.] = aggregatetarget
    }

    /* Both mode preserves the primary cohort-wise fit, then appends the
       independently selected aggregate-focused estimate and its family. */
    if (isfocused & focusstr=="both") {
        out.A[G+2,.] = jointA[oi,.]
        out.targets[G+2,.] = aggregatetarget
        out.selected[G+2,.] = mask
        out.outmap[G+2,.] = (3,0)
        for (i=1; i<=rows(C); i++) {
            pos++
            out.allA[pos,.] = jointA[i,.]
            out.allT[pos,.] = aggregatetarget
            out.allmeta[pos,.] = (G+2,1)
            out.selA[pos,.] = jointA[i,.]
            out.selT[pos,.] = aggregatetarget
        }
    }
    if (pos!=nfamily) jw_core_error("internal candidate family row count mismatch")
    return(out)
}

/* Preserve the six-argument numerical API and its cohort-wise WBIC default. */
struct jw_fit scalar jw_run(struct jw_model scalar M,
                            struct jw_family scalar F,
                            real matrix jointmatrix,
                            string scalar selectstr,
                            string scalar focusstr,
                            real scalar maxjoint)
{
    return(jw_run_scope(M,F,jointmatrix,selectstr,focusstr,maxjoint,"cohort"))
}
