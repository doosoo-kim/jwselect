*! jwselect 0.8.8 03nov2026
version 16.1
mata:

/* ===== jwselect_types.mata ===== */
/* jwselect 0.6.0: shared contracts. All V matrices use sampling covariance. */
struct jw_model {
    real colvector z
    real matrix V
    real matrix postmap
    real matrix hmap
    real colvector groups
    real colvector bases
    real matrix T
    real rowvector w
    real scalar p
    real scalar k
    real scalar nunit
    real scalar gap
}

struct jw_family {
    real matrix masks
    real colvector owners
    real colvector ids
}

struct jw_fit {
    real matrix A
    real matrix targets
    real matrix selected
    real matrix scores
    real matrix cellA
    real matrix cell_selected
    real matrix allA
    real matrix allT
    real matrix allmeta
    /* Actual selection criterion can have a different target from reporting. */
    real matrix selA
    real matrix selT
    real matrix outmap
}

/* ===== jwselect_design.mata ===== */
/* Build one canonical vector of baseline contrasts from unaggregated g,t margins. */
void jw_design_error(string scalar msg)
{
    errprintf("jwselect: %s\n", msg)
    _error(3498)
}

real scalar jw_findcell(real matrix cells, real scalar g, real scalar t,
                        real scalar gap)
{
    real colvector ix
    ix = vec(selectindex((cells[,1] :== g) :&
                     (abs(cells[,2] :- t) :<= 1e-8*gap)))
    if (rows(ix)==0) return(0)
    if (rows(ix)>1) jw_design_error("duplicate cohort-time cells")
    return(ix[1])
}

struct jw_model scalar jw_design(real matrix braw, real matrix Vraw,
    real matrix cells, real scalar nunit, real scalar gap, real scalar delta,
    string scalar baseline, real rowvector postwindow, string scalar aggtype,
    real matrix aggweights)
{
    struct jw_model scalar M
    real matrix C, V, P, H, L
    real colvector b, ord, gi, pi, evals
    real rowvector rr
    real scalar m, G, i, j, u, g, q, base, bi, e, npost, scale, tol, wi

    m = rows(cells)
    if (m==0 | cols(cells)!=2 | hasmissing(cells))
        jw_design_error("cells() must be a nonmissing matrix with columns cohort time")
    if (!(rows(braw)==1 | cols(braw)==1) | length(braw)!=m | hasmissing(braw))
        jw_design_error("bmat() must be a nonmissing vector matching cells()")
    if (rows(Vraw)!=m | cols(Vraw)!=m | hasmissing(Vraw))
        jw_design_error("vmat() must be a nonmissing square covariance matching bmat()")
    if (missing(nunit) | nunit<=1)
        jw_design_error("nunits() must exceed one")
    if (missing(gap) | gap<=0) jw_design_error("gap() must be positive")
    if (missing(delta) | delta<0 | delta!=floor(delta))
        jw_design_error("anticipation() must be a nonnegative integer")
    if (baseline!="cohort" & baseline!="common")
        jw_design_error("baseline() must be cohort or common")
    if (aggtype!="equalcohort" & aggtype!="equalcell")
        jw_design_error("aggregate() must be equalcohort or equalcell")
    if (length(postwindow)!=0) {
        if (length(postwindow)!=2 | hasmissing(postwindow))
            jw_design_error("postwindow() requires two nonnegative integer endpoints")
        if (any(postwindow:<0) | any(postwindow:!=floor(postwindow)) |
            postwindow[1]>postwindow[2])
            jw_design_error("postwindow() requires increasing nonnegative integer endpoints")
    }
    if (any(cells[,1]:<=0))
        jw_design_error("cells() must contain treated cohorts only, coded as positive numbers")
    ord = order(cells, (1,2))
    C = cells[ord,.]
    b = vec(braw)
    b = b[ord]
    V = Vraw[ord,ord]
    tol = 1e-8*gap
    for (i=2; i<=m; i++) {
        if (C[i,1]==C[i-1,1] & abs(C[i,2]-C[i-1,2])<=tol)
            jw_design_error("cells() contains duplicate cohort-time cells")
    }
    rr = ((C[,2]-C[,1])/gap)'
    if (any(abs(rr-round(rr)):>1e-8))
        jw_design_error("cohort and calendar times are not aligned with gap()")
    scale = max((1,max(abs(V))))
    if (max(abs(V-V'))>1e-10*scale)
        jw_design_error("vmat() is not symmetric")
    V = (V+V')/2
    evals = vec(symeigenvalues(V))
    if (min(evals)<-1e-10*scale)
        jw_design_error("vmat() is not positive semidefinite")

    M.groups = uniqrows(C[,1])
    G = rows(M.groups)
    q = M.groups[1]
    M.bases = J(G,1,.)
    P = J(0,m,0)
    H = J(0,m,0)
    M.postmap = J(0,2,.)
    M.hmap = J(0,2,.)
    for (i=1; i<=G; i++) {
        g = M.groups[i]
        base = g-gap*(delta+1)
        if (baseline=="common") base = q-gap*(delta+1)
        M.bases[i] = base
        bi = jw_findcell(C,g,base,gap)
        if (bi==0)
            jw_design_error(sprintf("baseline time %g is unavailable for cohort %g",base,g))
        gi = selectindex(C[,1]:==g)
        npost = 0
        if (length(postwindow)) {
            for (e=postwindow[1]; e<=postwindow[2]; e++) {
                j = jw_findcell(C,g,g+gap*e,gap)
                if (j==0)
                    jw_design_error(sprintf("event time %g is unavailable for cohort %g",e,g))
                rr = J(1,m,0)
                rr[j] = 1
                rr[bi] = -1
                P = P\rr
                M.postmap = M.postmap\C[j,.]
                npost++
            }
        }
        else {
            for (u=1; u<=rows(gi); u++) {
                j = gi[u]
                if (C[j,2]>=g-tol) {
                    rr = J(1,m,0)
                    rr[j] = 1
                    rr[bi] = -1
                    P = P\rr
                    M.postmap = M.postmap\C[j,.]
                    npost++
                }
            }
        }
        if (npost==0)
            jw_design_error(sprintf("cohort %g has no supported post-treatment cells",g))
        for (u=1; u<=rows(gi); u++) {
            j = gi[u]
            if (C[j,2]<g-gap*delta-tol & j!=bi) {
                rr = J(1,m,0)
                rr[j] = 1
                rr[bi] = -1
                H = H\rr
                M.hmap = M.hmap\C[j,.]
            }
        }
    }
    M.p = rows(P)
    M.k = rows(H)
    L = P\H
    M.z = L*b
    M.V = L*V*L'
    M.V = (M.V+M.V')/2
    M.T = J(G,M.p+M.k,0)
    M.w = J(1,M.p,0)
    for (i=1; i<=G; i++) {
        pi = selectindex(M.postmap[,1]:==M.groups[i])
        M.T[i,pi] = J(1,rows(pi),1/rows(pi))
        M.w[1,pi] = J(1,rows(pi),1/(G*rows(pi)))
    }
    if (aggtype=="equalcell") M.w = J(1,M.p,1/M.p)
    if (rows(aggweights)>0 | cols(aggweights)>0) {
        if (rows(aggweights)!=M.p | cols(aggweights)!=3 | hasmissing(aggweights))
            jw_design_error("aggweights() must list every target cell: cohort time weight")
        if (any(aggweights[,3]:<0) | abs(sum(aggweights[,3])-1)>1e-8)
            jw_design_error("aggweights() must be nonnegative and sum to one")
        M.w = J(1,M.p,0)
        rr = J(1,M.p,0)
        for (i=1; i<=rows(aggweights); i++) {
            wi = jw_findcell(M.postmap,aggweights[i,1],aggweights[i,2],gap)
            if (wi==0) jw_design_error("aggweights() contains a cell outside the target")
            if (rr[wi]) jw_design_error("aggweights() contains duplicate cells")
            rr[wi] = 1
            M.w[wi] = aggweights[i,3]
        }
        if (any(rr:==0)) jw_design_error("aggweights() must explicitly list every target cell")
        /* Preserve user weights exactly; no implicit normalization. */
    }
    M.nunit = nunit
    M.gap = gap
    return(M)
}

/* ===== jwselect_candidates.mata ===== */
/* Local candidate families. Their masks always index the canonical H vector. */
real scalar jw_mask_exists(real matrix masks, real rowvector mask)
{
    real scalar j
    if (cols(masks)!=cols(mask)) return(0)
    for (j=1; j<=rows(masks); j++) {
        if (sum(abs(masks[j,.]-mask))==0) return(1)
    }
    return(0)
}

struct jw_family scalar jw_candidates(struct jw_model scalar M,
    string scalar selectmethod, string scalar candidates, real matrix prematrix,
    real matrix candmatrix, real scalar maxcand)
{
    struct jw_family scalar F
    real matrix localm
    real colvector hi, ci, ip, ord, counts
    real rowvector mask
    real scalar G, i, j, s, g, h, kg, nc, bit, u, haspre, hascand, pos, used

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
    /* Count first: all-subsets families can have thousands of rows.
       Repeated concatenation copied the entire growing family each time. */
    G = rows(M.groups)
    counts = J(G,1,1)
    for (i=1; i<=G; i++) {
        g = M.groups[i]
        kg = sum(M.hmap[,1]:==g)
        if (selectmethod=="manual") continue
        if (hascand) {
            ci = vec(selectindex(candmatrix[,1]:==g))
            if (rows(ci)==0)
                jw_design_error(sprintf("candmatrix() has no rows for cohort %g",g))
            for (j=1; j<=rows(ci); j++) {
                if (M.k) {
                    if (sum(candmatrix[ci[j],2..(M.k+1)])>0) counts[i] = counts[i] + 1
                }
            }
        }
        else if (candidates=="nested") counts[i] = kg+1
        else {
            if (kg>floor(ln(maxcand)/ln(2)+1e-10))
                jw_design_error("all-subsets candidate count exceeds maxcand()")
            counts[i] = 2^kg
        }
        if (counts[i]>maxcand)
            jw_design_error("local candidate count exceeds maxcand()")
    }
    F.masks = J(sum(counts),M.k,0)
    F.owners = J(sum(counts),1,.)
    F.ids = J(sum(counts),1,.)
    pos = 0
    for (i=1; i<=G; i++) {
        g = M.groups[i]
        /* Scalar predicates can yield a 1x0 index; enforce column shape. */
        hi = vec(selectindex(M.hmap[,1]:==g))
        kg = rows(hi)
        mask = J(1,M.k,0)
        localm = J(counts[i],M.k,0)
        used = 1
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
                if (jw_mask_exists(localm[1..used,.],mask))
                    jw_design_error("candmatrix() contains duplicate nonbaseline candidates")
                used++
                localm[used,.] = mask
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
                    localm[j+1,.] = mask
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
                    localm[s+1,.] = mask
                }
            }
        }
        F.masks[(pos+1)..(pos+rows(localm)),.] = localm
        F.owners[(pos+1)..(pos+rows(localm)),1] = J(rows(localm),1,g)
        F.ids[(pos+1)..(pos+rows(localm)),1] = (1::rows(localm))
        pos = pos+rows(localm)
    }
    return(F)
}

/* ===== jwselect_blocks.mata ===== */
/* Conditional-index inputs. Each cohort/date is one selectable block;
   canonical H retains every identified coefficient in that block. */

real matrix jw_block_map(real matrix hmap)
{
    real matrix dates, blocks
    real scalar j

    if (rows(hmap)==0) return(J(0,3,.))
    if (cols(hmap)!=2 | hasmissing(hmap))
        jw_design_error("index hmap must contain finite cohort and time columns")
    dates = uniqrows(hmap[order(hmap,(1,2)),.])
    blocks = (dates,J(rows(dates),1,.))
    for (j=1; j<=rows(dates); j++) {
        blocks[j,3] = sum((hmap[,1]:==dates[j,1]):&(hmap[,2]:==dates[j,2]))
    }
    return(blocks)
}

void jw_block_masks_check(real matrix hmap, real matrix masks)
{
    real matrix blocks
    real colvector ix
    real scalar i, j

    if (cols(masks)!=rows(hmap))
        jw_design_error("index candidate mask must have one bit per coefficient in hmap")
    if (hasmissing(masks) | any(vec((masks:!=0):&(masks:!=1))))
        jw_design_error("index candidate masks must be binary")
    blocks = jw_block_map(hmap)
    for (j=1; j<=rows(blocks); j++) {
        ix = vec(selectindex((hmap[,1]:==blocks[j,1]):&(hmap[,2]:==blocks[j,2])))
        for (i=1; i<=rows(masks); i++) {
            if (min(masks[i,ix])!=max(masks[i,ix]))
                jw_design_error("a conditional-index date block must be selected in full; split-block masks are not allowed")
        }
    }
}

struct jw_model scalar jw_design_index(real matrix braw, real matrix Vraw,
    real matrix cells, real matrix hraw, real matrix Vhraw,
    real matrix Vdhraw, real matrix hmap, real scalar nunit,
    real scalar gap, real scalar delta, string scalar baseline,
    real rowvector postwindow, string scalar aggtype, real matrix aggweights)
{
    struct jw_model scalar M
    real matrix oldT, cross, blocks, expected, Vh, augmented
    real colvector h, evals, gi, seen
    real scalar k, i, j, u, bi, scale

    if (baseline!="cohort")
        jw_design_error("conditional-index blocks retain the source cohort baseline")
    M = jw_design(braw,Vraw,cells,nunit,gap,delta,baseline,postwindow,aggtype,aggweights)
    expected = M.hmap
    oldT = M.T[,1..M.p]
    k = rows(hmap)
    if (k==0) {
        if (rows(expected)>0)
            jw_design_error("index inputs omit eligible preperiod coefficient blocks")
        if (length(hraw)>0 | length(Vhraw)>0 | length(Vdhraw)>0)
            jw_design_error("empty index hmap has nonempty coefficient or covariance inputs")
        M.hmap = J(0,2,.)
        M.k = 0
        M.z = M.z[1..M.p]
        M.V = M.V[1..M.p,1..M.p]
        M.T = oldT
        return(M)
    }
    if (cols(hmap)!=2 | hasmissing(hmap))
        jw_design_error("index hmap must contain finite cohort and time columns")
    if (!(rows(hraw)==1 | cols(hraw)==1) | length(hraw)!=k | hasmissing(hraw))
        jw_design_error("index coefficient vector must match hmap")
    if (rows(Vhraw)!=k | cols(Vhraw)!=k | hasmissing(Vhraw))
        jw_design_error("index covariance must be a finite square matrix matching hmap")
    if (rows(Vdhraw)!=rows(cells) | cols(Vdhraw)!=k | hasmissing(Vdhraw))
        jw_design_error("ATT-index cross-covariance must have one row per source cell and one column per index coefficient")
    /* Require every eligible source date, including the whole lead block.
       This prevents partial support from silently changing candidate sets. */
    blocks = jw_block_map(hmap)
    if (rows(blocks)!=rows(expected))
        jw_design_error("index date blocks must match all eligible source preperiods")
    seen = J(rows(expected),1,0)
    for (i=1; i<=rows(blocks); i++) {
        j = jw_findcell(expected,blocks[i,1],blocks[i,2],gap)
        if (j==0)
            jw_design_error("index hmap contains an ineligible date, baseline, or unknown cohort")
        if (seen[j])
            jw_design_error("index hmap contains ambiguous duplicate dates within the time tolerance")
        seen[j] = 1
    }
    /* Normalize only the date labels to the exact source grid. Coefficient
       order stays unchanged, preserving the adapter's coefficient map. */
    M.hmap = hmap
    for (i=1; i<=k; i++) {
        j = jw_findcell(expected,hmap[i,1],hmap[i,2],gap)
        if (j==0) jw_design_error("index hmap has an unsupported source date")
        M.hmap[i,.] = expected[j,.]
    }
    h = vec(hraw)
    scale = max((1,max(abs(Vhraw))))
    if (max(abs(Vhraw-Vhraw'))>1e-10*scale)
        jw_design_error("index covariance is not symmetric")
    Vh = (Vhraw+Vhraw')/2
    evals = vec(symeigenvalues(Vh))
    if (min(evals)<=0 | min(evals)<=1e-10*max(evals))
        jw_design_error("index coefficient covariance is singular or ill-conditioned; complete identified blocks are required")
    /* Source outcome placebo margins may be deterministic delta-method
       combinations of H; PSD, rather than PD, is appropriate before the
       post-cell restriction below. */
    augmented = (Vraw,Vdhraw \ Vdhraw',Vh)
    scale = max((1,max(abs(augmented))))
    evals = vec(symeigenvalues((augmented+augmented')/2))
    if (min(evals)<-1e-10*scale)
        jw_design_error("joint ATT-index covariance is not positive semidefinite")
    cross = J(M.p,k,.)
    for (i=1; i<=M.p; i++) {
        j = jw_findcell(cells,M.postmap[i,1],M.postmap[i,2],gap)
        gi = vec(selectindex(M.groups:==M.postmap[i,1]))
        u = gi[1]
        bi = jw_findcell(cells,M.groups[u],M.bases[u],gap)
        cross[i,.] = Vdhraw[j,.]-Vdhraw[bi,.]
    }
    M.z = M.z[1..M.p] \ h
    M.V = (M.V[1..M.p,1..M.p],cross \ cross',Vh)
    M.V = (M.V+M.V')/2
    evals = vec(symeigenvalues(M.V))
    if (min(evals)<=0 | min(evals)<=1e-10*max(evals))
        jw_design_error("retained post-ATT/index covariance is singular or ill-conditioned")
    M.k = k
    M.T = (oldT,J(rows(M.groups),k,0))
    return(M)
}

struct jw_family scalar jw_index_candidates(struct jw_model scalar M,
    string scalar selectmethod, string scalar candidates, real matrix prematrix,
    real matrix candmatrix, real scalar maxcand)
{
    struct jw_model scalar P
    struct jw_family scalar C, F
    real matrix blocks, compressed
    real colvector ix
    real scalar j, i, hascand

    blocks = jw_block_map(M.hmap)
    P = M
    P.hmap = blocks[,1..2]
    P.k = rows(blocks)
    hascand = (rows(candmatrix)>0 | cols(candmatrix)>0)
    compressed = J(0,0,.)
    if (hascand) {
        if (cols(candmatrix)!=M.k+1 | hasmissing(candmatrix))
            jw_design_error("index candmatrix() requires cohort followed by coefficient bits in hmap order")
        compressed = J(rows(candmatrix),P.k+1,.)
        compressed[,1] = candmatrix[,1]
        if (M.k>0) {
            jw_block_masks_check(M.hmap,candmatrix[,2..(M.k+1)])
            for (j=1; j<=P.k; j++) {
                ix = vec(selectindex((M.hmap[,1]:==blocks[j,1]):&(M.hmap[,2]:==blocks[j,2])))
                compressed[,j+1] = candmatrix[,ix[1]+1]
            }
        }
    }
    C = jw_candidates(P,selectmethod,candidates,prematrix,compressed,maxcand)
    F = C
    F.masks = J(rows(C.masks),M.k,0)
    for (j=1; j<=P.k; j++) {
        ix = vec(selectindex((M.hmap[,1]:==blocks[j,1]):&(M.hmap[,2]:==blocks[j,2])))
        for (i=1; i<=rows(ix); i++) F.masks[,ix[i]] = C.masks[,j]
    }
    return(F)
}

/* ===== jwselect_core.mata ===== */
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

/* Joint zero-precontrast test with the SAME covariance as the MD fit.
   J is reported through the Wald/MD identity, not a second fitted model.
   An empty restriction set has no chi-squared test (p is missing). */
real rowvector jw_pretrend_test(real colvector z, real matrix V,
                               real rowvector mask, real scalar k)
{
    real scalar r, W, pv

    if (missing(k) | k<0 | k!=floor(k) | k>rows(z) |
        cols(z)!=1 | rows(V)!=rows(z) | cols(V)!=rows(z) |
        cols(mask)!=k | rows(mask)!=1) {
        jw_core_error("incompatible joint pretrend-test dimensions")
    }
    if (hasmissing(z) | hasmissing(V) | hasmissing(mask)) {
        jw_core_error("joint pretrend-test inputs must be finite")
    }
    if (sum((mask:!=0):&(mask:!=1))>0) {
        jw_core_error("joint pretrend-test mask must be binary")
    }
    r = sum(mask)
    W = jw_wald(V,z,mask,k)
    pv = .
    if (r>0) pv = chi2tail(r,W)
    return((r,W,W,pv))
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

/* ===== jwselect_select.mata ===== */
/* Candidate selection and linear maps for every displayed estimate.
   Focused exact ties favor fewer restrictions, then original family order.
   WBIC retains first-in-family ties. Inference families are stored
   in focused tie-priority order, including reporting/selection maps. */

real colvector jw_focused_order(real matrix masks)
{
    real scalar n
    real colvector rank
    n = rows(masks)
    if (n==0) return(J(0,1,.))
    rank = J(n,1,0)
    if (cols(masks)>0) rank = rowsum(masks)
    return(order((rank,(1::n)),(1,2)))
}

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
    real colvector ii, localchosen, rr, jointorder
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
    jointorder = jw_focused_order(C)

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
            if (score<bestscore | (isfocused & score==bestscore &
                sum(F.masks[i,.])<sum(F.masks[bestrow,.]))) {
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
            if (score<bestscore | (isfocused & score==bestscore &
                sum(C[i,.])<sum(C[bestrow,.]))) {
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
            ii = ii[jw_focused_order(F.masks[ii,.])]
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
            for (j=1; j<=rows(C); j++) {
                i = jointorder[j]
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
            ii = ii[jw_focused_order(F.masks[ii,.])]
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
        for (j=1; j<=rows(C); j++) {
            i = jointorder[j]
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
        for (j=1; j<=rows(C); j++) {
            i = jointorder[j]
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

/* ===== jwselect_infer.mata ===== */
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

/* ===== jwselect_reselect.mata ===== */
/* Reselection-based two-step focused inference, in sampling-variance units.

   At every trial nuisance mean mu and every common Gaussian draw E, select
   the first minimizer of the ORIGINAL focused criterion in the stored
   fewer-restrictions-first tie order (then original order). Reporting uses
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

/* ===== jwselect_wbic.mata ===== */
/* Selection-adjusted WBIC inference, sampling-covariance units.
   A finite COVER of the nuisance ellipsoid, not a grid approximation.
   Interval bounds retain every possible WBIC winner inside each cell.
   The per-draw maximum therefore bounds the continuum of studentized
   selected errors. Exact binomial bounds account for simulation error.
   See docs/wbic_inference_implementation.txt for the coverage argument. */

struct jw_wbic_result {
    real matrix table
    real matrix diagnostics
}

struct jw_wbic_score {
    real rowvector idx
    real matrix W
    real rowvector boxnorm
    real scalar penalty
}

struct jw_wbic_plan {
    real matrix masks
    real matrix tuples
    real colvector owners
    struct jw_wbic_score colvector scores
    real matrix L
    real colvector h
    real scalar radius
    real scalar p
    real scalar k
    real scalar d
    real scalar nout
    real rowvector outputs
    real matrix actions
    real matrix sigmas
    real matrix shifts
    real matrix supports
    real matrix directions
}

/* Deterministic longest-side partition of [-1,1]^d. Dropping boxes is
   allowed only when their minimum squared distance exceeds one.
   Row layout is (center,width), with closed cells on shared faces. */
real matrix jw_wbic_cells(real scalar d, real scalar budget)
{
    real matrix B
    real rowvector child, closest
    real scalar i, j, q, r, best, size, dim
    B = (J(1,d,0),J(1,d,1))
    while (rows(B)<budget) {
        best = -1
        q = 1
        for (i=1; i<=rows(B); i++) {
            size = sum(B[i,(d+1)..(2*d)]:^2)
            if (size>best) {
                best=size
                q=i
            }
        }
        dim = 1
        for (j=2; j<=d; j++) {
            if (B[q,d+j]>B[q,d+dim]) dim=j
        }
        child = B[q,.]
        child[d+dim] = child[d+dim]/2
        B[q,.] = child
        B[q,dim] = child[dim]-child[d+dim]
        child[dim] = child[dim]+child[d+dim]
        B = B \ child
        /* Closed boundary cells are retained with outward roundoff slack. */
        for (i=rows(B); i>=1; i--) {
            closest = J(1,d,0)
            for (j=1; j<=d; j++) closest[j]=max((0,abs(B[i,j])-B[i,d+j]))
            if (sum(closest:^2)>1+1e-12) {
                if (i==1) B=B[2..rows(B),.]
                else if (i==rows(B)) B=B[1..(i-1),.]
                else B=B[1..(i-1),.] \ B[(i+1)..rows(B),.]
            }
        }
    }
    return(B)
}

struct jw_wbic_plan scalar jw_wbic_prepare(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar gamma,
    real rowvector outputs)
{
    struct jw_wbic_plan scalar P
    real matrix smasks, Q, eig, A, D
    real rowvector active, idx, target, empty, ownercols
    real colvector ci
    real scalar i,j,g,o,t,ns,nq,m,v,found

    if (scope=="cohort" & rows(joint)>0) {
        jw_core_error("cohort WBIC inference requires the full Cartesian cohort family, not jointmatrix()")
    }
    P.masks=jw_joint_masks(M,C,joint,maxjoint)
    ns=rows(P.masks)
    P.p=M.p
    P.k=M.k
    P.outputs=outputs
    P.nout=cols(outputs)
    m=rows(M.z)
    active=selectindex(colsum(P.masks):>0)
    P.d=cols(active)
    P.h=M.z[(M.p+1)..m,1]
    P.radius=sqrt(invchi2(P.d,1-gamma))
    Q=M.V[M.p:+active,M.p:+active]
    eig=symeigenvalues((Q+Q')/2)
    if (hasmissing(eig) | min(eig)<=1e-12*max(abs(eig))) {
        jw_core_error("WBIC inference requires nonsingular active precontrast covariance")
    }
    P.L=J(M.k,P.d,0)
    P.L[active,.]=P.radius*jw_infer_sqrt(Q)
    if (scope=="joint") {
        smasks=P.masks
        P.owners=J(ns,1,1)
        P.tuples=(1::ns)
    }
    else {
        smasks=C.masks
        P.owners=C.owners
        P.tuples=J(ns,rows(M.groups),.)
        for (t=1; t<=ns; t++) {
            for (g=1; g<=rows(M.groups); g++) {
                ci=selectindex(C.owners:==M.groups[g])
                ownercols=(M.hmap[.,1]':==M.groups[g])
                found=0
                for (j=1; j<=rows(ci); j++) {
                    if (sum(abs(C.masks[ci[j],.]-P.masks[t,.]:*ownercols))==0) {
                        P.tuples[t,g]=ci[j]
                        found=1
                        break
                    }
                }
                if (!found) jw_core_error("WBIC inference tuple is not in the cohort candidate family")
            }
        }
    }
    nq=rows(smasks)
    P.scores=jw_wbic_score(nq)
    for (i=1; i<=nq; i++) {
        idx=selectindex(smasks[i,.]:==1)
        P.scores[i].idx=idx
        P.scores[i].penalty=cols(idx)*ln(M.nunit)
        if (cols(idx)>0) {
            Q=M.V[M.p:+idx,M.p:+idx]
            P.scores[i].W=jw_infer_sqrt(cholsolve(Q,I(cols(idx))))'
            D=P.scores[i].W*P.L[idx,.]
            P.scores[i].boxnorm=sqrt(colsum(D:^2))
        }
    }
    P.actions=J(ns*P.nout,m,.)
    P.sigmas=J(ns,P.nout,.)
    P.shifts=J(ns,P.nout,.)
    P.supports=J(ns,P.nout,.)
    P.directions=J(ns*P.nout,P.d,.)
    for (t=1; t<=ns; t++) {
        A=jw_loadings(M.V,F.targets[outputs,.],P.masks[t,.],M.k)
        for (o=1; o<=P.nout; o++) {
            i=(o-1)*ns+t
            P.actions[i,.]=A[o,.]
            v=jw_fixed_variance(M.V,A[o,.])
            if (v<=1e-12*max((jw_fixed_variance(M.V,F.targets[outputs[o],.]),1e-300))) {
                jw_core_error("WBIC inference requires nondegenerate variance for every candidate target")
            }
            P.sigmas[t,o]=sqrt(v)
            P.shifts[t,o]=A[o,(M.p+1)..m]*P.h
            P.directions[i,.]=A[o,(M.p+1)..m]*P.L
            P.supports[t,o]=sqrt(sum(P.directions[i,.]:^2))
        }
    }
    return(P)
}

/* Per-draw upper envelope; exposed for independent pathwise validation.
   E contains draws in ORIGINAL full z order. The result has one column
   per requested output. No data, RNG, or Stata state is read here. */
real matrix jw_wbic_cell_envelope(struct jw_wbic_plan scalar P,
    real matrix E, real rowvector center, real rowvector width)
{
    real matrix lower, upper, possible, U, X, Eh, Y
    real colvector normc,normg,lc,uc,lg,ug,viable,bestupper,v,uglob,ubox,covered
    real rowvector idx,mu,dir
    real colvector owners,ii
    real scalar i,j,g,t,o,ns,nq,nb,r,rad,err,shift,supbox,pos
    ns=rows(P.masks)
    nq=rows(P.owners)
    nb=rows(E)
    lower=J(nb,nq,0)
    upper=J(nb,nq,0)
    mu=(P.h+P.L*center')'
    Eh=E[.,(P.p+1)..(P.p+P.k)]
    for (i=1; i<=nq; i++) {
        idx=P.scores[i].idx
        if (cols(idx)==0) continue
        X=(Eh[.,idx]+J(nb,1,1)*mu[idx])*P.scores[i].W'
        normc=sqrt(rowsum(X:^2))
        X=(Eh[.,idx]+J(nb,1,1)*P.h[idx]')*P.scores[i].W'
        normg=sqrt(rowsum(X:^2))
        /* The global ball's restricted whitened operator norm is radius.
           Inflate all norm bounds outwards for floating-point arithmetic. */
        rad=P.scores[i].boxnorm*width'
        rad=rad+1e-10*(1+rad)
        r=P.radius+1e-10*(1+P.radius)
        lc=rowmax((J(nb,1,0),normc:-rad)):^2
        uc=(normc:+rad):^2
        lg=rowmax((J(nb,1,0),normg:-r)):^2
        ug=(normg:+r):^2
        lower[.,i]=rowmax((lc,lg)):-P.scores[i].penalty
        upper[.,i]=rowmin((uc,ug)):-P.scores[i].penalty
    }
    if (hasmissing(lower) | hasmissing(upper)) jw_core_error("nonfinite WBIC score bounds")
    possible=J(nb,nq,0)
    owners=uniqrows(sort(P.owners,1))
    for (g=1; g<=rows(owners); g++) {
        ii=selectindex(P.owners:==owners[g])
        bestupper=rowmin(upper[.,ii])
        for (j=1; j<=rows(ii); j++) {
            i=ii[j]
            possible[.,i]=lower[.,i]:<=bestupper:+1e-9:*(1:+abs(bestupper):+abs(lower[.,i]))
        }
    }
    U=J(nb,P.nout,0)
    covered=J(nb,1,0)
    for (t=1; t<=ns; t++) {
        viable=J(nb,1,1)
        for (g=1; g<=cols(P.tuples); g++) viable=viable:*possible[.,P.tuples[t,g]]
        if (sum(viable)==0) continue
        covered=covered:+viable
        for (o=1; o<=P.nout; o++) {
            pos=(o-1)*ns+t
            dir=P.directions[pos,.]
            shift=P.shifts[t,o]+dir*center'
            supbox=abs(dir)*width'
            v=E*P.actions[pos,.]'
            ubox=abs(v:+shift):+supbox
            uglob=abs(v:+P.shifts[t,o]):+P.supports[t,o]
            if (hasmissing(ubox) | hasmissing(uglob)) jw_core_error("nonfinite WBIC candidate error bounds")
            v=rowmin((ubox,uglob)):/P.sigmas[t,o]
            if (hasmissing(v)) jw_core_error("nonfinite WBIC studentized error bound")
            v=(v:+1e-10:*(1:+abs(v))):*viable
            U[.,o]=rowmax((U[.,o],v))
        }
    }
    if (min(covered)<1) jw_core_error("WBIC cell bounds retained no possible selection; check covariance conditioning")
    if (hasmissing(U)) jw_core_error("nonfinite WBIC cell envelope")
    return(U)
}

struct jw_wbic_result scalar jw_infer_wbic(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps,
    real scalar ncells, real rowvector outputs, real scalar progress)
{
    struct jw_wbic_result scalar R
    struct jw_wbic_plan scalar P
    real matrix boxes,S,E,U,draws
    real rowvector active,rr,critical,pvalue,ranks
    real colvector ordered
    real scalar alpha,tail,delta,i,o,b,start,stop,nout,est,se,kq
    nout=rows(F.A)
    if (level<=0 | level>=100 | missing(level) | gamma<=0 | mcerr<=0 |
        missing(gamma) | missing(mcerr) | reps<1 | reps!=floor(reps) |
        missing(reps) | ncells<1 | ncells>1024 | ncells!=floor(ncells)) {
        jw_core_error("invalid WBIC inference probability or simulation options")
    }
    alpha=1-level/100
    if (gamma+mcerr>=alpha) jw_core_error("gamma()+mcerr() must be smaller than 1-level()/100")
    if (scope!="cohort" & scope!="joint") jw_core_error("unknown WBIC inference scope")
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs))) {
        jw_core_error("invalid WBIC inference output indices")
    }
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,8,.)
    for (o=1; o<=nout; o++) {
        R.table[o,1]=F.A[o,.]*M.z
        R.table[o,2]=sqrt(jw_fixed_variance(M.V,F.A[o,.]))
        R.diagnostics[o,.]=(o,.,.,.,0,.,.,1)
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        for (i=1; i<=cols(outputs); i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.diagnostics[o,.]=(o,0,1,0,0,invnormal(1-alpha/2),.,2)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    boxes=jw_wbic_cells(P.d,ncells)
    /* Family-wise MC failure allocation covers every reported output
       and cell. Statistical confidence intervals remain MARGINAL. */
    delta=mcerr/(rows(boxes)*cols(outputs))
    tail=1-alpha+gamma+mcerr
    ranks=jw_rs_quantile_ranks(reps,tail,delta)
    kq=ranks[2]
    if (kq>reps) jw_core_error("reps() is too small for WBIC level(), gamma(), mcerr(), and wbiccells(); increase reps()")
    critical=J(1,cols(outputs),0)
    pvalue=J(1,cols(outputs),0)
    S=jw_infer_sqrt(M.V)
    for (b=1; b<=rows(boxes); b++) {
        if (progress) printf("WBIC inference: cell %g of %g; %g joint candidates\n",b,rows(boxes),rows(P.masks))
        draws=J(reps,cols(outputs),.)
        for (start=1; start<=reps; start=start+256) {
            stop=min((reps,start+255))
            E=rnormal(stop-start+1,rows(M.z),0,1)*S'
            U=jw_wbic_cell_envelope(P,E,boxes[b,1..P.d],boxes[b,(P.d+1)..(2*P.d)])
            draws[|start,1\stop,cols(outputs)|]=U
        }
        for (i=1; i<=cols(outputs); i++) {
            ordered=sort(draws[.,i],1)
            critical[i]=max((critical[i],ordered[kq]))
            o=outputs[i]
            est=abs(R.table[o,1])/R.table[o,2]
            /* >= handles atoms conservatively. This p-value inverts the
               same order-statistic construction, apart from boundary ties. */
            pvalue[i]=max((pvalue[i],jw_rs_probability_upper(sum(draws[.,i]:>=est),reps,delta)))
        }
    }
    for (i=1; i<=cols(outputs); i++) {
        o=outputs[i]
        est=R.table[o,1]
        se=R.table[o,2]
        R.table[o,3]=est-critical[i]*se
        R.table[o,4]=est+critical[i]*se
        R.table[o,5]=min((1,gamma+mcerr+pvalue[i]))
        R.diagnostics[o,.]=(o,P.d,rows(P.masks),rows(boxes),reps,critical[i],kq,0)
    }
    return(R)
}

/* Sharper, still conservative, WBIC calculation.  The legacy engine above
   is deliberately retained without changes.  The following three changes
   are pathwise tightenings, not an unconstrained nuisance optimization:
   signed box/ball intersection supports, winner-implied Wald caps, and
   zeroing draws outside the noise event already charged to gamma. */
struct jw_wbic_sharp_plan {
    real rowvector active
    real matrix noiseW
    real matrix adjustment_scales
}

struct jw_wbic_sharp_cell {
    real rowvector center
    real rowvector width
    real matrix support
    real rowvector score_radius
}

/* Upper bound to sup d*v over the closed box intersected with ||v||<=1.
   Every positive-lambda Lagrange-dual evaluation is an upper bound, so
   finite bisection cannot produce an optimistic support value. */
real scalar jw_wbic_intersection_support(real rowvector d,
    real rowvector center, real rowvector width)
{
    real rowvector lo,hi,z
    real scalar ans,box,ball,lam,left,right,value,normsq,j,it,scale
    if (cols(d)!=cols(center) | cols(d)!=cols(width) | any(width:<0) |
        hasmissing(d) | hasmissing(center) | hasmissing(width)) {
        jw_core_error("invalid WBIC signed-support inputs")
    }
    box=d*center'+abs(d)*width'
    ball=sqrt(sum(d:^2))
    ans=min((box,ball))
    if (ball==0) return(0)
    lo=center-width
    hi=center+width
    /* If the ball maximizer is in the box, the whole-ball support is
       attained in the intersection and no dual search is needed. */
    z=d/ball
    if (all(z:>=lo) & all(z:<=hi)) return(ball+1e-12*(1+ball))
    z=J(1,cols(d),0)
    for (j=1; j<=cols(d); j++) {
        if (d[j]>0) z[j]=hi[j]
        else if (d[j]<0) z[j]=lo[j]
        else z[j]=min((hi[j],max((lo[j],0))))
    }
    if (sum(z:^2)<=1) return(box+1e-12*(1+abs(box)+ball))
    left=0
    right=max((ball,1e-100))
    /* Find a finite bracket.  Retain each dual value even if the
       intersection touches the sphere and no finite root exists. */
    for (it=1; it<=80; it++) {
        lam=right
        for (j=1; j<=cols(d); j++) z[j]=min((hi[j],max((lo[j],d[j]/(2*lam)))))
        value=lam+sum(d:*z-lam:*(z:^2))
        scale=1+abs(lam)+sum(abs(d:*z))+lam*sum(z:^2)
        value=value+1e-12*scale
        if (!missing(value)) ans=min((ans,value))
        normsq=sum(z:^2)
        if (normsq<=1) break
        left=right
        right=2*right
    }
    for (it=1; it<=48; it++) {
        lam=(left+right)/2
        if (lam<=0 | missing(lam)) break
        for (j=1; j<=cols(d); j++) z[j]=min((hi[j],max((lo[j],d[j]/(2*lam)))))
        value=lam+sum(d:*z-lam:*(z:^2))
        scale=1+abs(lam)+sum(abs(d:*z))+lam*sum(z:^2)
        value=value+1e-12*scale
        if (!missing(value)) ans=min((ans,value))
        if (sum(z:^2)>1) left=lam
        else right=lam
    }
    return(ans+1e-12*(1+abs(ans)+ball))
}

struct jw_wbic_sharp_plan scalar jw_wbic_sharp_prepare(struct jw_wbic_plan scalar P)
{
    struct jw_wbic_sharp_plan scalar S
    real matrix Q
    real rowvector idx,d
    real scalar pos,t,g,ns
    ns=rows(P.masks)
    S.active=selectindex(rowsum(P.L:^2):>0)'
    Q=P.L[S.active,.]*P.L[S.active,.]'/(P.radius^2)
    S.noiseW=jw_infer_sqrt(cholsolve(Q,I(P.d)))'
    S.adjustment_scales=J(rows(P.actions),cols(P.tuples),0)
    for (pos=1; pos<=rows(P.actions); pos++) {
        t=mod(pos-1,ns)+1
        for (g=1; g<=cols(P.tuples); g++) {
            idx=P.scores[P.tuples[t,g]].idx
            if (cols(idx)==0) continue
            d=P.actions[pos,P.p:+idx]*P.L[idx,.]
            S.adjustment_scales[pos,g]=sqrt(sum(d:^2))/P.radius
        }
    }
    if (hasmissing(S.noiseW) | hasmissing(S.adjustment_scales)) {
        jw_core_error("nonfinite WBIC sharp preparation")
    }
    return(S)
}

struct jw_wbic_sharp_cell scalar jw_wbic_sharp_cell_prepare(
    struct jw_wbic_plan scalar P, real rowvector center, real rowvector width)
{
    struct jw_wbic_sharp_cell scalar B
    real matrix D
    real rowvector idx
    real scalar pos,i,oldrad,quadrad
    B.center=center
    B.width=width
    B.score_radius=J(1,rows(P.owners),0)
    for (i=1; i<=rows(P.owners); i++) {
        idx=P.scores[i].idx
        if (cols(idx)==0) continue
        D=P.scores[i].W*P.L[idx,.]
        oldrad=P.scores[i].boxnorm*width'
        /* For |delta|<=width, ||D delta||^2 is at most
           width*abs(D'D)*width'.  This improves the sum-of-column-norms
           bound, especially for nearly orthogonal directions. */
        quadrad=sqrt(max((0,width*abs(D'*D)*width')))
        B.score_radius[i]=min((oldrad,quadrad))
    }
    B.support=J(rows(P.directions),2,.)
    for (pos=1; pos<=rows(P.directions); pos++) {
        B.support[pos,1]=jw_wbic_intersection_support(P.directions[pos,.],center,width)
        B.support[pos,2]=jw_wbic_intersection_support(-P.directions[pos,.],center,width)
    }
    return(B)
}

/* Noise event included in the coverage decomposition.  The outward
   tolerance only enlarges this event and therefore preserves the bound. */
real colvector jw_wbic_sharp_noise_event(struct jw_wbic_plan scalar P,
    struct jw_wbic_sharp_plan scalar S, real matrix E)
{
    real matrix X
    X=E[.,P.p:+S.active]*S.noiseW'
    return(rowsum(X:^2):<=P.radius^2+1e-10*(1+P.radius^2))
}

/* Reusable prepared version; no RNG or Stata state is accessed. */
real matrix jw_wbic_sharp_envelope_batch(struct jw_wbic_plan scalar P,
    struct jw_wbic_sharp_plan scalar S, struct jw_wbic_sharp_cell scalar B,
    real matrix E)
{
    real matrix lower,upper,possible,rootcap,baselines,tuplecaps,U,X,Eh
    real colvector normc,normg,lc,uc,lg,ug,viable,bestupper,v,uglob,ubox,
        covered,noise,base,caps,capq,lowerror,upperror,oldvalue,newvalue
    real rowvector idx,mu,dir
    real colvector owners,ii
    real scalar i,j,g,t,o,ns,nq,nb,r,rad,shift,supbox,pos
    ns=rows(P.masks)
    nq=rows(P.owners)
    nb=rows(E)
    lower=J(nb,nq,0)
    upper=J(nb,nq,0)
    mu=(P.h+P.L*B.center')'
    Eh=E[.,(P.p+1)..(P.p+P.k)]
    for (i=1; i<=nq; i++) {
        idx=P.scores[i].idx
        if (cols(idx)==0) continue
        X=(Eh[.,idx]+J(nb,1,1)*mu[idx])*P.scores[i].W'
        normc=sqrt(rowsum(X:^2))
        X=(Eh[.,idx]+J(nb,1,1)*P.h[idx]')*P.scores[i].W'
        normg=sqrt(rowsum(X:^2))
        rad=B.score_radius[i]
        rad=rad+1e-10*(1+rad)
        r=P.radius+1e-10*(1+P.radius)
        lc=rowmax((J(nb,1,0),normc:-rad)):^2
        uc=(normc:+rad):^2
        lg=rowmax((J(nb,1,0),normg:-r)):^2
        ug=(normg:+r):^2
        lower[.,i]=rowmax((lc,lg)):-P.scores[i].penalty
        upper[.,i]=rowmin((uc,ug)):-P.scores[i].penalty
    }
    if (hasmissing(lower) | hasmissing(upper)) jw_core_error("nonfinite WBIC score bounds")
    possible=J(nb,nq,0)
    rootcap=J(nb,nq,.)
    owners=uniqrows(sort(P.owners,1))
    for (g=1; g<=rows(owners); g++) {
        ii=selectindex(P.owners:==owners[g])
        bestupper=rowmin(upper[.,ii])
        for (j=1; j<=rows(ii); j++) {
            i=ii[j]
            possible[.,i]=lower[.,i]:<=bestupper:+1e-9:*(1:+abs(bestupper):+abs(lower[.,i]))
            /* Cache the winner's Wald-radius cap once per score and
               Gaussian draw, instead of once per candidate and target. */
            capq=bestupper:+P.scores[i].penalty
            capq=capq:+1e-9:*(1:+abs(bestupper):+P.scores[i].penalty)
            rootcap[.,i]=sqrt(rowmax((J(nb,1,0),capq)))
        }
    }
    U=J(nb,P.nout,0)
    covered=J(nb,1,0)
    noise=jw_wbic_sharp_noise_event(P,S,E)
    /* MD alters only precontrast coefficients.  The post-treatment
       loading is identical across candidates for each reporting target. */
    baselines=J(nb,P.nout,.)
    for (o=1; o<=P.nout; o++) {
        pos=(o-1)*ns+1
        baselines[.,o]=E[.,1..P.p]*P.actions[pos,1..P.p]'
    }
    for (t=1; t<=ns; t++) {
        viable=J(nb,1,1)
        for (g=1; g<=cols(P.tuples); g++) viable=viable:*possible[.,P.tuples[t,g]]
        if (sum(viable)==0) continue
        covered=covered:+viable
        tuplecaps=rootcap[.,P.tuples[t,.]]
        for (o=1; o<=P.nout; o++) {
            pos=(o-1)*ns+t
            dir=P.directions[pos,.]
            shift=P.shifts[t,o]+dir*B.center'
            supbox=abs(dir)*B.width'
            v=E*P.actions[pos,.]'
            /* Retain the exact legacy expression as a numerical cap. */
            ubox=abs(v:+shift):+supbox
            uglob=abs(v:+P.shifts[t,o]):+P.supports[t,o]
            oldvalue=rowmin((ubox,uglob)):/P.sigmas[t,o]
            oldvalue=oldvalue:+1e-10:*(1:+abs(oldvalue))
            lowerror=v:+P.shifts[t,o]:-B.support[pos,2]
            upperror=v:+P.shifts[t,o]:+B.support[pos,1]
            base=baselines[.,o]
            caps=tuplecaps*S.adjustment_scales[pos,.]'
            caps=caps:+1e-10:*(1:+abs(caps):+abs(base))
            lowerror=rowmax((lowerror,base-caps))
            upperror=rowmin((upperror,base+caps))
            /* A false candidate can have an empty signed intersection.
               Keeping both endpoints is harmless; actual winners have a
               nonempty intersection containing their realized error. */
            newvalue=rowmax((abs(lowerror),abs(upperror))):/P.sigmas[t,o]
            newvalue=newvalue:+1e-10:*(1:+abs(newvalue))
            v=rowmin((oldvalue,newvalue)):*viable:*noise
            if (hasmissing(v)) jw_core_error("nonfinite WBIC sharp error bound")
            U[.,o]=rowmax((U[.,o],v))
        }
    }
    if (min(covered)<1) jw_core_error("WBIC cell bounds retained no possible selection; check covariance conditioning")
    if (hasmissing(U)) jw_core_error("nonfinite WBIC sharp cell envelope")
    return(U)
}

/* Convenience API used by independent pathwise tests. */
real matrix jw_wbic_cell_envelope_sharp(struct jw_wbic_plan scalar P,
    real matrix E, real rowvector center, real rowvector width)
{
    struct jw_wbic_sharp_plan scalar S
    struct jw_wbic_sharp_cell scalar B
    S=jw_wbic_sharp_prepare(P)
    B=jw_wbic_sharp_cell_prepare(P,center,width)
    return(jw_wbic_sharp_envelope_batch(P,S,B,E))
}

struct jw_wbic_result scalar jw_infer_wbic_sharp(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps,
    real scalar ncells, real rowvector outputs, real scalar progress)
{
    struct jw_wbic_result scalar R
    struct jw_wbic_plan scalar P
    struct jw_wbic_sharp_plan scalar SP
    struct jw_wbic_sharp_cell scalar BC
    real matrix boxes,S,E,U,draws
    real rowvector active,rr,critical,pvalue,ranks
    real colvector ordered
    real scalar alpha,tail,delta,i,o,b,start,stop,nout,est,se,kq
    nout=rows(F.A)
    if (level<=0 | level>=100 | missing(level) | gamma<=0 | mcerr<=0 |
        missing(gamma) | missing(mcerr) | reps<1 | reps!=floor(reps) |
        missing(reps) | ncells<1 | ncells>1024 | ncells!=floor(ncells)) {
        jw_core_error("invalid WBIC inference probability or simulation options")
    }
    alpha=1-level/100
    if (gamma+mcerr>=alpha) jw_core_error("gamma()+mcerr() must be smaller than 1-level()/100")
    if (scope!="cohort" & scope!="joint") jw_core_error("unknown WBIC inference scope")
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs))) {
        jw_core_error("invalid WBIC inference output indices")
    }
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,8,.)
    for (o=1; o<=nout; o++) {
        R.table[o,1]=F.A[o,.]*M.z
        R.table[o,2]=sqrt(jw_fixed_variance(M.V,F.A[o,.]))
        R.diagnostics[o,.]=(o,.,.,.,0,.,.,1)
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        for (i=1; i<=cols(outputs); i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.diagnostics[o,.]=(o,0,1,0,0,invnormal(1-alpha/2),.,2)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    SP=jw_wbic_sharp_prepare(P)
    boxes=jw_wbic_cells(P.d,ncells)
    delta=mcerr/(rows(boxes)*cols(outputs))
    tail=1-alpha+gamma+mcerr
    ranks=jw_rs_quantile_ranks(reps,tail,delta)
    kq=ranks[2]
    if (kq>reps) jw_core_error("reps() is too small for WBIC level(), gamma(), mcerr(), and wbiccells(); increase reps()")
    critical=J(1,cols(outputs),0)
    pvalue=J(1,cols(outputs),0)
    S=jw_infer_sqrt(M.V)
    for (b=1; b<=rows(boxes); b++) {
        if (progress) printf("Sharp WBIC inference: cell %g of %g; %g joint candidates\n",b,rows(boxes),rows(P.masks))
        BC=jw_wbic_sharp_cell_prepare(P,boxes[b,1..P.d],boxes[b,(P.d+1)..(2*P.d)])
        draws=J(reps,cols(outputs),.)
        for (start=1; start<=reps; start=start+256) {
            stop=min((reps,start+255))
            E=rnormal(stop-start+1,rows(M.z),0,1)*S'
            U=jw_wbic_sharp_envelope_batch(P,SP,BC,E)
            draws[|start,1\stop,cols(outputs)|]=U
        }
        for (i=1; i<=cols(outputs); i++) {
            ordered=sort(draws[.,i],1)
            critical[i]=max((critical[i],ordered[kq]))
            o=outputs[i]
            est=abs(R.table[o,1])/R.table[o,2]
            pvalue[i]=max((pvalue[i],jw_rs_probability_upper(sum(draws[.,i]:>=est),reps,delta)))
        }
    }
    for (i=1; i<=cols(outputs); i++) {
        o=outputs[i]
        est=R.table[o,1]
        se=R.table[o,2]
        R.table[o,3]=est-critical[i]*se
        R.table[o,4]=est+critical[i]*se
        R.table[o,5]=min((1,gamma+mcerr+pvalue[i]))
        R.diagnostics[o,.]=(o,P.d,rows(P.masks),rows(boxes),reps,critical[i],kq,0)
    }
    return(R)
}

/* ===== jwselect_mccloskey.mata ===== */
/* Full outer size calibration of a fixed Gaussian quantile surface.
   Numerical mode validates searched points, not the nuisance continuum.
   Certified mode adds a genuine covering/total-variation bound and an
   analytic infinite-exterior bound.  Failure is explicit; the package
   wrapper may then use protected Berger--Boos inference, visibly marked.
   Inner nodes and all calibration decisions are independent of observed H.
   The realized inner Gaussian bank is conditioned on in outer calibration. */
struct jw_mc_result {
    real matrix table
    real matrix diagnostics
    real matrix calibration
    real scalar domain_radius
    real scalar exterior_radius
    real matrix bb_diagnostics
    real matrix cis_diagnostics
    real matrix cis_points
}

struct jw_mc_state {
    real rowvector active
    real matrix Hroot
    real matrix Hinv
    real matrix VH
    real matrix baseline
    real rowvector basevar
    real matrix cross
    real colvector owners
    real colvector rank
    real rowvector strides
    real colvector lookup
    real scalar boundary
    real scalar totalpenalty
}

struct jw_mc_draw {
    real colvector chosen
    real matrix statistic
    real matrix sigma
    real matrix adjustment
    real matrix standardized
}

struct jw_mc_state scalar jw_mc_prepare(struct jw_wbic_plan scalar P,
                                        real matrix V)
{
    struct jw_mc_state scalar S
    real rowvector ix,target,eig
    real colvector ii,code
    real scalar g,i,o,ns,mult,r,maxpen
    ns=rows(P.masks)
    S.active=selectindex(colsum(P.masks):>0)
    if (cols(S.active)!=P.d | P.d<1) jw_core_error("invalid McCloskey active dimensions")
    S.VH=V[P.p:+S.active,P.p:+S.active]
    eig=symeigenvalues((S.VH+S.VH')/2)
    if (hasmissing(eig) | min(eig)<=1e-12*max(abs(eig))) jw_core_error("singular McCloskey nuisance covariance")
    S.Hroot=jw_infer_sqrt(S.VH)
    S.Hinv=qrsolve(S.Hroot,I(P.d))
    if (hasmissing(S.Hinv)) jw_core_error("singular McCloskey nuisance covariance")
    S.baseline=J(P.nout,P.p+P.k,0)
    S.basevar=J(1,P.nout,.)
    for (o=1;o<=P.nout;o++) {
        S.baseline[o,1..P.p]=P.actions[(o-1)*ns+1,1..P.p]
        if (max(abs(P.actions[((o-1)*ns+1)..(o*ns),1..P.p]-
            J(ns,1,1)*S.baseline[o,1..P.p]))>
            1e-10*(1+max(abs(S.baseline[o,.])))) jw_core_error("McCloskey actions must share their baseline post loading")
        S.basevar[o]=jw_fixed_variance(V,S.baseline[o,.])
    }
    if (min(S.basevar)<=0 | hasmissing(S.basevar)) jw_core_error("degenerate McCloskey baseline variance")
    S.cross=V[P.p:+S.active,.]*S.baseline'
    S.owners=uniqrows(sort(P.owners,1))
    if (rows(S.owners)!=cols(P.tuples)) jw_core_error("incompatible McCloskey score families")
    S.rank=J(rows(P.owners),1,.)
    S.strides=J(1,rows(S.owners),.)
    mult=1
    S.boundary=0
    S.totalpenalty=0
    for (g=1;g<=rows(S.owners);g++) {
        ii=vec(selectindex(P.owners:==S.owners[g]))
        S.rank[ii]=(1::rows(ii))
        S.strides[g]=mult
        mult=mult*rows(ii)
        maxpen=0
        for (i=1;i<=rows(ii);i++) {
            r=cols(P.scores[ii[i]].idx)
            maxpen=max((maxpen,P.scores[ii[i]].penalty))
            if (r>0) S.boundary=max((S.boundary,sqrt(max((0,P.scores[ii[i]].penalty/r)))))
        }
        S.totalpenalty=S.totalpenalty+maxpen
    }
    if (mult>10000000) jw_core_error("McCloskey selection lookup is too large")
    S.lookup=J(mult,1,.)
    code=J(ns,1,1)
    for (g=1;g<=cols(P.tuples);g++) code=code+(S.rank[P.tuples[.,g]]:-1)*S.strides[g]
    for (i=1;i<=ns;i++) {
        if (!missing(S.lookup[code[i]])) jw_core_error("duplicate McCloskey candidate tuple")
        S.lookup[code[i]]=i
    }
    return(S)
}

/* Exact first-in-family WBIC decisions followed by a tuple lookup.
   No loop over the (possibly 32^3) union actions is needed per draw. */
real colvector jw_mc_choose(struct jw_wbic_plan scalar P,
                            struct jw_mc_state scalar S, real matrix H)
{
    real colvector code,best,choice,ii,score,changed
    real rowvector ix
    real matrix X
    real scalar g,j,q,n
    n=rows(H)
    code=J(n,1,1)
    for (g=1;g<=rows(S.owners);g++) {
        ii=vec(selectindex(P.owners:==S.owners[g]))
        best=J(n,1,.)
        choice=J(n,1,ii[1])
        for (j=1;j<=rows(ii);j++) {
            q=ii[j]
            ix=P.scores[q].idx
            if (cols(ix)==0) score=J(n,1,-P.scores[q].penalty)
            else {
                X=H[.,ix]*P.scores[q].W'
                score=rowsum(X:^2):-P.scores[q].penalty
            }
            changed=vec(selectindex(score:<best))
            if (rows(changed)>0) {
                best[changed]=score[changed]
                choice[changed]=J(rows(changed),1,q)
            }
        }
        code=code+(S.rank[choice]:-1)*S.strides[g]
    }
    choice=S.lookup[code]
    if (hasmissing(choice)) jw_core_error("selected tuple missing from McCloskey action family")
    return(choice)
}

struct jw_mc_draw scalar jw_mc_draws(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix E, real rowvector mu)
{
    struct jw_mc_draw scalar D
    real matrix H,base,A
    real scalar o,ns
    ns=rows(P.masks)
    H=E[.,(P.p+1)..(P.p+P.k)]+J(rows(E),1,1)*mu
    D.chosen=jw_mc_choose(P,S,H)
    D.standardized=H[.,S.active]*S.Hinv'
    D.statistic=J(rows(E),P.nout,.)
    D.sigma=J(rows(E),P.nout,.)
    D.adjustment=J(rows(E),P.nout,.)
    base=E*S.baseline'
    for (o=1;o<=P.nout;o++) {
        A=P.actions[(o-1)*ns:+D.chosen,(P.p+1)..(P.p+P.k)]
        D.adjustment[.,o]=rowsum(H:*A)
        D.sigma[.,o]=P.sigmas[D.chosen,o]
        D.statistic[.,o]=abs(base[.,o]+D.adjustment[.,o]):/D.sigma[.,o]
    }
    return(D)
}

/* Fixed design conditional on covariance, family, N and independent RNG.
   This routine deliberately accepts no observed pretrend argument. */
real matrix jw_mc_nodes(real scalar d, real scalar n, real scalar radius,
                         real scalar boundary)
{
    real matrix X
    real rowvector v
    real scalar i,j,sgn,len,r
    X=J(n,d,0)
    for (i=2;i<=n;i++) {
        if (i<=2*d+1) {
            j=ceil((i-1)/2)
            sgn=mod(i,2)==0 ? 1 : -1
            X[i,j]=sgn*min((radius,boundary))
        }
        else {
            v=rnormal(1,d,0,1)
            len=sqrt(sum(v:^2))
            if (len==0) v[1]=1
            len=sqrt(sum(v:^2))
            if (mod(i,4)==0) r=radius
            else if (mod(i,4)==1) r=min((radius,
                sqrt(1+mod(floor((i-2*d-2)/4),d))*boundary))
            else if (mod(i,4)==2) r=.5*radius
            else r=radius*runiform(1,1)^(1/d)
            X[i,.]=r*v/len
        }
    }
    return(X)
}

/* Actual random critical-value map; output-major quantile-level columns.
   Outside the numerical region / empty intersection uses the SAME
   baseline-containing safeguard in calibration and reported inference. */
real matrix jw_mc_critical(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, struct jw_mc_draw scalar D,
    real matrix nodes, real matrix quantiles, real scalar levels,
    real scalar radius, real scalar basecv, real colvector usedguard)
{
    real matrix K,X,Q,old
    real colvector inside,eligible,found,ix,guard
    real scalar j,o,n
    n=rows(D.standardized)
    K=J(n,cols(quantiles),0)
    inside=rowsum(D.standardized:^2):<=radius^2
    found=J(n,1,0)
    for (j=1;j<=rows(nodes);j++) {
        X=D.standardized-J(n,1,1)*nodes[j,.]
        eligible=inside:*(rowsum(X:^2):<=P.radius^2)
        ix=vec(selectindex(eligible))
        if (rows(ix)==0) continue
        found[ix]=J(rows(ix),1,1)
        Q=J(rows(ix),1,1)*quantiles[j,.]
        old=K[ix,.]
        K[ix,.]=old:*(old:>Q)+Q:*(old:<=Q)
    }
    guard=vec(selectindex(found:==0))
    usedguard=1:-found
    if (rows(guard)>0) {
        for (o=1;o<=P.nout;o++) {
            K[guard,((o-1)*levels+1)..(o*levels)]=
                ((abs(D.adjustment[guard,o]):+basecv*sqrt(S.basevar[o])):/
                D.sigma[guard,o])*J(1,levels,1)
        }
    }
    return(K)
}

real rowvector jw_mc_risk(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix E, real rowvector node,
    real matrix nodes, real matrix quantiles, real scalar levels,
    real scalar radius, real scalar basecv)
{
    struct jw_mc_draw scalar D
    real rowvector mu,counts
    real matrix K
    real colvector guard
    real scalar o,j,start,stop
    mu=J(1,P.k,0)
    mu[S.active]=node*S.Hroot'
    counts=J(1,P.nout*levels,0)
    for (start=1;start<=rows(E);start=start+256) {
        stop=min((rows(E),start+255))
        D=jw_mc_draws(P,S,E[start..stop,.],mu)
        K=jw_mc_critical(P,S,D,nodes,quantiles,levels,radius,basecv,guard)
        for (o=1;o<=P.nout;o++) {
            for (j=1;j<=levels;j++) {
                counts[(o-1)*levels+j]=counts[(o-1)*levels+j]+
                    sum(D.statistic[.,o]:>K[.,(o-1)*levels+j])
            }
        }
    }
    return(counts)
}

/* Gaussian TV uses the FULL mean-shift metric: d's null mean remains
   fixed, so conditioning out d changes the nuisance covariance. */
real matrix jw_mc_tv(struct jw_mc_state scalar S, real matrix boxes,
                      real scalar radius)
{
    real matrix TV,C,G
    real rowvector eig,w
    real scalar d,o,b,h
    d=cols(S.active)
    TV=J(rows(boxes),cols(S.basevar),1)
    for (o=1;o<=cols(S.basevar);o++) {
        C=S.VH-S.cross[.,o]*S.cross[.,o]'/S.basevar[o]
        eig=symeigenvalues((C+C')/2)
        if (hasmissing(eig) | min(eig)<=1e-12*max(abs(eig))) continue
        G=S.Hroot'*cholsolve(C,S.Hroot)
        G=(G+G')/2
        for (b=1;b<=rows(boxes);b++) {
            w=boxes[b,(d+1)..(2*d)]
            h=radius*sqrt(max((0,w*abs(G)*w')))
            h=h+1e-10*(1+h)
            TV[b,o]=min((1,2*normal(h/2)-1+1e-12))
        }
    }
    return(TV)
}

/* Generic API returns LOCAL output rows.  Status: 0 calibrated, 3 no
   feasible level, 4 cover budget insufficient.  No hidden BB fallback.
   Diagnostics: output nactive candidates inner_nodes inner_reps
   search_nodes outer_reps validation_nodes validation_reps pstar
   risk_upper tv_slack critical status certified rule_fallback.
   Calibration rows: level_index p, then for each output
   training_maxrisk validation_maxupper_including_TV max_TV_slack. */
struct jw_mc_result scalar jw_mc_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar gamma, real scalar mcerr, real scalar reps,
    string scalar mode, real scalar nodes, real scalar search,
    real scalar outerreps, real scalar validreps, real scalar steps,
    real scalar certcells, real scalar radius, real scalar pmax,
    real scalar progress)
{
    struct jw_mc_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    real matrix G,GP,GV,Q,E,ES,train,valid,boxes,TV,K
    real rowvector pgrid,mu,counts,maxtrain,maxvalid,maxtv,eligible
    real colvector sorted,worst,guard,choice,ix
    real scalar alpha,alpha0,beta,tau,B,outerradius,basecv,delta,n,o,j,i,c,
        ns,rank,k,est,se,upper,loc,selected
    if (mode!="numerical" & mode!="certified") jw_core_error("mcmode() must be numerical or certified")
    alpha=1-level/100
    beta=mcerr/2
    alpha0=alpha-beta
    if (hasmissing((level,gamma,mcerr,reps,nodes,search,outerreps,validreps,
        steps,certcells,radius,pmax)) | alpha<=0 | alpha>=1 | gamma<=0 | mcerr<=0 | gamma+mcerr>=alpha |
        reps<1 | reps!=floor(reps) | nodes<1 | nodes!=floor(nodes) |
        search<1 | search!=floor(search) | outerreps<1 | outerreps!=floor(outerreps) |
        validreps<1 | validreps!=floor(validreps) | steps<2 | steps!=floor(steps) |
        certcells<1 | certcells>65536 | certcells!=floor(certcells) |
        missing(radius) | radius<0) jw_core_error("invalid McCloskey calibration controls")
    if (rows(observedA)!=P.nout | cols(observedA)!=rows(z) |
        rows(V)!=rows(z) | cols(V)!=rows(z) | P.p+P.k!=rows(z) |
        P.nout<1 | cols(P.outputs)!=P.nout | cols(P.masks)!=P.k |
        rows(P.actions)!=rows(P.masks)*P.nout | cols(P.actions)!=rows(z) |
        rows(P.sigmas)!=rows(P.masks) | cols(P.sigmas)!=P.nout |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0 |
        hasmissing(V) | hasmissing(z) | hasmissing(observedA)) jw_core_error("incompatible McCloskey input dimensions")
    if (missing(P.radius) | abs(P.radius-sqrt(invchi2(P.d,1-gamma)))>
        1e-10*(1+P.radius)) jw_core_error("McCloskey plan radius disagrees with gamma()")
    S=jw_mc_prepare(P,V)
    ns=rows(P.masks)
    /* Exterior allocation is independent of the inner confidence-set gamma.
       A half-alpha guard leaves room for validation and covering error. */
    tau=alpha0/2
    B=radius
    if (B==0) B=sqrt(S.totalpenalty)+P.radius
    outerradius=B+sqrt(invchi2(P.d,1-tau))
    R.domain_radius=B
    R.exterior_radius=outerradius
    basecv=invnormal(1-(alpha0-tau)/2)
    if (pmax==0) pmax=1-alpha0+gamma
    if (pmax<1-alpha0 | pmax>=1 | missing(pmax)) jw_core_error("mcpmax() is outside the admissible probability range")
    pgrid=(1-alpha0):+((0..(steps-1))/(steps-1))*(pmax-(1-alpha0))
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,16,.)
    R.calibration=J(steps,2+3*P.nout,.)
    R.calibration[.,1]=(1::steps)
    R.calibration[.,2]=pgrid'
    for (o=1;o<=P.nout;o++) {
        est=observedA[o,.]*z
        se=sqrt(jw_fixed_variance(V,observedA[o,.]))
        R.table[o,1..2]=(est,se)
        R.diagnostics[o,.]=(P.outputs[o],P.d,ns,0,0,0,0,0,0,.,.,0,.,3,0,.)
    }
    eligible=J(1,P.nout,1)
    maxtv=J(1,P.nout,0)
    boxes=J(0,0,.)
    if (mode=="certified") {
        boxes=jw_wbic_cells(P.d,certcells)
        TV=jw_mc_tv(S,boxes,outerradius)
        maxtv=colmax(TV)
        eligible=maxtv:<alpha0
        for (o=1;o<=P.nout;o++) {
            R.diagnostics[o,8]=rows(boxes)
            R.diagnostics[o,12]=maxtv[o]
            R.calibration[.,2+3*o]=J(steps,1,maxtv[o])
            if (!eligible[o]) R.diagnostics[o,14]=4
        }
        if (sum(eligible)==0) return(R)
    }
    if (progress) printf("McCloskey: inner surface, %g nodes and %g Gaussian draws\n",nodes,reps)
    G=jw_mc_nodes(P.d,nodes,B,max((S.boundary,1e-8)))
    ES=jw_infer_sqrt(V)
    E=rnormal(reps,rows(z),0,1)*ES'
    Q=J(nodes,P.nout*steps,.)
    for (i=1;i<=nodes;i++) {
        mu=J(1,P.k,0)
        mu[S.active]=G[i,.]*S.Hroot'
        D=jw_mc_draws(P,S,E,mu)
        for (o=1;o<=P.nout;o++) {
            sorted=sort(D.statistic[.,o],1)
            for (j=1;j<=steps;j++) {
                rank=max((1,min((reps,ceil(reps*pgrid[j])))))
                Q[i,(o-1)*steps+j]=sorted[rank]
            }
        }
    }
    maxtrain=J(1,P.nout*steps,.)
    if (mode=="numerical") {
        GP=G \ jw_mc_nodes(P.d,search,outerradius,max((S.boundary,1e-8)))
        GP=uniqrows(sort(GP,1..P.d))
        if (progress) printf("McCloskey: pilot size search at %g nuisance nodes\n",rows(GP))
        E=rnormal(outerreps,rows(z),0,1)*ES'
        train=J(rows(GP),P.nout*steps,.)
        for (i=1;i<=rows(GP);i++) train[i,.]=jw_mc_risk(P,S,E,GP[i,.],G,Q,steps,B,basecv)/outerreps
        maxtrain=colmax(train)
        worst=J(P.nout*steps,1,.)
        for (j=1;j<=cols(train);j++) {
            ix=vec(selectindex(train[.,j]:==maxtrain[j]))
            worst[j]=ix[1]
        }
        GV=GP[worst,.] \ jw_mc_nodes(P.d,search,outerradius,max((S.boundary,1e-8)))
        GV=uniqrows(sort(GV,1..P.d))
        TV=J(rows(GV),P.nout,0)
    }
    else {
        GP=J(0,P.d,.)
        GV=outerradius*boxes[.,1..P.d]
    }
    delta=beta/(rows(GV)*P.nout*steps)
    if (progress) printf("McCloskey: independent validation at %g nodes, %g draws each\n",rows(GV),validreps)
    E=rnormal(validreps,rows(z),0,1)*ES'
    maxvalid=J(1,P.nout*steps,0)
    for (i=1;i<=rows(GV);i++) {
        counts=jw_mc_risk(P,S,E,GV[i,.],G,Q,steps,B,basecv)
        for (o=1;o<=P.nout;o++) {
            for (j=1;j<=steps;j++) {
                c=(o-1)*steps+j
                upper=jw_rs_probability_upper(counts[c],validreps,delta)+TV[i,o]
                maxvalid[c]=max((maxvalid[c],min((1,upper))))
            }
        }
    }
    /* Observed H enters only here, after the surface and size calibration. */
    E=J(1,rows(z),0)
    D=jw_mc_draws(P,S,E,z[(P.p+1)..rows(z)]')
    K=jw_mc_critical(P,S,D,G,Q,steps,B,basecv,guard)
    for (o=1;o<=P.nout;o++) {
        loc=(o-1)*steps
        R.calibration[.,3*o]=maxtrain[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+1]=maxvalid[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+2]=J(steps,1,maxtv[o])
        R.diagnostics[o,4..9]=(nodes,reps,rows(GP),(mode=="numerical" ? outerreps : 0),rows(GV),validreps)
        R.diagnostics[o,11]=maxvalid[loc+steps]
        R.diagnostics[o,12]=maxtv[o]
        R.diagnostics[o,16]=guard[1]
        if (!eligible[o]) continue
        selected=0
        for (j=1;j<=steps;j++) {
            if (maxvalid[loc+j]<=alpha0) {
                selected=j
                break
            }
        }
        if (!selected) continue
        if (max(abs(P.actions[(o-1)*ns+D.chosen[1],.]-observedA[o,.]))>
            1e-8*(1+max(abs(observedA[o,.])))) jw_core_error("observed McCloskey action differs from the actual WBIC selector")
        R.diagnostics[o,10]=pgrid[selected]
        R.diagnostics[o,11]=maxvalid[loc+selected]
        R.diagnostics[o,13]=K[1,loc+selected]
        R.diagnostics[o,14]=0
        R.diagnostics[o,15]=mode=="certified"
        est=R.table[o,1]
        se=R.table[o,2]
        R.table[o,3..4]=(est-K[1,loc+selected]*se,est+K[1,loc+selected]*se)
    }
    return(R)
}

/* Package API: transport local generic results into all reporting rows.
   A failed calibration can use a FRESH protected BB calculation at
   alpha0=alpha-mcerr/2; this is explicitly status5 or6, never status0. */
struct jw_mc_result scalar jw_infer_mc(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps, real scalar ncells,
    real rowvector outputs, real scalar progress, string scalar mode,
    real scalar nodes, real scalar search, real scalar outerreps,
    real scalar validreps, real scalar steps, real scalar certcells,
    real scalar radius, real scalar pmax, string scalar fallback)
{
    struct jw_mc_result scalar R,Q
    struct jw_wbic_result scalar BB
    struct jw_wbic_plan scalar P
    real rowvector active,failed
    real scalar i,o,nout,alpha0,status
    if (fallback!="bb" & fallback!="error") jw_core_error("mcfallback() must be bb or error")
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs))) jw_core_error("invalid McCloskey output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.calibration=J(0,0,.)
    R.domain_radius=0
    R.exterior_radius=0
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,14]=1
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,2..9]=(0,1,0,0,0,0,0,0)
            R.diagnostics[o,13..16]=(invnormal(1-(1-level/100)/2),2,1,0)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    Q=jw_mc_from_plan(P,M.V,M.z,F.A[outputs,.],level,gamma,mcerr,reps,
        mode,nodes,search,outerreps,validreps,steps,certcells,radius,pmax,progress)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    R.calibration=Q.calibration
    R.domain_radius=Q.domain_radius
    R.exterior_radius=Q.exterior_radius
    failed=J(1,0,.)
    for (i=1;i<=cols(outputs);i++) {
        if (Q.diagnostics[i,14]>=3) failed=(failed,outputs[i])
    }
    if (cols(failed)>0) {
        if (fallback=="error") jw_core_error("McCloskey calibration did not validate; increase its budget or use mcfallback(bb)")
        alpha0=1-level/100-mcerr/2
        if (progress) printf("McCloskey: protected BB fallback for %g outputs; this is NOT calibrated-M success\n",cols(failed))
        BB=jw_infer_wbic_sharp(M,C,F,joint,scope,maxjoint,100*(1-alpha0),
            gamma,mcerr/2,reps,ncells,failed,progress)
        for (i=1;i<=cols(failed);i++) {
            o=failed[i]
            status=R.diagnostics[o,14]
            R.table[o,.]=BB.table[o,.]
            R.table[o,5]=.
            R.diagnostics[o,13]=BB.diagnostics[o,6]
            R.diagnostics[o,14]=status==4 ? 6 : 5
            R.diagnostics[o,15]=0
        }
    }
    return(R)
}

/* ===== jwselect_bb.mata ===== */
/* Fixed-nuisance Berger--Boos quantiles, version 0.8.3.
   The target is sup_{mu in C(H)} q_(1-alpha+gamma+mcerr)(|T(mu)|).
   Each draw reselects WBIC. Nuisance maximization is OUTSIDE the quantile.

   Certified mode uses fresh independent normal draws at each adaptively
   chosen feasible center. Three one-sided order-statistic events per
   center/output receive a fixed union-bound allocation. A Gaussian
   Neyman--Pearson shift bound extends each center bound over its cell.
   An analytic normal union bound supplies a finite, visibly flagged cap.
   Budget exhaustion returns a conservative upper bound, NOT a claim that
   the sharp extremum has been found. Floating-point probability routines
   are assumed accurate; this is not an interval-arithmetic certificate.

   Numerical scalar mode evaluates the actual selected Gaussian law by
   bivariate-normal integration. Other numerical cases use pointwise
   Monte Carlo upper quantiles. Neither numerical mode is certified. */

struct jw_bb_result {
    real matrix table
    real matrix diagnostics
    real matrix cis_diagnostics
    real matrix cis_points
}

struct jw_bb_cell {
    real rowvector box
    real rowvector point
    real rowvector lower
    real rowvector pointupper
    real rowvector pointcap
    real rowvector upper
    real rowvector radius
    real rowvector cap
    real scalar alive
}

/* A feasible point in box intersect unit ball. Radial projection of the
   box center alone need not stay in the box, hence the segment below. */
real rowvector jw_bb_feasible(real rowvector box)
{
    real scalar d,j,lo,hi,t,it
    real rowvector c,w,a,x
    d=cols(box)/2
    c=box[1..d]
    w=box[(d+1)..(2*d)]
    if (sum(c:^2)<=1) return(c)
    a=J(1,d,0)
    for (j=1;j<=d;j++) {
        if (c[j]>w[j]) a[j]=c[j]-w[j]
        else if (c[j]<-w[j]) a[j]=c[j]+w[j]
    }
    if (sum(a:^2)>1) {
        /* Only an outward-rounded boundary cell can reach this branch.
           Projection is still feasible; its increased cell distance is
           explicitly included in jw_bb_radius(). */
        return(a/sqrt(sum(a:^2)))
    }
    lo=0
    hi=1
    for (it=1;it<=55;it++) {
        t=(lo+hi)/2
        x=a+t*(c-a)
        if (sum(x:^2)<=1) lo=t
        else hi=t
    }
    return(a+lo*(c-a))
}

real scalar jw_bb_intersects(real rowvector box)
{
    real scalar d,j,s
    d=cols(box)/2
    s=0
    for (j=1;j<=d;j++) s=s+max((0,abs(box[j])-box[d+j]))^2
    return(s<=1+1e-12)
}

/* Rows contain vec(G_o)' where G_o=L' Cov(H|baseline_o)^(-1) L.
   Missing rows signal a singular full-joint location metric. */
real matrix jw_bb_metric(struct jw_wbic_plan scalar P,
                         struct jw_mc_state scalar S)
{
    real matrix out,C,G,L
    real rowvector eig
    real scalar o
    out=J(P.nout,P.d^2,.)
    L=P.L[S.active,.]
    for (o=1;o<=P.nout;o++) {
        C=S.VH-S.cross[.,o]*S.cross[.,o]'/S.basevar[o]
        C=(C+C')/2
        eig=symeigenvalues(C)
        if (hasmissing(eig) | min(eig)<=1e-12*max(abs(eig))) continue
        G=L'*cholsolve(C,L)
        G=(G+G')/2
        out[o,.]=vec(G)'
    }
    return(out)
}

real rowvector jw_bb_radius(real matrix metric, real rowvector box,
                            real rowvector point)
{
    real scalar d,o,r
    real matrix G
    real rowvector w,out
    d=cols(point)
    w=box[(d+1)..(2*d)]+abs(box[1..d]-point)
    out=J(1,rows(metric),.)
    for (o=1;o<=rows(metric);o++) {
        if (hasmissing(metric[o,.])) continue
        G=colshape(metric[o,.],d)
        r=sqrt(max((0,w*abs(G)*w')))
        out[o]=r+1e-10*(1+r)
    }
    return(out)
}

/* Analytic cap for the actual selected law: union bound over the finite
   candidate normal errors, with exact affine bias support on a box and
   a unit-ball support cap. This cap is not substituted for the target. */
real rowvector jw_bb_normal_cap(struct jw_wbic_plan scalar P,
                                real rowvector box, real scalar tail)
{
    real scalar ns,o,t,i,z,b1,b2,b
    real rowvector out,c,w,direction
    ns=rows(P.masks)
    c=box[1..P.d]
    w=box[(P.d+1)..(2*P.d)]
    z=invnormal(1-tail/(2*ns))
    if (missing(z)) z=sqrt(2*ln(2*ns/tail))
    out=J(1,P.nout,0)
    for (o=1;o<=P.nout;o++) {
        for (t=1;t<=ns;t++) {
            i=(o-1)*ns+t
            direction=P.actions[i,(P.p+1)..(P.p+P.k)]*P.L
            b=P.actions[i,(P.p+1)..(P.p+P.k)]*P.h
            b1=abs(b+direction*c')+sum(abs(direction):*w)
            b2=abs(b)+sqrt(sum(direction:^2))
            out[o]=max((out[o],min((b1,b2))/P.sigmas[t,o]+z))
        }
    }
    return(out:+1e-10*(1:+out))
}

/* Scalar selection regions. All WBIC scores are affine functions of H^2;
   every pairwise crossing is included before identifying each winner. */
real matrix jw_bb_scalar_regions(struct jw_wbic_plan scalar P,
                                 struct jw_mc_state scalar S)
{
    real colvector roots,bounds,ix,chosen
    real matrix regions,H
    real rowvector a
    real scalar g,i,j,q,t,x,mid,n
    roots=J(0,1,.)
    for (g=1;g<=rows(S.owners);g++) {
        ix=vec(selectindex(P.owners:==S.owners[g]))
        a=J(1,rows(ix),0)
        for (i=1;i<=rows(ix);i++) {
            q=ix[i]
            if (cols(P.scores[q].idx)>0) a[i]=sum(P.scores[q].W:^2)
        }
        for (i=1;i<=rows(ix);i++) for (j=i+1;j<=rows(ix);j++) {
            if (a[i]==a[j]) continue
            x=(P.scores[ix[i]].penalty-P.scores[ix[j]].penalty)/(a[i]-a[j])
            if (x>0 & !missing(x)) roots=roots\sqrt(x)
        }
    }
    if (rows(roots)>0) roots=uniqrows(sort(roots,1))
    bounds=(-8e307)\(-roots)\roots\(8e307)
    bounds=uniqrows(sort(bounds,1))
    n=rows(bounds)-1
    regions=J(n,3,.)
    for (t=1;t<=n;t++) {
        if (t==1) mid=bounds[t+1]-max((1,abs(bounds[t+1])))
        else if (t==n) mid=bounds[t]+max((1,abs(bounds[t])))
        else mid=(bounds[t]+bounds[t+1])/2
        if (n==1) mid=0
        H=J(1,P.k,0)
        H[S.active]=mid
        chosen=jw_mc_choose(P,S,H)
        regions[t,.]=(bounds[t],bounds[t+1],chosen[1])
    }
    return(regions)
}

real scalar jw_bb_bvn(real scalar x, real scalar y, real scalar rho)
{
    if (hasmissing((x,y,rho))) jw_core_error("nonfinite scalar Gaussian integration argument")
    if (x<=-38 | y<=-38) return(0)
    if (x>=38) return(normal(y))
    if (y>=38) return(normal(x))
    if (rho>=1-1e-14) return(normal(min((x,y))))
    if (rho<=-1+1e-14) return(max((0,normal(x)-normal(-y))))
    return(binormal(x,y,rho))
}

/* Actual scalar selected-law tail, analytically reduced to Gaussian
   rectangle probabilities. Exposed for independent benchmark tests.
   It is used only in numerical mode: no quadrature certificate is claimed. */
real scalar jw_bb_scalar_tail(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix regions, real scalar mu,
    real scalar critical, real scalar output)
{
    real scalar t,k,ns,ah,sigma,sd,rho,lo,hi,yl,yu,p
    ns=rows(P.masks)
    sd=sqrt(S.VH[1,1])
    p=0
    for (t=1;t<=rows(regions);t++) {
        k=regions[t,3]
        ah=P.actions[(output-1)*ns+k,P.p+S.active[1]]
        sigma=P.sigmas[k,output]
        rho=(S.cross[1,output]+ah*S.VH[1,1])/(sd*sigma)
        rho=max((-1,min((1,rho))))
        /* Do not divide the infinity sentinels: small sd would overflow
           to Mata missing, which compares as a large POSITIVE number. */
        lo=regions[t,1]<-1e300 ? -40 : (regions[t,1]-mu)/sd
        hi=regions[t,2]>1e300 ? 40 : (regions[t,2]-mu)/sd
        yl=-critical-ah*mu/sigma
        yu=critical-ah*mu/sigma
        p=p+max((0,jw_bb_bvn(hi,yl,rho)-jw_bb_bvn(lo,yl,rho)))+
            max((0,jw_bb_bvn(-lo,-yu,rho)-jw_bb_bvn(-hi,-yu,rho)))
    }
    return(min((1,max((0,p)))))
}

real rowvector jw_bb_scalar_quantile(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix regions, real scalar mu,
    real scalar tail, real scalar output, real scalar cap)
{
    real scalar lo,hi,mid,it
    lo=0
    hi=cap
    if (tail<=0 | missing(tail)) return((0,.))
    while (jw_bb_scalar_tail(P,S,regions,mu,hi,output)>tail) {
        hi=2*hi+1
        if (hi>1e12) return((0,.))
    }
    for (it=1;it<=65;it++) {
        mid=(lo+hi)/2
        if (jw_bb_scalar_tail(P,S,regions,mu,mid,output)>tail) lo=mid
        else hi=mid
        if (hi-lo<=1e-9*(1+hi)) break
    }
    return((lo,hi))
}

struct jw_bb_cell scalar jw_bb_evaluate(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix root, real matrix metric,
    real matrix regions, real rowvector box, real scalar tail,
    real scalar delta, real scalar reps, real scalar scalaranalytic)
{
    struct jw_bb_cell scalar B
    struct jw_mc_draw scalar D
    real matrix E
    real rowvector mu,cap,rank,pointbox,qbounds
    real colvector sorted
    real scalar o,q,n
    B.box=box
    B.point=jw_bb_feasible(box)
    B.radius=jw_bb_radius(metric,box,B.point)
    B.alive=1
    B.lower=J(1,P.nout,0)
    B.pointupper=J(1,P.nout,.)
    B.pointcap=J(1,P.nout,0)
    B.upper=jw_bb_normal_cap(P,box,tail)
    B.cap=J(1,P.nout,1)
    mu=P.h'+(P.L*B.point')'
    pointbox=(B.point,J(1,P.d,0))
    cap=jw_bb_normal_cap(P,pointbox,tail)
    if (!scalaranalytic) {
        E=rnormal(reps,cols(root),0,1)*root'
        D=jw_mc_draws(P,S,E,mu)
    }
    for (o=1;o<=P.nout;o++) {
        if (scalaranalytic) {
            qbounds=jw_bb_scalar_quantile(P,S,regions,mu[S.active[1]],tail,o,cap[o])
            B.lower[o]=qbounds[1]
            B.pointupper[o]=qbounds[2]
        }
        else {
            sorted=sort(D.statistic[.,o],1)
            rank=jw_rs_quantile_ranks(reps,1-tail,delta)
            if (rank[1]>0) B.lower[o]=sorted[rank[1]]
            B.pointupper[o]=cap[o]
            B.pointcap[o]=1
            if (rank[2]<=reps) {
                if (sorted[rank[2]]<cap[o]) {
                    B.pointupper[o]=sorted[rank[2]]
                    B.pointcap[o]=0
                }
            }
        }
        if (missing(B.radius[o])) continue
        q=normal(invnormal(tail)-B.radius[o])
        /* A tiny q cannot be resolved by finite order statistics. */
        if (q<=0 | 1-q>=1) continue
        if (scalaranalytic) {
            qbounds=jw_bb_scalar_quantile(P,S,regions,mu[S.active[1]],q,o,B.upper[o])
            n=qbounds[2]
        }
        else {
            rank=jw_rs_quantile_ranks(reps,1-q,delta)
            n=.
            if (rank[2]<=reps) n=sorted[rank[2]]
        }
        if (!missing(n) & n<B.upper[o]) {
            B.upper[o]=n+1e-10*(1+n)
            B.cap[o]=0
        }
    }
    return(B)
}

/* Diagnostics columns: output nactive candidates nodes draws target_tail
   critical_lo critical_hi gap critical status certified accuracy_certified
   cap_used max_radius fixed_point_cv. Status: 0 tolerance met; 1 numerical;
   2 no-restriction normal; 4 budget exhausted; 9 not requested. */
struct jw_bb_result scalar jw_bb_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar gamma, real scalar mcerr, real scalar reps,
    string scalar mode, real scalar cells, real scalar tol, real scalar progress)
{
    struct jw_bb_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    struct jw_bb_cell colvector B
    real matrix root,metric,regions,G
    real rowvector lower,upper,pointupper,pointcap,gap,cap,radius,box,child,w,eig
    real colvector live
    real scalar alpha,tail,maxnodes,delta,n,i,o,q,j,dim,score,best,ns,
        critical,est,se,scalaranalytic,alive,done,weight,stat
    if (mode!="numerical" & mode!="certified") jw_core_error("bbmode() must be numerical or certified")
    alpha=1-level/100
    tail=alpha-gamma-mcerr
    if (hasmissing((level,gamma,mcerr,reps,cells,tol)) | alpha<=0 | alpha>=1 |
        gamma<=0 | mcerr<=0 | tail<=0 | reps<1 | reps!=floor(reps) |
        cells<1 | cells>65536 | cells!=floor(cells) | tol<=0) jw_core_error("invalid BB controls")
    if (P.d<1 | rows(observedA)!=P.nout | cols(observedA)!=rows(z) |
        rows(V)!=rows(z) | cols(V)!=rows(z) | P.p+P.k!=rows(z) |
        P.nout<1 | cols(P.outputs)!=P.nout | cols(P.masks)!=P.k |
        rows(P.actions)!=rows(P.masks)*P.nout | cols(P.actions)!=rows(z) |
        rows(P.sigmas)!=rows(P.masks) | cols(P.sigmas)!=P.nout |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0 |
        hasmissing(V) | hasmissing(z) | hasmissing(observedA)) jw_core_error("incompatible BB input dimensions")
    if (missing(P.radius) | abs(P.radius-sqrt(invchi2(P.d,1-gamma)))>
        1e-10*(1+P.radius)) jw_core_error("BB plan radius disagrees with gamma()")
    S=jw_mc_prepare(P,V)
    ns=rows(P.masks)
    root=jw_infer_sqrt(V)
    metric=jw_bb_metric(P,S)
    scalaranalytic=(P.d==1 & mode=="numerical")
    regions=J(0,3,.)
    if (scalaranalytic) regions=jw_bb_scalar_regions(P,S)
    maxnodes=2*cells-1
    /* Fixed maximum and fresh draws make this valid under adaptive splits. */
    delta=mcerr/(3*maxnodes*P.nout)
    B=jw_bb_cell(maxnodes)
    n=1
    B[1]=jw_bb_evaluate(P,S,root,metric,regions,
        (J(1,P.d,0),J(1,P.d,1)),tail,delta,reps,scalaranalytic)
    lower=B[1].lower
    pointupper=B[1].pointupper
    pointcap=B[1].pointcap
    if (progress) printf("BB: actual selected-law quantiles, at most %g adaptive nodes\n",maxnodes)
    while (1) {
        upper=J(1,P.nout,0)
        cap=J(1,P.nout,0)
        radius=J(1,P.nout,0)
        alive=0
        for (i=1;i<=n;i++) if (B[i].alive) {
            alive=alive+1
            for (o=1;o<=P.nout;o++) {
                if (B[i].upper[o]>upper[o]) {
                    upper[o]=B[i].upper[o]
                    cap[o]=B[i].cap[o]
                }
                radius[o]=max((radius[o],B[i].radius[o]))
            }
        }
        gap=upper-lower
        done=all(gap:<=tol)
        if (done | n+2>maxnodes | alive>=cells) break
        q=0
        best=-1
        for (i=1;i<=n;i++) if (B[i].alive) {
            score=max(B[i].upper-lower)
            if (score>best) {
                best=score
                q=i
            }
        }
        box=B[q].box
        dim=1
        best=-1
        for (j=1;j<=P.d;j++) {
            weight=1
            for (o=1;o<=P.nout;o++) if (!hasmissing(metric[o,.])) {
                G=colshape(metric[o,.],P.d)
                weight=max((weight,G[j,j]))
            }
            score=box[P.d+j]^2*weight
            if (score>best) {
                best=score
                dim=j
            }
        }
        B[q].alive=0
        for (j=1;j<=2;j++) {
            child=box
            child[P.d+dim]=box[P.d+dim]/2
            child[dim]=box[dim]+(2*j-3)*child[P.d+dim]
            if (!jw_bb_intersects(child)) continue
            n=n+1
            B[n]=jw_bb_evaluate(P,S,root,metric,regions,child,tail,delta,reps,scalaranalytic)
            for (o=1;o<=P.nout;o++) {
                lower[o]=max((lower[o],B[n].lower[o]))
                if (B[n].pointupper[o]>pointupper[o]) {
                    pointupper[o]=B[n].pointupper[o]
                    pointcap[o]=B[n].pointcap[o]
                }
                /* The parent bound remains valid on either child. */
                if (B[q].upper[o]<B[n].upper[o]) {
                    B[n].upper[o]=B[q].upper[o]
                    B[n].cap[o]=B[q].cap[o]
                }
            }
        }
        if (progress & mod(n,31)==0) printf("BB: %g nodes, current maximum critical bracket gap %9.5g\n",n,max(gap))
    }
    /* Validate correspondence between reported and reselected action. */
    D=jw_mc_draws(P,S,J(1,rows(z),0),z[(P.p+1)..rows(z)]')
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,16,.)
    for (o=1;o<=P.nout;o++) {
        if (max(abs(P.actions[(o-1)*ns+D.chosen[1],.]-observedA[o,.]))>
            1e-8*(1+max(abs(observedA[o,.])))) jw_core_error("observed BB action differs from actual WBIC selector")
        est=observedA[o,.]*z
        se=sqrt(jw_fixed_variance(V,observedA[o,.]))
        critical=mode=="certified" ? upper[o] : min((pointupper[o],upper[o]))
        stat=mode=="numerical" ? 1 : (gap[o]<=tol ? 0 : 4)
        R.table[o,.]=(est,se,est-critical*se,est+critical*se,.)
        R.diagnostics[o,.]=(P.outputs[o],P.d,ns,n,
            (scalaranalytic ? 0 : n*reps),tail,lower[o],upper[o],gap[o],critical,
            stat,(mode=="certified"),(mode=="certified" & gap[o]<=tol),
            (mode=="certified" ? cap[o] : (pointupper[o]<=upper[o] ? pointcap[o] : cap[o])),
            radius[o],pointupper[o])
    }
    return(R)
}

struct jw_bb_result scalar jw_infer_bb(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps, real scalar cells,
    real scalar tol, real rowvector outputs, real scalar progress, string scalar mode)
{
    struct jw_bb_result scalar R,Q
    struct jw_wbic_plan scalar P
    real rowvector active
    real scalar o,i,nout,cv,tail
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout |
        any(outputs:!=floor(outputs))) jw_core_error("invalid BB output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,11]=9
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        cv=invnormal(1-(1-level/100)/2)
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,.]=(o,0,1,0,0,1-level/100,cv,cv,0,cv,2,1,1,0,0,cv)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    Q=jw_bb_from_plan(P,M.V,M.z,F.A[outputs,.],level,gamma,mcerr,reps,mode,cells,tol,progress)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    return(R)
}

/* ===== jwselect_mc_refined.mata ===== */
/* McCloskey adjusted calibration for the actual Gaussian WBIC rule.
   This is a finite-Gaussian adaptation, not the published 2020 localized
   asymptotic algorithm.  Every fixed candidate has a scalar bias interval;
   its inner folded-normal quantile supremum is computed at the endpoints.
   The random candidate and random interval are accounted for JOINTLY by
   outer calibration.  mcbeta is independent of the BB fallback gamma.
   All tuning and validation precede use of the observed precontrasts.
   Dependencies: jwselect_mccloskey.mata and jwselect_bb.mata. */

/* CDF of |N(mean,1)|.  The expression avoids subtracting two near-one CDFs. */
real scalar jw_mcr_foldcdf(real scalar x, real scalar mean)
{
    if (x<=0) return(0)
    mean=abs(mean)
    return(max((0,min((1,normal(x-mean)-normal(-x-mean))))))
}

/* Return the upper endpoint of a deterministic quantile bracket.
   The interval [0, |mean|+z_(1+p)/2] contains the desired quantile. */
real scalar jw_mcr_foldquantile(real scalar mean, real scalar sigma,
                                  real scalar p)
{
    real scalar lo,hi,mid,m,i
    if (hasmissing((mean,sigma,p)) | sigma<=0 | p<=0 | p>=1)
        jw_core_error("invalid folded-normal quantile inputs")
    m=abs(mean)/sigma
    lo=0
    hi=m+invnormal((1+p)/2)
    for (i=1;i<=80;i++) {
        mid=(lo+hi)/2
        if (jw_mcr_foldcdf(mid,m)<p) lo=mid
        else hi=mid
        if (hi-lo<=2e-12*(1+hi)) break
    }
    return(sigma*hi)
}

real matrix jw_mcr_bias_sd(struct jw_wbic_plan scalar P, real matrix V)
{
    real matrix B,VH
    real rowvector a
    real scalar ns,o,s,v
    ns=rows(P.masks)
    VH=V[(P.p+1)..(P.p+P.k),(P.p+1)..(P.p+P.k)]
    B=J(ns,P.nout,0)
    for (o=1;o<=P.nout;o++) {
        for (s=1;s<=ns;s++) {
            a=P.actions[(o-1)*ns+s,(P.p+1)..(P.p+P.k)]
            v=a*VH*a'
            if (missing(v) | v< -1e-10*(1+max(abs(VH))))
                jw_core_error("invalid candidate bias variance")
            B[s,o]=sqrt(max((0,v)))
        }
    }
    return(B)
}

/* Sharp Gaussian probability transport for a fixed event and mean shift.
   If the reference probability is at most u, its shifted probability is
   at most Phi(Phi^-1(u)+distance). */
real scalar jw_mcr_np_upper(real scalar u, real scalar distance)
{
    if (missing(u) | missing(distance) | distance<0) return(1)
    if (u>=1) return(1)
    if (u<=0) return(0)
    return(min((1,normal(invnormal(u)+distance)+1e-12)))
}

/* FULL (baseline,H) Gaussian shift metric, separately for every target.
   A cell is radius*(center +/- width); absolute Gram entries bound the
   maximum quadratic norm throughout its box, hence also its ball slice. */
real matrix jw_mcr_cell_distance(struct jw_mc_state scalar S,
                                  real matrix boxes, real scalar radius)
{
    real matrix D,C,G
    real rowvector eig,w
    real scalar d,o,b,h
    d=cols(S.active)
    D=J(rows(boxes),cols(S.basevar),.)
    for (o=1;o<=cols(S.basevar);o++) {
        C=S.VH-S.cross[.,o]*S.cross[.,o]'/S.basevar[o]
        eig=symeigenvalues((C+C')/2)
        if (hasmissing(eig) | min(eig)<=1e-12*max(abs(eig))) continue
        G=S.Hroot'*cholsolve(C,S.Hroot)
        G=(G+G')/2
        for (b=1;b<=rows(boxes);b++) {
            w=boxes[b,(d+1)..(2*d)]
            h=radius*sqrt(max((0,w*abs(G)*w')))
            D[b,o]=h+1e-10*(1+h)
        }
    }
    return(D)
}

/* Exact replay of the WHOLE implemented rule.  For observations inside B,
   T>q_p(|bhat|+z_beta*tau_s,sigma_s) iff its folded-normal CDF exceeds p.
   Thus outer simulation needs no quantile approximation or inner draws. */
real rowvector jw_mcr_risk(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix biasSD, real matrix E,
    real rowvector node, real rowvector pgrid, real scalar B,
    real scalar basecv, real scalar zbeta)
{
    struct jw_mc_draw scalar D
    real rowvector mu,counts
    real colvector inside,bias,pivot,guardreject,ix
    real scalar o,j,start,stop,levels
    levels=cols(pgrid)
    mu=J(1,P.k,0)
    mu[S.active]=node*S.Hroot'
    counts=J(1,P.nout*levels,0)
    for (start=1;start<=rows(E);start=start+256) {
        stop=min((rows(E),start+255))
        D=jw_mc_draws(P,S,E[start..stop,.],mu)
        inside=rowsum(D.standardized:^2):<=B^2
        for (o=1;o<=P.nout;o++) {
            bias=(abs(D.adjustment[.,o]):+zbeta*biasSD[D.chosen,o]):/D.sigma[.,o]
            pivot=normal(D.statistic[.,o]-bias)-normal(-D.statistic[.,o]-bias)
            guardreject=D.statistic[.,o]:>
                ((abs(D.adjustment[.,o]):+basecv*sqrt(S.basevar[o])):/D.sigma[.,o])
            for (j=1;j<=levels;j++) {
                counts[(o-1)*levels+j]=counts[(o-1)*levels+j]+
                    sum(inside:*(pivot:>pgrid[j])+(1:-inside):*guardreject)
            }
        }
    }
    return(counts)
}

/* Same result layout as the retained legacy MC engine.  Columns4/5 are
   zero: there are no inner nuisance nodes or inner simulation draws.
   Column12 is now the maximum Gaussian mean-shift radius, NOT TV slack.
   Numerical mode covers searched points only. Certified mode covers the
   entire nuisance ball plus its analytic infinite exterior. */
struct jw_mc_result scalar jw_mcr_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar gamma, real scalar mcerr, real scalar reps,
    string scalar mode, real scalar nodes, real scalar search,
    real scalar outerreps, real scalar validreps, real scalar steps,
    real scalar certcells, real scalar radius, real scalar pmax,
    real scalar mcbeta, real scalar progress)
{
    struct jw_mc_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    real matrix biasSD,ES,E,GP,GV,train,boxes,distances
    real rowvector pgrid,maxtrain,maxvalid,maxdistance,eligible,counts,node
    real colvector ix
    real scalar alpha,alpha0,validationerror,tau,B,outerradius,basecv,zbeta,
        delta,n,o,j,i,c,ns,loc,selected,est,se,upper,minupper,k,worst,
        coordinate,sign,step,len,guard,mean,width
    alpha=1-level/100
    validationerror=mcerr/2
    alpha0=alpha-validationerror
    if (mode!="numerical" & mode!="certified") jw_core_error("mcmode() must be numerical or certified")
    if (hasmissing((level,gamma,mcerr,reps,nodes,search,outerreps,validreps,
        steps,certcells,radius,pmax,mcbeta)) | alpha<=0 | alpha>=1 |
        gamma<=0 | mcerr<=0 | gamma+mcerr>=alpha | mcbeta<=0 | mcbeta>=1 |
        reps<1 | reps!=floor(reps) | nodes<1 | nodes!=floor(nodes) |
        search<1 | search!=floor(search) | outerreps<1 | outerreps!=floor(outerreps) |
        validreps<1 | validreps!=floor(validreps) | steps<2 | steps!=floor(steps) |
        certcells<1 | certcells>65536 | certcells!=floor(certcells) | radius<0)
        jw_core_error("invalid refined McCloskey calibration controls")
    if (rows(observedA)!=P.nout | cols(observedA)!=rows(z) |
        rows(V)!=rows(z) | cols(V)!=rows(z) | P.p+P.k!=rows(z) |
        P.nout<1 | cols(P.outputs)!=P.nout | cols(P.masks)!=P.k |
        rows(P.actions)!=rows(P.masks)*P.nout | cols(P.actions)!=rows(z) |
        rows(P.sigmas)!=rows(P.masks) | cols(P.sigmas)!=P.nout |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0 |
        hasmissing(V) | hasmissing(z) | hasmissing(observedA))
        jw_core_error("incompatible refined McCloskey input dimensions")
    S=jw_mc_prepare(P,V)
    biasSD=jw_mcr_bias_sd(P,V)
    ns=rows(P.masks)
    zbeta=invnormal(1-mcbeta/2)
    tau=alpha0/2
    B=radius
    /* Fixed numerical domain; it does not depend on observed H or mcbeta. */
    if (B==0) B=sqrt(S.totalpenalty)+sqrt(invchi2(P.d,.99))
    outerradius=B+sqrt(invchi2(P.d,1-tau))
    basecv=invnormal(1-(alpha0-tau)/2)
    if (pmax==0) pmax=1-alpha0/100
    if (pmax<1-alpha0 | pmax>=1) jw_core_error("mcpmax() is outside the admissible probability range")
    pgrid=(1-alpha0):+((0..(steps-1))/(steps-1))*(pmax-(1-alpha0))
    R.domain_radius=B
    R.exterior_radius=outerradius
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,16,.)
    R.calibration=J(steps,2+3*P.nout,.)
    R.calibration[.,1]=(1::steps)
    R.calibration[.,2]=pgrid'
    for (o=1;o<=P.nout;o++) {
        est=observedA[o,.]*z
        se=sqrt(jw_fixed_variance(V,observedA[o,.]))
        R.table[o,1..2]=(est,se)
        R.diagnostics[o,.]=(P.outputs[o],P.d,ns,0,0,0,0,0,0,.,.,0,.,3,0,.)
    }
    eligible=J(1,P.nout,1)
    maxdistance=J(1,P.nout,0)
    maxtrain=J(1,P.nout*steps,.)
    ES=jw_infer_sqrt(V)
    if (mode=="numerical") {
        GP=jw_mc_nodes(P.d,nodes,B,max((S.boundary,1e-8))) \
            jw_mc_nodes(P.d,search,outerradius,max((S.boundary,1e-8)))
        GP=uniqrows(sort(GP,1..P.d))
        if (progress) printf("McCloskey: exact scalar-bias inner extrema; pilot at %g nodes\n",rows(GP))
        E=rnormal(outerreps,rows(z),0,1)*ES'
        train=J(rows(GP),P.nout*steps,.)
        for (i=1;i<=rows(GP);i++)
            train[i,.]=jw_mcr_risk(P,S,biasSD,E,GP[i,.],pgrid,B,basecv,zbeta)/outerreps
        /* Adapt only to the pilot bank. All resulting points are frozen
           before the independent validation bank is generated below. */
        for (k=1;k<=search;k++) {
            c=1+mod(k-1,P.nout*steps)
            ix=vec(selectindex(train[.,c]:==max(train[.,c])))
            worst=ix[1]
            coordinate=1+mod(floor((k-1)/2),P.d)
            sign=mod(k,2)==1 ? 1 : -1
            step=outerradius*2^(-1-floor((k-1)/(2*P.d)))
            node=GP[worst,.]
            node[coordinate]=node[coordinate]+sign*step
            len=sqrt(sum(node:^2))
            if (len>outerradius) node=node*(outerradius/len)
            GP=GP \ node
            train=train \ (jw_mcr_risk(P,S,biasSD,E,node,pgrid,B,basecv,zbeta)/outerreps)
        }
        maxtrain=colmax(train)
        GV=uniqrows(sort(GP,1..P.d))
        distances=J(rows(GV),P.nout,0)
    }
    else {
        GP=J(0,P.d,.)
        boxes=jw_wbic_cells(P.d,certcells)
        GV=outerradius*boxes[.,1..P.d]
        distances=jw_mcr_cell_distance(S,boxes,outerradius)
        maxdistance=colmax(distances)
        delta=validationerror/(rows(GV)*P.nout*steps)
        minupper=jw_rs_probability_upper(0,validreps,delta)
        for (o=1;o<=P.nout;o++) {
            eligible[o]=jw_mcr_np_upper(minupper,maxdistance[o])<=alpha0
            R.diagnostics[o,8]=rows(GV)
            R.diagnostics[o,12]=maxdistance[o]
            R.calibration[.,2+3*o]=J(steps,1,maxdistance[o])
            if (!eligible[o]) R.diagnostics[o,14]=4
        }
        if (sum(eligible)==0) return(R)
    }
    delta=validationerror/(rows(GV)*P.nout*steps)
    if (progress) printf("McCloskey: independent whole-rule validation, %g nodes and %g draws\n",rows(GV),validreps)
    E=rnormal(validreps,rows(z),0,1)*ES'
    maxvalid=J(1,P.nout*steps,0)
    for (i=1;i<=rows(GV);i++) {
        counts=jw_mcr_risk(P,S,biasSD,E,GV[i,.],pgrid,B,basecv,zbeta)
        for (o=1;o<=P.nout;o++) {
            for (j=1;j<=steps;j++) {
                c=(o-1)*steps+j
                upper=jw_rs_probability_upper(counts[c],validreps,delta)
                upper=jw_mcr_np_upper(upper,distances[i,o])
                maxvalid[c]=max((maxvalid[c],upper))
            }
        }
    }
    /* Only now use the observed precontrasts. */
    E=J(1,rows(z),0)
    D=jw_mc_draws(P,S,E,z[(P.p+1)..rows(z)]')
    guard=sum(D.standardized[1,.]:^2)>B^2
    for (o=1;o<=P.nout;o++) {
        loc=(o-1)*steps
        R.calibration[.,3*o]=maxtrain[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+1]=maxvalid[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+2]=J(steps,1,maxdistance[o])
        R.diagnostics[o,4..9]=(0,0,rows(GP),(mode=="numerical" ? outerreps : 0),rows(GV),validreps)
        R.diagnostics[o,11]=maxvalid[loc+steps]
        R.diagnostics[o,12]=maxdistance[o]
        R.diagnostics[o,16]=guard
        if (!eligible[o]) continue
        selected=0
        for (j=1;j<=steps;j++) {
            if (maxvalid[loc+j]<=alpha0) {
                selected=j
                break
            }
        }
        if (!selected) continue
        if (max(abs(P.actions[(o-1)*ns+D.chosen[1],.]-observedA[o,.]))>
            1e-8*(1+max(abs(observedA[o,.]))))
            jw_core_error("observed refined McCloskey action differs from WBIC selector")
        est=R.table[o,1]
        se=R.table[o,2]
        if (guard) width=abs(D.adjustment[1,o])+basecv*sqrt(S.basevar[o])
        else {
            mean=abs(D.adjustment[1,o])+zbeta*biasSD[D.chosen[1],o]
            width=jw_mcr_foldquantile(mean,D.sigma[1,o],pgrid[selected])
        }
        R.diagnostics[o,10]=pgrid[selected]
        R.diagnostics[o,11]=maxvalid[loc+selected]
        R.diagnostics[o,13]=width/se
        R.diagnostics[o,14]=0
        R.diagnostics[o,15]=mode=="certified"
        R.table[o,3..4]=(est-width,est+width)
    }
    return(R)
}

struct jw_mc_result scalar jw_infer_mc_refined(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps, real scalar ncells,
    real rowvector outputs, real scalar progress, string scalar mode,
    real scalar nodes, real scalar search, real scalar outerreps,
    real scalar validreps, real scalar steps, real scalar certcells,
    real scalar radius, real scalar pmax, real scalar mcbeta,
    string scalar fallback)
{
    struct jw_mc_result scalar R,Q
    struct jw_bb_result scalar BB
    struct jw_wbic_plan scalar P
    real rowvector active,failed
    real scalar i,o,nout,alpha0,status
    if (fallback!="bb" & fallback!="error") jw_core_error("mcfallback() must be bb or error")
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs)))
        jw_core_error("invalid McCloskey output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.calibration=J(0,0,.)
    R.domain_radius=0
    R.exterior_radius=0
    R.bb_diagnostics=J(0,0,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,14]=1
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,2..9]=(0,1,0,0,0,0,0,0)
            R.diagnostics[o,13..16]=(invnormal(1-(1-level/100)/2),2,1,0)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    Q=jw_mcr_from_plan(P,M.V,M.z,F.A[outputs,.],level,gamma,mcerr,reps,
        mode,nodes,search,outerreps,validreps,steps,certcells,radius,pmax,mcbeta,progress)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    R.calibration=Q.calibration
    R.domain_radius=Q.domain_radius
    R.exterior_radius=Q.exterior_radius
    failed=J(1,0,.)
    for (i=1;i<=cols(outputs);i++) {
        if (Q.diagnostics[i,14]>=3) failed=(failed,outputs[i])
    }
    if (cols(failed)>0) {
        if (fallback=="error") jw_core_error("McCloskey calibration did not validate; increase its budget or use mcfallback(bb)")
        alpha0=1-level/100-mcerr/2
        if (progress) printf("McCloskey: fresh certified BB fallback for %g outputs; not calibrated-M success\n",cols(failed))
        BB=jw_infer_bb(M,C,F,joint,scope,maxjoint,100*(1-alpha0),gamma,
            mcerr/2,reps,ncells,.05,failed,progress,"certified")
        R.bb_diagnostics=BB.diagnostics
        for (i=1;i<=cols(failed);i++) {
            o=failed[i]
            status=R.diagnostics[o,14]
            R.table[o,.]=BB.table[o,.]
            R.table[o,5]=.
            R.diagnostics[o,13]=BB.diagnostics[o,10]
            R.diagnostics[o,14]=status==4 ? 6 : 5
            R.diagnostics[o,15]=0
        }
    }
    return(R)
}

/* ===== jwselect_cis.mata ===== */
/* Conditional importance sampling for the actual WBIC selected law.
   H is whitened once. Conditional on H, the scalar baseline error is
   Gaussian; integrating it out leaves a smooth nuisance objective with
   the exact (possibly discontinuous) WBIC selector cached in the bank.

   Importance estimates and their ordinary standard errors are numerical
   search diagnostics, NOT probability certificates. Independent direct
   conditional samples are bounded in [0,1]; the KL bounds below are valid
   one-sided concentration bounds for those fractional observations.
   No exact-binomial calculation is applied to fractional tail counts. */

struct jw_cis_bank {
    real matrix X
    real colvector logq
    real matrix adjustment
    real matrix sigma
    real colvector chosen
    real matrix b
    real rowvector residualsd
    real matrix m0
}

struct jw_cis_bank scalar jw_cis_bank_fill(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real matrix X, real colvector logq)
{
    struct jw_cis_bank scalar K
    real matrix H,A,C
    real rowvector rv
    real scalar o,n,ns
    n=rows(X)
    ns=rows(P.masks)
    if (n<1 | cols(X)!=P.d | rows(logq)!=n | cols(logq)!=1 |
        hasmissing(X) | hasmissing(logq)) jw_core_error("invalid CIS bank dimensions")
    K.X=X
    K.logq=logq
    H=J(n,P.k,0)
    H[.,S.active]=X*S.Hroot'
    K.chosen=jw_mc_choose(P,S,H)
    K.adjustment=J(n,P.nout,.)
    K.sigma=J(n,P.nout,.)
    C=cholsolve(S.VH,S.cross)
    K.b=S.Hroot'*C
    rv=S.basevar-colsum(S.cross:*C)
    if (hasmissing(rv) | any(rv:<=1e-12*S.basevar))
        jw_core_error("CIS requires positive conditional baseline variance")
    K.residualsd=sqrt(rv)
    for (o=1;o<=P.nout;o++) {
        A=P.actions[(o-1)*ns:+K.chosen,(P.p+1)..(P.p+P.k)]
        K.adjustment[.,o]=rowsum(H:*A)
        K.sigma[.,o]=P.sigmas[K.chosen,o]
    }
    K.m0=K.adjustment+X*K.b
    return(K)
}

/* IID equally weighted Gaussian-mixture proposal. A random component is
   selected independently for each draw, so the ordinary variance formula
   describes this sampling scheme. log-sum-exp avoids mixture underflow. */
struct jw_cis_bank scalar jw_cis_bank_make(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real scalar n, real matrix centers,
    real scalar scale)
{
    real matrix X,D
    real colvector component,logq,ld,mx
    real scalar j,nc,constant
    if (hasmissing((n,scale)) | n<1 | n!=floor(n) | scale<=0 |
        rows(centers)<1 | cols(centers)!=P.d | hasmissing(centers))
        jw_core_error("invalid CIS proposal controls")
    nc=rows(centers)
    component=1:+floor(nc*runiform(n,1))
    X=rnormal(n,P.d,0,scale)+centers[component,.]
    constant=-.5*P.d*ln(2*pi())-P.d*ln(scale)
    D=(X-J(n,1,1)*centers[1,.])/scale
    logq=constant:-.5*rowsum(D:^2)
    for (j=2;j<=nc;j++) {
        D=(X-J(n,1,1)*centers[j,.])/scale
        ld=constant:-.5*rowsum(D:^2)
        mx=rowmax((logq,ld))
        logq=mx+ln(exp(logq-mx)+exp(ld-mx))
    }
    logq=logq:-ln(nc)
    return(jw_cis_bank_fill(P,S,X,logq))
}

/* Fresh direct samples at a frozen nuisance point; unlike pilot IS,
   their conditional rejection probabilities are bounded observations. */
struct jw_cis_bank scalar jw_cis_bank_at(struct jw_wbic_plan scalar P,
    struct jw_mc_state scalar S, real scalar n, real rowvector theta)
{
    real matrix Z
    real colvector logq
    if (missing(n) | n<1 | n!=floor(n) | cols(theta)!=P.d | hasmissing(theta))
        jw_core_error("invalid direct CIS simulation inputs")
    Z=rnormal(n,P.d,0,1)
    logq=-.5*P.d*ln(2*pi()):-.5*rowsum(Z:^2)
    return(jw_cis_bank_fill(P,S,Z+J(n,1,1)*theta,logq))
}

real colvector jw_cis_conditional(struct jw_cis_bank scalar K,
    real rowvector theta, real matrix widths, real scalar output)
{
    real colvector m,g
    real scalar s,n
    n=rows(K.X)
    if (cols(theta)!=cols(K.X) | hasmissing(theta) | rows(widths)!=n |
        cols(widths)!=cols(K.sigma) | output<1 | output>cols(K.sigma) |
        output!=floor(output) | hasmissing(widths[.,output]) |
        min(widths[.,output])<0) jw_core_error("invalid CIS conditional inputs")
    m=K.m0[.,output]:-theta*K.b[.,output]
    s=K.residualsd[output]
    g=normal((-widths[.,output]-m)/s)+normal((m-widths[.,output])/s)
    return(g)
}

/* Columns: risk, ordinary Monte Carlo SE, weight ESS, mean weight,
   maximum weight, then the gradient with respect to whitened nuisance.
   The IS mean is deliberately NOT self-normalized or clipped at one.
   Both changes would alter the objective and its analytic gradient. */
real rowvector jw_cis_evaluate(struct jw_cis_bank scalar K,
    real rowvector theta, real matrix widths, real scalar output)
{
    real matrix D
    real colvector g,lw,w,y,m,dg
    real rowvector gradient
    real scalar n,d,r,se,sw,ss,s
    n=rows(K.X)
    d=cols(K.X)
    g=jw_cis_conditional(K,theta,widths,output)
    D=K.X-J(n,1,1)*theta
    lw=(-.5*d*ln(2*pi()):-.5*rowsum(D:^2))-K.logq
    /* Squaring weights for ESS/variance must also remain representable. */
    if (hasmissing(lw) | max(lw)>300) jw_core_error("CIS importance weight overflow; rebuild proposal near search region")
    w=exp(lw)
    y=w:*g
    r=mean(y)
    se=.
    if (n>1) se=sqrt(sum((y:-r):^2)/(n*(n-1)))
    sw=sum(w)
    ss=sum(w:^2)
    m=K.m0[.,output]:-theta*K.b[.,output]
    s=K.residualsd[output]
    dg=(normalden((m-widths[.,output])/s)-
        normalden((-widths[.,output]-m)/s))/s
    gradient=(D'*(w:*g)/n-K.b[.,output]*mean(w:*dg))'
    return((r,se,(ss>0 ? sw^2/ss : 0),sw/n,max(w),gradient))
}

/* Bernoulli KL is the sharp Chernoff bound available from only the range
   [0,1] and sample mean; the Bernoulli mgf dominates every bounded law.
   With n IID bounded observations, P{true mean > upper} <= delta.
   Numerical bisection returns the outward endpoint. */
real scalar jw_cis_kl(real scalar x, real scalar p)
{
    if (x<=0) return(-ln(1-p))
    if (x>=1) return(-ln(p))
    return(x*ln(x/p)+(1-x)*ln((1-x)/(1-p)))
}

real scalar jw_cis_kl_upper(real scalar x, real scalar n, real scalar delta)
{
    real scalar lo,hi,mid,target,i
    if (hasmissing((x,n,delta)) | x<0 | x>1 | n<1 | n!=floor(n) |
        delta<=0 | delta>=1) jw_core_error("invalid bounded CIS confidence inputs")
    if (x>=1) return(1)
    if (x==0) return(min((1,1-exp(ln(delta)/n)+1e-14)))
    lo=x
    hi=1
    target=-ln(delta)/n
    for (i=1;i<=70;i++) {
        mid=(lo+hi)/2
        if (mid>=1) break
        if (jw_cis_kl(x,mid)>target) hi=mid
        else lo=mid
        if (hi-lo<1e-14) break
    }
    return(min((1,hi+1e-14)))
}

real scalar jw_cis_kl_lower(real scalar x, real scalar n, real scalar delta)
{
    return(max((0,1-jw_cis_kl_upper(1-x,n,delta)-1e-14)))
}

real rowvector jw_cis_project(real rowvector theta, real rowvector center,
    real scalar radius)
{
    real scalar length
    length=sqrt(sum((theta-center):^2))
    if (length>radius) return(center+(theta-center)*(radius/length))
    return(theta)
}

/* Projected gradient multistart on the frozen, smooth pilot objective.
   maxeval is a TOTAL bound across starts. Starts and all accepted points
   are in the requested ball. Trial steps with weight ESS below max(50,
   .01*n) or mean weight outside [.5,1.5] are rejected as poorly supported.
   This trust rule is a search diagnostic, NOT a probability guarantee.
   If no start passes it, the best unsupported start is returned with its
   diagnostic values. Callers must recenter proposals and validate anew.
   This is a local numerical search, never a proof of a global extremum.
   Return: theta[1..d], risk, SE, ESS, meanweight, maxweight, evaluations. */
real rowvector jw_cis_maximize(struct jw_cis_bank scalar K,
    real matrix widths, real scalar output, real rowvector center,
    real scalar radius, real matrix starts, real scalar maxeval)
{
    real rowvector theta,trial,value,candidate,best,besttheta,gradient,stepvec
    real scalar d,nstart,i,k,used,remaining,budget,step,gn,move,accepted,localused,
        trusted,besttrusted,miness
    d=cols(K.X)
    nstart=rows(starts)
    if (cols(center)!=d | hasmissing(center) | missing(radius) | radius<0 |
        nstart<1 | cols(starts)!=d | hasmissing(starts) |
        missing(maxeval) | maxeval<1 | maxeval!=floor(maxeval))
        jw_core_error("invalid CIS optimization controls")
    used=0
    best=J(1,5+d,.)
    best[1]=-1
    besttheta=center
    besttrusted=0
    miness=max((50,.01*rows(K.X)))
    for (i=1;i<=nstart & used<maxeval;i++) {
        remaining=maxeval-used
        budget=max((1,floor(remaining/(nstart-i+1))))
        theta=jw_cis_project(starts[i,.],center,radius)
        value=jw_cis_evaluate(K,theta,widths,output)
        used=used+1
        localused=1
        trusted=(value[3]>=miness & value[4]>=.5 & value[4]<=1.5)
        if ((trusted & !besttrusted) | (trusted==besttrusted & value[1]>best[1])) {
            best=value
            besttheta=theta
            besttrusted=trusted
        }
        if (!trusted) continue
        step=max((.05,min((1,radius/3))))
        while (localused<budget & used<maxeval & radius>0) {
            gradient=value[6..(5+d)]
            gn=sqrt(sum(gradient:^2))
            if (missing(gn) | gn<=1e-10) break
            accepted=0
            for (k=1;k<=12 & localused<budget & used<maxeval;k++) {
                trial=jw_cis_project(theta+(step/gn)*gradient,center,radius)
                stepvec=trial-theta
                move=sqrt(sum(stepvec:^2))
                if (move<=1e-7*(1+radius)) break
                candidate=jw_cis_evaluate(K,trial,widths,output)
                used=used+1
                localused=localused+1
                trusted=(candidate[3]>=miness & candidate[4]>=.5 & candidate[4]<=1.5)
                if (trusted & candidate[1]>best[1]) {
                    best=candidate
                    besttheta=trial
                }
                if (trusted & candidate[1]>=value[1]+1e-4*sum(gradient:*stepvec)) {
                    theta=trial
                    value=candidate
                    accepted=1
                    step=min((max((radius,1)),step*1.5))
                    break
                }
                step=step/2
            }
            if (!accepted) break
        }
    }
    return((besttheta,best[1..5],used))
}

/* ===== jwselect_bb_cis.mata ===== */
/* Conditional importance-sampling sharp Berger--Boos engine (0.8.4).
   The pilot integrates the baseline Gaussian error analytically, caches
   selection on fixed H banks, and optimizes smooth nuisance risks.
   Fresh direct conditional banks validate frozen nuisance witnesses.
   A numerical maximum is never a global certificate. Certified output
   uses either an independent full-metric NP cover or the analytic cap. */

real scalar jw_bbc_quantile(struct jw_cis_bank scalar K,
    real rowvector theta, real scalar output, real scalar tail,
    real scalar cap, real scalar weighted)
{
    real scalar lo,hi,mid,i,r
    real rowvector ev
    lo=0
    hi=cap
    for (i=1;i<=40;i++) {
        mid=(lo+hi)/2
        if (weighted) {
            ev=jw_cis_evaluate(K,theta,mid*K.sigma,output)
            r=ev[1]
        }
        else r=mean(jw_cis_conditional(K,theta,mid*K.sigma,output))
        if (r>tail) lo=mid
        else hi=mid
        if (hi-lo<1e-6*(1+hi)) break
    }
    return((lo+hi)/2)
}

/* All starts are feasible, and include the observed center and the
   projection of zero. A score-boundary scale supplements random starts. */
real matrix jw_bbc_starts(struct jw_mc_state scalar S,
    real rowvector center, real scalar radius, real scalar n)
{
    real matrix X
    real rowvector v
    real scalar i,j,len,r
    X=J(n,1,1)*center
    if (n>1) X[2,.]=jw_cis_project(J(1,cols(center),0),center,radius)
    for (i=3;i<=n;i++) {
        v=rnormal(1,cols(center),0,1)
        len=sqrt(sum(v:^2))
        if (len==0) {
            v[1]=1
            len=1
        }
        r=mod(i,3)==0 ? radius : min((radius,max((S.boundary,.5*radius))))
        X[i,.]=center+r*v/len
    }
    return(X)
}

struct jw_bb_result scalar jw_bb_cis_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar gamma, real scalar mcerr, real scalar reps,
    string scalar mode, real scalar cells, real scalar tol, real scalar progress,
    real scalar isreps, real scalar isstarts, real scalar isiter)
{
    struct jw_bb_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    struct jw_cis_bank scalar K,KV
    real matrix starts,witnesses,bankcenters,metric,boxes,coverpoints,coverdist,
        pointsets,allpoints,gridall,covergrid,coverupper,directions
    real rowvector center,theta,fit,cap,low,pointupper,numerical,upper,
        grid,pcap,rad,box,u,ev,besttheta,ng,good,covermax,eligible,direction
    real colvector vals,ix
    real scalar alpha,tail,o,i,j,k,stage,ns,d,cv,lo,hi,target,
        pilotdraws,validationdraws,evals,miness,nw,delta,q,pointlo,
        pointup,prob,pl,pu,est,se,critical,status,capused,nsplit,
        lowerbound,gridstep,maxdist,step,number,offset,ncover,coverdraws,
        gridn,observed,iteration,maxq
    alpha=1-level/100
    tail=alpha-gamma-mcerr
    if (mode!="numerical" & mode!="certified") jw_core_error("bbmode() must be numerical or certified")
    if (hasmissing((level,gamma,mcerr,reps,cells,tol,isreps,isstarts,isiter)) |
        alpha<=0 | alpha>=1 | gamma<=0 | mcerr<=0 | tail<=0 |
        reps<1 | reps!=floor(reps) | cells<1 | cells>65536 |
        cells!=floor(cells) | tol<=0 | isreps<64 | isreps!=floor(isreps) |
        isstarts<1 | isstarts!=floor(isstarts) | isiter<isstarts |
        isiter!=floor(isiter)) jw_core_error("invalid conditional BB controls")
    if (P.d<1 | rows(observedA)!=P.nout | cols(observedA)!=rows(z) |
        rows(V)!=rows(z) | cols(V)!=rows(z) | P.p+P.k!=rows(z) |
        rows(P.actions)!=rows(P.masks)*P.nout | cols(P.actions)!=rows(z) |
        rows(P.sigmas)!=rows(P.masks) | cols(P.sigmas)!=P.nout |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0 |
        hasmissing(V) | hasmissing(z) | hasmissing(observedA))
        jw_core_error("incompatible conditional BB input dimensions")
    if (missing(P.radius) | abs(P.radius-sqrt(invchi2(P.d,1-gamma)))>
        1e-10*(1+P.radius)) jw_core_error("BB plan radius disagrees with gamma()")
    S=jw_mc_prepare(P,V)
    d=P.d
    ns=rows(P.masks)
    center=P.h[S.active]'*S.Hinv'
    box=(J(1,d,0),J(1,d,1))
    cap=jw_bb_normal_cap(P,box,tail)
    metric=jw_bb_metric(P,S)
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,16,.)
    R.cis_diagnostics=J(P.nout,12,.)
    R.cis_points=J(P.nout,d,.)
    low=J(1,P.nout,0)
    pointupper=cap
    numerical=J(1,P.nout,.)
    upper=cap
    allpoints=J(0,d,.)
    ng=J(1,P.nout,0)
    /* Each output gets a proposal adapted to its own least-favorable
       candidates. This avoids diluting an aggregate proposal across
       unrelated cohort targets. */
    for (o=1;o<=P.nout;o++) {
        starts=jw_bbc_starts(S,center,P.radius,isstarts)
        /* Bias directions are often much more useful than random axes in
           aggregate problems. Retain the center/zero starts and some
           random starts, but seed both signs of the strongest loading. */
        directions=P.actions[((o-1)*ns+1)..(o*ns),P.p:+S.active]*S.Hroot
        vals=rowsum(directions:^2)
        ix=vec(selectindex(vals:==max(vals)))
        if (max(vals)>0) {
            direction=directions[ix[1],.]/sqrt(max(vals))
            for (i=3;i<=min((6,isstarts));i++) {
                step=i<=4 ? min((P.radius,S.boundary)) : P.radius
                starts[i,.]=center+(mod(i,2)==1 ? 1 : -1)*step*direction
            }
        }
        bankcenters=starts
        witnesses=starts[1,.]
        cv=invnormal(1-tail/2)
        pilotdraws=0
        evals=0
        miness=isreps
        theta=center
        for (stage=1;stage<=2;stage++) {
            K=jw_cis_bank_make(P,S,isreps,bankcenters,1)
            pilotdraws=pilotdraws+isreps
            lo=0
            hi=cap[o]
            /* Repeated risk maximization and scalar critical bisection.
               All decisions here use the disposable pilot bank. */
            for (iteration=1;iteration<=8;iteration++) {
                fit=jw_cis_maximize(K,cv*K.sigma,o,center,P.radius,starts,isiter)
                theta=fit[1..d]
                evals=evals+fit[d+6]
                miness=min((miness,fit[d+3]))
                witnesses=witnesses\theta
                if (fit[d+1]>tail) lo=cv
                else hi=cv
                starts[1,.]=theta
                cv=(lo+hi)/2
            }
            /* Recenter the next bank on the current maximizer while
               retaining diverse proposal support. */
            bankcenters=theta\starts
            cv=jw_bbc_quantile(K,theta,o,tail,cap[o],1)
        }
        witnesses=uniqrows(sort(witnesses,1..d))
        /* Freeze all witness locations and all confidence-grid thresholds
           before generating validation data. Critical roots themselves
           are descriptive estimates, never confidence certificates. */
        gridstep=max((.025,min((.05,tol/2))))
        /* The coarse whole-range grid still yields informative witness
           brackets when a noisy high-dimensional pilot mislocates the
           root. The additional local grid resolves well-behaved pilots.
           Both grids are frozen before any validation observation. */
        grid=J(1,86,0)
        grid[1..65]=(0..64)*(cap[o]/64)
        for (j=1;j<=21;j++) grid[j+65]=max((0,min((cap[o],cv+(j-11)*gridstep))))
        grid=uniqrows(sort(grid',1))'
        delta=mcerr/(4*rows(witnesses)*P.nout*cols(grid))
        maxq=-1
        pointup=0
        for (i=1;i<=rows(witnesses);i++) {
            theta=witnesses[i,.]
            KV=jw_cis_bank_at(P,S,reps,theta)
            q=jw_bbc_quantile(KV,theta,o,tail,cap[o],0)
            if (q>maxq) {
                maxq=q
                besttheta=theta
            }
            pl=0
            pu=cap[o]
            for (j=1;j<=cols(grid);j++) {
                prob=mean(jw_cis_conditional(KV,theta,grid[j]*KV.sigma,o))
                if (jw_cis_kl_lower(prob,reps,delta)>tail) pl=max((pl,grid[j]))
                if (jw_cis_kl_upper(prob,reps,delta)<=tail) pu=min((pu,grid[j]))
            }
            low[o]=max((low[o],pl))
            pointup=max((pointup,pu))
        }
        numerical[o]=maxq
        pointupper[o]=pointup
        R.cis_points[o,.]=besttheta
        ng[o]=rows(witnesses)
        allpoints=allpoints\witnesses
        validationdraws=rows(witnesses)*reps
        R.cis_diagnostics[o,.]=(P.outputs[o],pilotdraws,validationdraws,
            isstarts,evals,miness,low[o],pointupper[o],numerical[o],cap[o],0,0)
        if (progress) printf("Conditional BB output %g: critical %8.4f; %g frozen witnesses; pilot ESS >= %g\n",
            P.outputs[o],numerical[o],rows(witnesses),miness)
    }
    /* Optional genuine global cover. Its nodes, thresholds and allocation
       are frozen before its fresh conditional draws. The metric is that
       of the FULL (baseline,H) law, not the marginal nuisance law. */
    boxes=jw_wbic_cells(d,cells)
    ncover=rows(boxes)
    coverpoints=J(ncover,d,.)
    coverdist=J(ncover,P.nout,.)
    for (i=1;i<=ncover;i++) {
        theta=jw_bb_feasible(boxes[i,.])
        coverpoints[i,.]=center+(P.L[S.active,.]*theta')'*S.Hinv'
        coverdist[i,.]=jw_bb_radius(metric,boxes[i,.],theta)
    }
    covergrid=J(3,P.nout,.)
    for (o=1;o<=P.nout;o++) {
        step=max((.075,tol/2))
        covergrid[.,o]=min((cap[o],max((pointupper[o],numerical[o]+step))))\
            min((cap[o],max((pointupper[o],numerical[o]+2*step))))\
            min((cap[o],max((pointupper[o],numerical[o]+4*step))))
    }
    delta=mcerr/(4*ncover*P.nout*3)
    eligible=J(1,P.nout,0)
    for (o=1;o<=P.nout;o++) {
        maxdist=max(coverdist[.,o])
        eligible[o]=(mode=="certified" & !missing(maxdist) &
            jw_mcr_np_upper(jw_cis_kl_upper(0,reps,delta),maxdist)<=tail)
    }
    coverupper=J(3,P.nout,0)
    coverdraws=0
    if (sum(eligible)>0) {
        for (i=1;i<=ncover;i++) {
            theta=coverpoints[i,.]
            KV=jw_cis_bank_at(P,S,reps,theta)
            coverdraws=coverdraws+reps
            for (o=1;o<=P.nout;o++) if (eligible[o]) {
                for (j=1;j<=3;j++) {
                    prob=mean(jw_cis_conditional(KV,theta,covergrid[j,o]*KV.sigma,o))
                    prob=jw_cis_kl_upper(prob,reps,delta)
                    prob=jw_mcr_np_upper(prob,coverdist[i,o])
                    coverupper[j,o]=max((coverupper[j,o],prob))
                }
            }
        }
        for (o=1;o<=P.nout;o++) if (eligible[o]) {
            for (j=1;j<=3;j++) if (coverupper[j,o]<=tail)
                upper[o]=min((upper[o],covergrid[j,o]))
        }
    }
    D=jw_mc_draws(P,S,J(1,rows(z),0),z[(P.p+1)..rows(z)]')
    for (o=1;o<=P.nout;o++) {
        if (max(abs(P.actions[(o-1)*ns+D.chosen[1],.]-observedA[o,.]))>
            1e-8*(1+max(abs(observedA[o,.])))) jw_core_error("observed conditional BB action differs from actual WBIC selector")
        est=observedA[o,.]*z
        se=sqrt(jw_fixed_variance(V,observedA[o,.]))
        critical=mode=="certified" ? upper[o] : numerical[o]
        status=mode=="numerical" ? 1 : (upper[o]-low[o]<=tol ? 0 : 4)
        capused=(mode=="certified" & upper[o]>=cap[o]-1e-12)
        maxdist=max(coverdist[.,o])
        R.table[o,.]=(est,se,est-critical*se,est+critical*se,.)
        R.diagnostics[o,.]=(P.outputs[o],d,ns,ng[o]+(eligible[o] ? ncover : 0),
            R.cis_diagnostics[o,2]+R.cis_diagnostics[o,3]+(eligible[o] ? coverdraws : 0),
            tail,low[o],upper[o],upper[o]-low[o],critical,status,
            (mode=="certified"),(mode=="certified" & upper[o]-low[o]<=tol),
            capused,maxdist,numerical[o])
        R.cis_diagnostics[o,3]=R.cis_diagnostics[o,3]+(eligible[o] ? coverdraws : 0)
        R.cis_diagnostics[o,10..12]=(upper[o],mode=="certified",
            mode=="certified" & upper[o]-low[o]<=tol)
    }
    return(R)
}

struct jw_bb_result scalar jw_infer_bb_cis(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps, real scalar cells,
    real scalar tol, real rowvector outputs, real scalar progress, string scalar mode,
    real scalar isreps, real scalar isstarts, real scalar isiter)
{
    struct jw_bb_result scalar R,Q
    struct jw_wbic_plan scalar P
    real rowvector active
    real scalar o,i,nout,cv
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout |
        any(outputs:!=floor(outputs))) jw_core_error("invalid conditional BB output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.cis_diagnostics=J(nout,12,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,11]=9
        R.cis_diagnostics[o,1]=o
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        R.cis_points=J(nout,0,.)
        cv=invnormal(1-(1-level/100)/2)
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,.]=(o,0,1,0,0,1-level/100,cv,cv,0,cv,2,1,1,0,0,cv)
            R.cis_diagnostics[o,.]=(o,0,0,0,0,.,cv,cv,cv,cv,1,1)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    Q=jw_bb_cis_from_plan(P,M.V,M.z,F.A[outputs,.],level,gamma,mcerr,reps,
        mode,cells,tol,progress,isreps,isstarts,isiter)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    R.cis_diagnostics[outputs,.]=Q.cis_diagnostics
    R.cis_points=J(nout,P.d,.)
    R.cis_points[outputs,.]=Q.cis_points
    return(R)
}

/* ===== jwselect_mc_cis.mata ===== */
/* Conditional importance-sampling calibration of the McCloskey-inspired
   WBIC rule. The rule itself is unchanged from jwselect_mc_refined.mata.
   Integration conditions on the complete vector of selection contrasts.
   A fixed H bank caches actual WBIC selections, and the independent
   conditional validation observations lie in [0,1]. Their confidence
   bounds use bounded-observation Chernoff/KL inequalities, not binomial
   confidence intervals for fractional counts.

   Numerical mode validates frozen searched points only. Certified mode
   covers the nuisance ball with the full-joint Gaussian NP bound and the
   existing analytic exterior argument. An optimized point is NEVER a
   certificate for the continuum. The reported calibration is the finite
   Gaussian WBIC adaptation, not the literal published McCloskey algorithm.
*/

/* Vectorized deterministic upper quantile bracket for every cached H. */
real matrix jw_mcc_widths(struct jw_cis_bank scalar K,
    struct jw_mc_state scalar S, real matrix biasSD, real scalar p,
    real scalar B, real scalar basecv, real scalar zbeta)
{
    real matrix m,lo,hi,mid,cdf,W
    real colvector guard,ix
    real scalar o,it
    m=(abs(K.adjustment):+zbeta*biasSD[K.chosen,.]):/K.sigma
    lo=J(rows(m),cols(m),0)
    hi=m:+invnormal((1+p)/2)
    for (it=1;it<=44;it++) {
        mid=(lo+hi)/2
        cdf=normal(mid-m)-normal(-mid-m)
        lo=(cdf:<p):*mid+(cdf:>=p):*lo
        hi=(cdf:>=p):*mid+(cdf:<p):*hi
    }
    W=hi:*K.sigma
    guard=rowsum(K.X:^2):>B^2
    ix=vec(selectindex(guard))
    if (rows(ix)>0) {
        for (o=1;o<=cols(W);o++)
            W[ix,o]=abs(K.adjustment[ix,o]):+basecv*sqrt(S.basevar[o])
    }
    return(W)
}

struct jw_mc_result scalar jw_mcc_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar gamma, real scalar mcerr, real scalar reps,
    string scalar mode, real scalar nodes, real scalar search,
    real scalar outerreps, real scalar validreps, real scalar steps,
    real scalar certcells, real scalar radius, real scalar pmax,
    real scalar mcbeta, real scalar progress, real scalar isreps,
    real scalar isstarts, real scalar isiter)
{
    struct jw_mc_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    struct jw_cis_bank scalar K,KL,KV
    real matrix biasSD,GP,GV,boxes,distances,W,WL,starts,values,pointindex
    real rowvector pgrid,maxtrain,maxvalid,pointvalid,minvalid,maxdistance,eligible,
        answer,node,v,evals,miness,pilotdraws,validationdraws
    real colvector ix,orderix,g
    real scalar alpha,alpha0,validationerror,tau,B,outerradius,basecv,zbeta,
        delta,o,j,i,c,ns,loc,selected,est,se,upper,lower,minupper,d,k,
        guard,mean,width,nstart,trust,npoint,dirradius,q
    alpha=1-level/100
    validationerror=mcerr/2
    alpha0=alpha-validationerror
    if (mode!="numerical" & mode!="certified") jw_core_error("mcmode() must be numerical or certified")
    if (hasmissing((level,gamma,mcerr,reps,nodes,search,outerreps,validreps,
        steps,certcells,radius,pmax,mcbeta,isreps,isstarts,isiter)) |
        alpha<=0 | alpha>=1 | gamma<=0 | mcerr<=0 | gamma+mcerr>=alpha |
        mcbeta<=0 | mcbeta>=1 | reps<1 | reps!=floor(reps) |
        nodes<1 | nodes!=floor(nodes) | search<1 | search!=floor(search) |
        outerreps<1 | outerreps!=floor(outerreps) | validreps<1 |
        validreps!=floor(validreps) | steps<2 | steps!=floor(steps) |
        certcells<1 | certcells>65536 | certcells!=floor(certcells) | radius<0 |
        isreps<100 | isreps!=floor(isreps) | isstarts<1 |
        isstarts!=floor(isstarts) | isiter<1 | isiter!=floor(isiter))
        jw_core_error("invalid conditional McCloskey controls")
    if (rows(observedA)!=P.nout | cols(observedA)!=rows(z) |
        rows(V)!=rows(z) | cols(V)!=rows(z) | P.p+P.k!=rows(z) |
        P.nout<1 | cols(P.outputs)!=P.nout | cols(P.masks)!=P.k |
        rows(P.actions)!=rows(P.masks)*P.nout | cols(P.actions)!=rows(z) |
        rows(P.sigmas)!=rows(P.masks) | cols(P.sigmas)!=P.nout |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0 |
        hasmissing(V) | hasmissing(z) | hasmissing(observedA))
        jw_core_error("incompatible conditional McCloskey dimensions")
    S=jw_mc_prepare(P,V)
    biasSD=jw_mcr_bias_sd(P,V)
    d=P.d
    ns=rows(P.masks)
    zbeta=invnormal(1-mcbeta/2)
    tau=alpha0/2
    B=radius
    if (B==0) B=sqrt(S.totalpenalty)+sqrt(invchi2(d,.99))
    outerradius=B+sqrt(invchi2(d,1-tau))
    basecv=invnormal(1-(alpha0-tau)/2)
    if (pmax==0) pmax=1-alpha0/100
    if (pmax<1-alpha0 | pmax>=1) jw_core_error("mcpmax() is outside the admissible probability range")
    pgrid=(1-alpha0):+((0..(steps-1))/(steps-1))*(pmax-(1-alpha0))
    R.domain_radius=B
    R.exterior_radius=outerradius
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,16,.)
    R.calibration=J(steps,2+3*P.nout,.)
    R.bb_diagnostics=J(0,0,.)
    R.cis_diagnostics=J(P.nout,12,.)
    R.cis_points=J(P.nout,d,.)
    R.calibration[.,1]=(1::steps)
    R.calibration[.,2]=pgrid'
    for (o=1;o<=P.nout;o++) {
        R.table[o,1..2]=(observedA[o,.]*z,sqrt(jw_fixed_variance(V,observedA[o,.])))
        R.diagnostics[o,.]=(P.outputs[o],d,ns,0,0,0,0,0,0,.,.,0,.,3,0,.)
        R.cis_diagnostics[o,.]=(P.outputs[o],0,0,0,0,.,.,.,.,1,0,0)
    }
    eligible=J(1,P.nout,1)
    maxdistance=J(1,P.nout,0)
    maxtrain=J(1,P.nout*steps,.)
    minvalid=J(1,P.nout*steps,0)
    maxvalid=J(1,P.nout*steps,0)
    pointvalid=J(1,P.nout*steps,0)
    pointindex=J(P.nout,steps,1)
    evals=J(1,P.nout,0)
    miness=J(1,P.nout,.)
    pilotdraws=J(1,P.nout,0)
    validationdraws=J(1,P.nout,0)
    GP=J(0,d,.)
    if (mode=="certified") {
        boxes=jw_wbic_cells(d,certcells)
        GV=outerradius*boxes[.,1..d]
        distances=jw_mcr_cell_distance(S,boxes,outerradius)
        maxdistance=colmax(distances)
        delta=validationerror/(2*rows(GV)*P.nout*steps)
        minupper=jw_cis_kl_upper(0,validreps,delta)
        for (o=1;o<=P.nout;o++) {
            eligible[o]=jw_mcr_np_upper(minupper,maxdistance[o])<=alpha0
            R.diagnostics[o,8]=rows(GV)
            R.diagnostics[o,12]=maxdistance[o]
            R.calibration[.,2+3*o]=J(steps,1,maxdistance[o])
            if (!eligible[o]) R.diagnostics[o,14]=4
        }
        /* This is a proof of insufficiency of THIS certificate budget,
           not evidence that the true rule has excessive rejection. */
        if (sum(eligible)==0) return(R)
    }
    else {
        nstart=max((nodes,2*d+1,isstarts))
        GP=jw_mc_nodes(d,nstart,min((outerradius,1.5*sqrt(d)*max((S.boundary,.1)))),max((S.boundary,1e-8)))
        /* Coordinate starts alone miss diagonal aggregate bias. Include
           both signs of each target's largest candidate-bias direction.
           These depend only on candidate loadings and covariance. */
        dirradius=min((outerradius,sqrt(d)*max((S.boundary,.1))))
        for (o=1;o<=P.nout;o++) {
            orderix=order(biasSD[.,o],-1)
            q=orderix[1]
            v=P.actions[(o-1)*ns+q,P.p:+S.active]*S.Hroot
            if (sqrt(sum(v:^2))>1e-12) {
                v=dirradius*v/sqrt(sum(v:^2))
                GP=GP\v\(-v)\(.5*v)\(-.5*v)
            }
        }
        GP=uniqrows(sort(GP,1..d))
        K=jw_cis_bank_make(P,S,isreps,GP,1)
        pilotdraws=J(1,P.nout,isreps)
        GV=J(0,d,.)
        values=J(rows(GP),1,.)
        if (progress) printf("McCloskey CIS: cached selector bank %g draws; %g starting locations\n",isreps,rows(GP))
        for (j=1;j<=steps;j++) {
            W=jw_mcc_widths(K,S,biasSD,pgrid[j],B,basecv,zbeta)
            for (o=1;o<=P.nout;o++) {
                for (i=1;i<=rows(GP);i++) {
                    v=jw_cis_evaluate(K,GP[i,.],W,o)
                    values[i]=v[1]
                }
                orderix=order(values,-1)
                starts=GP[orderix[1..min((rows(GP),isstarts))],.]
                answer=jw_cis_maximize(K,W,o,J(1,d,0),outerradius,starts,isiter)
                node=answer[1..d]
                c=(o-1)*steps+j
                maxtrain[c]=answer[d+1]
                miness[o]=min((miness[o],answer[d+3]))
                evals[o]=evals[o]+answer[d+6]+rows(GP)
                GV=GV\node\starts[1,.]
                /* A new local bank makes the smooth objective reliable
                   near the representative extremum. Its one-unit trust
                   ball is contained in the global nuisance ball. */
                {
                    KL=jw_cis_bank_at(P,S,isreps,node)
                    WL=jw_mcc_widths(KL,S,biasSD,pgrid[j],B,basecv,zbeta)
                    trust=min((1,max((0,outerradius-sqrt(sum(node:^2))))))
                    if (trust>1e-8) {
                        answer=jw_cis_maximize(KL,WL,o,node,trust,node,isiter)
                        GV=GV\answer[1..d]
                        maxtrain[c]=max((maxtrain[c],answer[d+1]))
                        miness[o]=min((miness[o],answer[d+3]))
                        evals[o]=evals[o]+answer[d+6]
                    }
                    pilotdraws[o]=pilotdraws[o]+isreps
                }
            }
        }
        GV=uniqrows(sort(GV,1..d))
        distances=J(rows(GV),P.nout,0)
    }
    /* GV and the p-grid are now frozen. No validation draw is reused for
       nuisance optimization, proposal construction, or adding p levels. */
    npoint=rows(GV)
    delta=validationerror/(2*npoint*P.nout*steps)
    if (progress) printf("McCloskey CIS: fresh bounded conditional validation at %g frozen locations\n",npoint)
    for (i=1;i<=npoint;i++) {
        KV=jw_cis_bank_at(P,S,validreps,GV[i,.])
        for (j=1;j<=steps;j++) {
            W=jw_mcc_widths(KV,S,biasSD,pgrid[j],B,basecv,zbeta)
            for (o=1;o<=P.nout;o++) {
                if (!eligible[o]) continue
                c=(o-1)*steps+j
                g=jw_cis_conditional(KV,GV[i,.],W,o)
                mean=mean(g)
                upper=jw_cis_kl_upper(mean,validreps,delta)
                lower=jw_cis_kl_lower(mean,validreps,delta)
                if (upper>pointvalid[c]) {
                    pointvalid[c]=upper
                    pointindex[o,j]=i
                }
                upper=jw_mcr_np_upper(upper,distances[i,o])
                maxvalid[c]=max((maxvalid[c],upper))
                minvalid[c]=max((minvalid[c],lower))
            }
        }
    }
    validationdraws=eligible:*(npoint*validreps)
    /* Only now use the observed precontrasts to evaluate the calibrated
       rule. Holding covariance fixed, changing observed H cannot affect
       the selected p level or any calibration/search decision. */
    D=jw_mc_draws(P,S,J(1,rows(z),0),z[(P.p+1)..rows(z)]')
    guard=sum(D.standardized[1,.]:^2)>B^2
    for (o=1;o<=P.nout;o++) {
        loc=(o-1)*steps
        R.calibration[.,3*o]=maxtrain[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+1]=maxvalid[(loc+1)..(loc+steps)]'
        R.calibration[.,3*o+2]=J(steps,1,maxdistance[o])
        R.diagnostics[o,4..9]=(0,0,rows(GP),(mode=="numerical" ? isreps : 0),npoint,(eligible[o] ? validreps : 0))
        R.diagnostics[o,11]=maxvalid[loc+steps]
        R.diagnostics[o,12]=maxdistance[o]
        R.diagnostics[o,16]=guard
        R.cis_diagnostics[o,2..6]=(pilotdraws[o],validationdraws[o],(mode=="numerical" ? isstarts : 0),evals[o],miness[o])
        if (!eligible[o]) continue
        selected=0
        for (j=1;j<=steps;j++) if (maxvalid[loc+j]<=alpha0) {
            selected=j
            break
        }
        if (!selected) {
            R.cis_diagnostics[o,7..8]=(minvalid[loc+steps],pointvalid[loc+steps])
            continue
        }
        if (max(abs(P.actions[(o-1)*ns+D.chosen[1],.]-observedA[o,.]))>
            1e-8*(1+max(abs(observedA[o,.]))))
            jw_core_error("observed conditional McCloskey action differs from WBIC selector")
        est=R.table[o,1]
        se=R.table[o,2]
        if (guard) width=abs(D.adjustment[1,o])+basecv*sqrt(S.basevar[o])
        else {
            mean=abs(D.adjustment[1,o])+zbeta*biasSD[D.chosen[1],o]
            width=jw_mcr_foldquantile(mean,D.sigma[1,o],pgrid[selected])
        }
        R.diagnostics[o,10]=pgrid[selected]
        R.diagnostics[o,11]=maxvalid[loc+selected]
        R.diagnostics[o,13]=width/se
        R.diagnostics[o,14]=0
        R.diagnostics[o,15]=mode=="certified"
        R.table[o,3..4]=(est-width,est+width)
        R.cis_diagnostics[o,7..12]=(minvalid[loc+selected],pointvalid[loc+selected],
            pgrid[selected],(mode=="certified" ? max((maxvalid[loc+selected],alpha0)) : 1),mode=="certified",0)
        R.cis_points[o,.]=GV[pointindex[o,selected],.]
    }
    return(R)
}

struct jw_mc_result scalar jw_infer_mc_cis(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real scalar reps, real scalar ncells,
    real rowvector outputs, real scalar progress, string scalar mode,
    real scalar nodes, real scalar search, real scalar outerreps,
    real scalar validreps, real scalar steps, real scalar certcells,
    real scalar radius, real scalar pmax, real scalar mcbeta,
    string scalar fallback, real scalar isreps, real scalar isstarts,
    real scalar isiter)
{
    struct jw_mc_result scalar R,Q
    struct jw_bb_result scalar BB
    struct jw_wbic_plan scalar P
    real rowvector active,failed
    real scalar i,o,nout,alpha0,status
    if (fallback!="bb" & fallback!="error") jw_core_error("mcfallback() must be bb or error")
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs)))
        jw_core_error("invalid conditional McCloskey output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.calibration=J(0,0,.)
    R.domain_radius=0
    R.exterior_radius=0
    R.bb_diagnostics=J(0,0,.)
    R.cis_diagnostics=J(nout,12,.)
    R.cis_points=J(nout,0,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,14]=1
        R.cis_diagnostics[o,1]=o
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,2..9]=(0,1,0,0,0,0,0,0)
            R.diagnostics[o,13..16]=(invnormal(1-(1-level/100)/2),2,1,0)
            R.cis_diagnostics[o,.]=(o,0,0,0,0,.,1-level/100,1-level/100,level/100,1-level/100,1,1)
        }
        return(R)
    }
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,gamma,outputs)
    Q=jw_mcc_from_plan(P,M.V,M.z,F.A[outputs,.],level,gamma,mcerr,reps,
        mode,nodes,search,outerreps,validreps,steps,certcells,radius,pmax,mcbeta,progress,
        isreps,isstarts,isiter)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    R.calibration=Q.calibration
    R.domain_radius=Q.domain_radius
    R.exterior_radius=Q.exterior_radius
    R.cis_diagnostics[outputs,.]=Q.cis_diagnostics
    R.cis_points=J(nout,P.d,.)
    R.cis_points[outputs,.]=Q.cis_points
    failed=J(1,0,.)
    for (i=1;i<=cols(outputs);i++) if (Q.diagnostics[i,14]>=3) failed=(failed,outputs[i])
    if (cols(failed)>0) {
        if (fallback=="error") jw_core_error("conditional McCloskey calibration did not validate; increase budget or use mcfallback(bb)")
        alpha0=1-level/100-mcerr/2
        if (progress) printf("McCloskey CIS: certified BB fallback for %g outputs; not calibrated-M success\n",cols(failed))
        BB=jw_infer_bb_cis(M,C,F,joint,scope,maxjoint,100*(1-alpha0),gamma,
            mcerr/2,reps,ncells,.05,failed,progress,"certified",isreps,isstarts,isiter)
        R.bb_diagnostics=BB.diagnostics
        for (i=1;i<=cols(failed);i++) {
            o=failed[i]
            status=R.diagnostics[o,14]
            R.table[o,.]=BB.table[o,.]
            R.table[o,5]=.
            R.diagnostics[o,13]=BB.diagnostics[o,10]
            R.diagnostics[o,14]=status==4 ? 6 : 5
            R.diagnostics[o,15]=0
        }
    }
    return(R)
}

/* ===== jwselect_wbic_pointwise.mata ===== */
/* Equal-tailed Gaussian reselection inference in sampling-covariance units.
   The bootstrap draws H* around observed h and repeats the WBIC selector.
   It subtracts the observed pretrend adjustment and studentizes
   by the variance of the candidate selected on each draw. */
struct jw_pw_result {
    real matrix table
    real matrix diagnostics
}

void jw_pw_controls(real scalar level, real scalar reps, string scalar engine)
{
    if (missing(level) | level<=0 | level>=100 |
        missing(reps) | reps<1 | reps!=floor(reps))
        jw_core_error("invalid level() or reps() for WBIC reselection inference")
    if (engine!="direct" & engine!="integrated")
        jw_core_error("unknown pointwise WBIC simulation engine")
}

/* For a continuous finite normal mixture, the smallest/largest component
   p-quantiles bracket its p-quantile even when studentizers vary by draw. */
real scalar jw_pw_mixture_quantile(real colvector location,
                                    real colvector spread, real scalar prob)
{
    real colvector component
    real scalar lower,upper,middle,i
    component=location:+invnormal(prob)*spread
    lower=min(component)
    upper=max(component)
    for (i=1;i<=64;i++) {
        middle=lower/2+upper/2
        if (middle<=lower | middle>=upper) break
        if (mean(normal((middle:-location):/spread))<prob) lower=middle
        else upper=middle
    }
    return(upper)
}

/* A caller with its own complete action plan (e.g. an aggregate of separate
   cohort MD fits) may reuse this kernel. Check the action actually selected
   at observed h against the point estimator that will be reported. */
struct jw_pw_result scalar jw_pw_from_plan(struct jw_wbic_plan scalar P,
    real matrix V, real colvector z, real matrix observedA,
    real scalar level, real scalar reps, string scalar engine)
{
    struct jw_pw_result scalar R
    struct jw_mc_state scalar S
    struct jw_mc_draw scalar D
    real matrix root,slopes,E,base,draws,spread,values,candidateA
    real rowvector h,adjustment,residualvar,residualsd
    real colvector ordered
    real scalar total,ns,first,last,n,o,chosen,lo,hi,tail,est,se
    jw_pw_controls(level,reps,engine)
    total=rows(z)
    ns=rows(P.masks)
    if (P.d<1 | P.p+P.k!=total | P.nout<1 |
        cols(P.outputs)!=P.nout | rows(V)!=total | cols(V)!=total |
        rows(observedA)!=P.nout | cols(observedA)!=total |
        rows(P.actions)!=P.nout*ns | cols(P.actions)!=total |
        rows(P.sigmas)!=ns | cols(P.sigmas)!=P.nout |
        hasmissing(z) | hasmissing(V) | hasmissing(observedA) |
        hasmissing(P.actions) | hasmissing(P.sigmas) | min(P.sigmas)<=0)
        jw_core_error("incompatible WBIC reselection inputs")
    S=jw_mc_prepare(P,V)
    h=z[(P.p+1)..total]'
    D=jw_mc_draws(P,S,J(1,total,0),h)
    chosen=D.chosen[1]
    adjustment=D.adjustment
    for (o=1;o<=P.nout;o++) {
        candidateA=P.actions[(o-1)*ns+chosen,.]
        if (max(abs(candidateA-observedA[o,.]))>1e-8*(1+max(abs(observedA[o,.]))))
            jw_core_error("observed reselection action differs from the WBIC estimator")
    }

    draws=J(reps,P.nout,.)
    residualsd=J(1,P.nout,0)
    if (engine=="direct") {
        root=jw_infer_sqrt(V)
    }
    else {
        slopes=cholsolve(S.VH,S.cross)
        residualvar=S.basevar-colsum(S.cross:*slopes)
        if (hasmissing(residualvar) | any(residualvar:<-1e-10*S.basevar))
            jw_core_error("nonpositive-semidefinite conditional covariance in WBIC reselection")
        for (o=1;o<=P.nout;o++) residualsd[o]=sqrt(max((0,residualvar[o])))
        spread=J(reps,P.nout,.)
    }
    tail=(1-level/100)/2
    if ((engine=="direct" | any(residualsd:==0)) & reps*tail < 1-1e-12)
        jw_core_error("increase reps() to include an empirical draw in each requested WBIC tail")
    /* Bound temporary primitive noise and action matrices. Shared Gaussian
       draws preserve the covariance across selection, estimates and cohorts. */
    for (first=1;first<=reps;first=last+1) {
        last=min((reps,first+255))
        n=last-first+1
        if (engine=="direct") {
            E=rnormal(n,total,0,1)*root'
            base=E*S.baseline'
        }
        else {
            E=J(n,total,0)
            E[.,P.p:+S.active]=rnormal(n,P.d,0,1)*S.Hroot'
            base=E[.,P.p:+S.active]*slopes
        }
        D=jw_mc_draws(P,S,E,h)
        values=base+D.adjustment
        values=(values-J(n,1,1)*adjustment):/D.sigma
        draws[|first,1\last,P.nout|]=values
        if (engine=="integrated") {
            values=J(n,1,1)*residualsd
            values=values:/D.sigma
            spread[|first,1\last,P.nout|]=values
        }
    }
    if (hasmissing(draws)) jw_core_error("nonfinite signed WBIC reselection root")
    R.table=J(P.nout,5,.)
    R.diagnostics=J(P.nout,8,.)
    for (o=1;o<=P.nout;o++) {
        if (engine=="direct" | residualsd[o]==0) {
            /* Left empirical quantiles. At zero conditional variance the
               integrated mixture consists entirely of these point masses. */
            ordered=sort(draws[.,o],1)
            lo=ordered[max((1,ceil(reps*tail)))]
            hi=ordered[min((reps,ceil(reps*(1-tail))))]
        }
        else {
            lo=jw_pw_mixture_quantile(draws[.,o],spread[.,o],tail)
            hi=jw_pw_mixture_quantile(draws[.,o],spread[.,o],1-tail)
        }
        est=observedA[o,.]*z
        se=sqrt(jw_fixed_variance(V,observedA[o,.]))
        R.table[o,.]=(est,se,est-hi*se,est-lo*se,.)
        R.diagnostics[o,.]=(P.outputs[o],P.d,ns,reps,lo,hi,adjustment[o],0)
    }
    return(R)
}

struct jw_pw_result scalar jw_infer_wbic_pointwise(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level, real scalar reps,
    real rowvector outputs, string scalar engine)
{
    struct jw_pw_result scalar R,Q
    struct jw_wbic_plan scalar P
    real rowvector active
    real scalar i,o,nout,critical
    jw_pw_controls(level,reps,engine)
    nout=rows(F.A)
    if (cols(outputs)==0) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout |
        any(outputs:!=floor(outputs))) jw_core_error("invalid WBIC reselection output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,8,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,1]=o
        R.diagnostics[o,8]=9
    }
    active=J(1,0,.)
    if (M.k>0) {
        if (scope=="joint" & rows(joint)>0) active=selectindex(colsum(joint):>0)
        else active=selectindex(colsum(C.masks):>0)
    }
    if (cols(active)==0) {
        critical=invnormal(1-(1-level/100)/2)
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,.]=jw_infer_normal(R.table[o,1],R.table[o,2]^2,level)
            R.table[o,5]=.
            R.diagnostics[o,.]=(o,0,1,0,-critical,critical,0,2)
        }
        return(R)
    }
    /* The shared plan takes an ellipsoid level for its BB-specific fields.
       Pointwise methods read only its actions, selector and variances. */
    P=jw_wbic_prepare(M,C,F,joint,scope,maxjoint,.05,outputs)
    Q=jw_pw_from_plan(P,M.V,M.z,F.A[outputs,.],level,reps,engine)
    R.table[outputs,.]=Q.table
    R.diagnostics[outputs,.]=Q.diagnostics
    return(R)
}

/* ===== jwselect_grid_io.mata ===== */
/* Private exchange format for the native BB fallback. All filenames
   and Stata matrix names come from caller-owned tempfiles and tempnames. */
void jw_grid_write_input(string scalar path, string rowvector names)
{
    real scalar fh,i,r
    real matrix A
    string rowvector labels
    labels=("h","V","c","actions","masks","sigma","penalty","control","score_cov")
    fh=fopen(path,"w")
    for (i=1;i<=cols(labels);i++) {
        A=st_matrix(names[i])
        fput(fh,labels[i]+" "+strofreal(rows(A))+" "+strofreal(cols(A)))
        for (r=1;r<=rows(A);r++) fput(fh,invtokens(strofreal(A[r,.],"%25.17e")))
    }
    fclose(fh)
}

/* ===== jwselect_bb_grid.mata ===== */
/* Adaptive BB integration. Cohort-wise selection still fits joint MD on
   the union of the selected restrictions. A cohort is removed only when
   its scores cannot affect selection of retained cohorts and its entire
   covariance block is independent of both their restrictions and B. */
struct jw_grid_result {
    real matrix table
    real matrix diagnostics
    real matrix grid
}

real rowvector jw_grid_relevant(struct jw_model scalar M, real rowvector active,
    real colvector cross, string scalar scope)
{
    real rowvector keep, ix, groups
    real scalar j, old
    if (scope=="joint") return(active)
    groups=M.hmap[active,1]'
    keep=J(1,cols(active),0)
    for (j=1;j<=cols(active);j++) if (cross[j]!=0)
        keep=keep:|(groups:==groups[j])
    if (!sum(keep)) return(J(1,0,.))
    do {
        old=sum(keep)
        ix=M.p:+active[selectindex(keep)]
        for (j=1;j<=cols(active);j++) if (!keep[j]) {
            if (any(M.V[M.p+active[j],ix]:!=0)) keep=keep:|(groups:==groups[j])
        }
    } while (sum(keep)>old)
    return(active[selectindex(keep)])
}

/* Resource fallback without Cartesian enumeration. With a whitened
   baseline regression b and conditional residual SD s, every efficient
   MD action has ||adjustment|| <= ||b|| and SE >= s. All its centered
   studentized errors are unit projections of a single N(0,I_(d+1)).
   Both a chi bound and a log-cardinality Bonferroni bound apply. */
real rowvector jw_grid_projection_cap(real matrix V, real colvector h,
    real colvector c, real scalar v, real scalar radius, real scalar tail,
    real scalar logcount, real rowvector groups, real scalar maxpen)
{
    real scalar gain,s2,normh,z,upper,lower,baserisk,basecv,j,lambda,K
    real matrix q,L,block,D
    q=cholsolve(V,(h,c))
    if (hasmissing(q)) jw_core_error("adaptive BB cannot bound a singular non-enumerable family; reduce its candidates or increase maxjoint()")
    normh=sqrt(max((0,h'*q[.,1])))
    gain=max((0,c'*q[.,2]))
    s2=v-gain
    if (s2<=1e-12*v) jw_core_error("adaptive BB needs a nondegenerate target for a non-enumerable family; reduce its candidates or increase maxjoint()")
    z=min((sqrt(invchi2(rows(h)+1,1-tail)),sqrt(2*(ln(2/tail)+logcount))))
    upper=(z+(normh+radius)*sqrt(gain/s2))*(1+1e-10)+1e-10
    lower=0
    baserisk=.
    /* No union listing is needed for the following selection bound.
       Conditional on winning, the sum of retained owner-Wald statistics
       is <= maxpen.  If D=blockdiag(V_gg), Cauchy-Schwarz gives
       |a'H|^2 <= maxpen * a'Da. For every MD action a'Va<=gain.
       A Gershgorin bound on L^-1 D L^-T (V=LL') controls a'Da/a'Va. */
    L=cholesky(V)
    if (!hasmissing(L)) {
        block=V
        for (j=1;j<=rows(V);j++) block[j,.]=block[j,.]:*(groups:==groups[j])
        D=lusolve(L,block)
        D=lusolve(L,D')
        if (!hasmissing(D)) {
            lambda=max(rowsum(abs((D+D')/2)))*(1+1e-10)+1e-10
            K=sqrt(lambda*gain*maxpen)*(1+1e-10)+1e-10*sqrt(v)
            basecv=invnormal(1-tail/2)
            z=(basecv*sqrt(v)+K)/sqrt(s2)*(1+1e-10)+1e-10
            if (z<upper) {
                upper=z
                baserisk=tail
            }
            lower=max((0,basecv-K/sqrt(v)-1e-10))
        }
    }
    return((lower,upper,baserisk))
}

/* Dimension-independent global BB upper endpoint. For action j, E=B+a'H
   and a'H are independent at every fixed drift (efficient joint MD).
   Winning requires |a'H|<=K_j. The drift magnitude |a'mu| varies over
   an interval given by the support of the BB ellipsoid. On each scalar
   interval, the rejection tail is largest at the right endpoint and the
   retention probability is largest at the left. Summing the products of
   these bounds across possible winners protects the entire confidence set.
   Bounds are expressed in units of each fixed action's standard error. */
real rowvector jw_grid_action_profile(real matrix V, real colvector h,
    real colvector c, real scalar v, real matrix A, real matrix masks, real colvector s,
    real colvector pens, real rowvector groups, string scalar scope,
    real scalar rad, real scalar tail, real scalar cap, real scalar cells)
{
    real matrix left,right,ret,fold
    real colvector b,d0,gain,k,lo,hi,tv,r,nr,risk,gg
    real rowvector ix
    real scalar j,g,n,bins,a,bound,mid,worst,certifiedrisk,penunit,
        basecv,baseupper,baselower
    n=rows(A)
    gain=rowsum((A*V):*A)
    b=sqrt(rowmax((gain,J(n,1,0))))
    if (scope=="joint") k=b:*sqrt(pens)
    else {
        k=J(n,1,0)
        gg=uniqrows(sort(groups',1))
        penunit=max(pens:/rowmax((rowsum(masks),J(n,1,1))))
        for (g=1;g<=rows(gg);g++) {
            ix=selectindex(groups:==gg[g])
            nr=rowsum(masks[.,ix])
            r=sqrt(rowmax((rowsum((A[.,ix]*V[ix,ix]):*A[.,ix]),J(n,1,0))))
            /* Every nonempty model uses the common penalty per restriction. */
            k=k:+r:*sqrt(nr*penunit)
        }
    }
    /* Roundoff in the residual cross-covariance can be covered by TV
       between the actual bivariate law and the product of its marginals. */
    ix=selectindex(b:>0)'
    tv=J(n,1,0)
    if (cols(ix)) {
        r=abs(A[ix,.]*c+gain[ix]):/(s[ix]:*b[ix])
        r=rowmin((r:+1e-12,J(rows(ix'),1,1)))
        tv[ix]=sqrt(-ln(1:-r:^2):/4)
        tv=editmissing(tv,1)
    }
    d0=abs(A*h):/s
    b=b:/s
    k=(k:/s):*(1+1e-10):+1e-10
    lo=rowmax((J(n,1,0),d0-rad*b:*(1+1e-10):-1e-10))
    hi=(d0+rad*b):*(1+1e-10):+1e-10
    /* The same Gaussian baseline error applies to every winner. Triangle
       inequalities alone give two endpoints valid for all nuisance drifts.
       These can be much tighter than a many-candidate union bound when
       the estimated adjustment is small or nearly target-orthogonal. */
    basecv=invnormal(1-tail/2)*sqrt(v)
    baseupper=max(basecv:/s+k)*(1+1e-10)+1e-10
    baselower=max((0,min(basecv:/s-k)*(1-1e-10)-1e-10))
    bins=max((1,min((128,cells,floor(262144/n)))))
    left=lo*J(1,bins,1)+(hi-lo)*((0..(bins-1))/bins)
    right=lo*J(1,bins,1)+(hi-lo)*((1..bins)/bins)
    ret=J(n,bins,1)
    if (cols(ix)) {
        ret[ix,.]=normal((k[ix]:-left[ix,.]):/b[ix])-normal((-k[ix]:-left[ix,.]):/b[ix]):+2e-10
    }
    a=0
    bound=min((cap,baseupper))
    certifiedrisk=baseupper<=cap ? tail : .
    for (j=1;j<=32;j++) {
        mid=(a+bound)/2
        fold=normal(right:-mid)+normal(-right:-mid):+2e-10
        risk=rowmax(fold:*ret):+tv
        worst=sum(rowmin((risk,J(n,1,1))))+2e-10
        if (worst<=tail) {
            bound=mid
            certifiedrisk=worst
        }
        else a=mid
        if (bound-a<=1e-7*(1+bound)) break
    }
    return((bound,certifiedrisk,baselower))
}

/* A finite-witness LOWER endpoint for the portable/profile routes. Replay
   the true owner scores, then use the joint-MD loading for their selected
   union. Even a Cartesian family too large to list can be replayed this
   way. Witnesses and draw counts are fixed before simulation. At a fixed
   witness |T|_(k) is a lower bound on its (1-tail) quantile except with
   probability BinomialCDF(n,k-1,1-tail). The three error allocations cover
   the maximum of these lower bounds even with common random numbers. */
real rowvector jw_grid_fallback_witness(real matrix V, real colvector h,
    real colvector c, real scalar v, real matrix masks, real colvector owners,
    real scalar logN, real scalar radius, real scalar tail, real scalar error,
    real scalar reps, real colvector direction)
{
    real matrix L,noise,X,chosen,Q
    real colvector b,base,best,choice,score,piv,ord,starts,ends,who,x,gids,ii,ok,y
    real rowvector ix
    real scalar d,n,j,g,i,nb,first,last,s2,sd,lo,hi,k,mid,lower,point,fac
    d=rows(h)
    /* Limit the cost of scoring actual candidates (not their nominal
       Cartesian product). A small sample without a useful rank is skipped. */
    n=min((reps,4096,floor(3e7/(3*(d^2+sum(rowsum(masks):^2)+rows(masks))))))
    if (n<32) return((0,0,.))
    L=cholesky(V)
    b=cholsolve(V,c)
    if (hasmissing(L) | hasmissing(b)) return((0,0,.))
    s2=v-c'*b
    if (s2<=1e-12*v) return((0,0,.))
    /* Conservative order statistic with exact binomial rank calculation. */
    lo=0
    hi=n+1
    while (hi-lo>1) {
        mid=floor((hi+lo)/2)
        if (binomial(n,mid-1,1-tail)<=error/3) lo=mid
        else hi=mid
    }
    k=lo
    if (k==0) return((0,0,.))
    sd=sqrt(max((0,direction'*V*direction)))
    if (sd>0) {
        fac=direction'*h<0 ? -1 : 1
        direction=fac*(1-1e-10)*radius*V*direction/sd
        /* Verify a computed boundary witness is actually in the ellipsoid.
           If cancellation or an ill-conditioned metric prevents this, the
           center alone remains genuinely feasible. */
        x=(h+direction)-h
        if (x'*cholsolve(V,x)>radius^2) direction=J(d,1,0)
        x=(h+.5*direction)-h
        if (x'*cholsolve(V,x)>radius^2) direction=J(d,1,0)
    }
    else direction=J(d,1,0)
    noise=rnormal(n,d,0,1)*L'
    base=noise*b+sqrt(s2)*rnormal(n,1,0,1)
    X=noise:+h'
    X=X\(noise:+(h+.5*direction)')\(noise:+(h+direction)')
    base=base\base\base
    nb=rows(X)
    chosen=J(nb,d,0)
    gids=uniqrows(sort(owners,1))
    for (g=1;g<=rows(gids);g++) {
        ii=vec(selectindex(owners:==gids[g]))
        best=J(nb,1,.)
        choice=J(nb,1,ii[1])
        for (i=1;i<=rows(ii);i++) {
            j=ii[i]
            ix=selectindex(masks[j,.])
            if (cols(ix)>0) {
                Q=cholsolve(V[ix,ix],I(cols(ix)))
                if (hasmissing(Q)) return((0,0,.))
                score=rowsum((X[.,ix]*Q):*X[.,ix]):-cols(ix)*logN
            }
            else score=J(nb,1,0)
            ok=vec(selectindex(score:<best))
            if (rows(ok)) {
                best[ok]=score[ok]
                choice[ok]=J(rows(ok),1,j)
            }
        }
        chosen=chosen:|masks[choice,.]
    }
    ord=order(chosen,1..d)
    chosen=chosen[ord,.]
    starts=1\(vec(selectindex(rowsum(chosen[2..nb,.]:!=chosen[1..(nb-1),.]):>0)):+1)
    ends=starts:-1
    ends=(ends\nb)[2..(rows(ends)+1)]
    piv=J(nb,1,.)
    for (i=1;i<=rows(starts);i++) {
        first=starts[i]
        last=ends[i]
        who=ord[first..last]
        ix=selectindex(chosen[first,.])
        if (!cols(ix)) piv[who]=abs(base[who])/sqrt(v)
        else {
            b=cholsolve(V[ix,ix],c[ix])
            sd=sqrt(v-c[ix]'*b)
            if (missing(sd) | sd<=0) return((0,0,.))
            piv[who]=abs(base[who]-X[who,ix]*b)/sd
        }
    }
    lower=point=0
    for (j=1;j<=3;j++) {
        y=sort(piv[((j-1)*n+1)..(j*n)],1)
        lower=max((lower,y[k]))
        point=max((point,y[ceil((1-tail)*n)]))
    }
    return((max((0,lower*(1-1e-10)-1e-10)),nb,point))
}

real rowvector jw_grid_callback(real colvector h, real matrix V, real colvector c,
    real matrix actions, real matrix masks, real colvector sigmas, real colvector pens,
    real rowvector controls, real matrix scorecov, real scalar route)
{
    string rowvector n
    real scalar rc
    real rowvector result
    n=st_tempname(10)
    st_matrix(n[1],h)
    st_matrix(n[2],V)
    st_matrix(n[3],c)
    st_matrix(n[4],actions)
    st_matrix(n[5],masks)
    st_matrix(n[6],sigmas)
    st_matrix(n[7],pens)
    st_matrix(n[8],controls)
    st_matrix(n[10],scorecov)
    rc=_stata("_jwselect_grid_call "+invtokens(n))
    if (rc==497 | rc==498 | rc==499) {
        route=(rc==497 ? 8 : (rc==498 ? 5 : 6))
        return((0,controls[6],.,0,0,0,.,4,1))
    }
    if (rc) _error(rc)
    result=st_matrix(n[9])
    if (cols(result)!=9 | rows(result)!=1) jw_core_error("invalid adaptive BB backend result")
    if (hasmissing(result[(1,2,4,5,6,8,9)]) | result[1]<0 | result[2]<result[1] |
        result[2]>controls[6]+1e-9 | min(result[(4,5,6)])<0 |
        !any(result[8]:==(0,4)) | !any(result[9]:==(0,1)))
        jw_core_error("invalid adaptive BB backend bracket")
    if (!result[9] & (missing(result[7]) | result[7]>controls[3] | result[7]<0))
        jw_core_error("adaptive BB global risk does not clear the protected tail")
    if (result[9] & abs(result[2]-controls[6])>1e-10*(1+controls[6]))
        jw_core_error("adaptive BB fallback endpoint disagrees with the analytic cap")
    return(result)
}

struct jw_grid_result scalar jw_infer_bb_grid(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level, real scalar gamma,
    real scalar mcerr, real scalar reps, real scalar cells, real scalar tol,
    real rowvector outputs, real scalar progress, real scalar pilotreps, real scalar seed)
{
    struct jw_grid_result scalar R
    struct jw_model scalar restricted
    real matrix masks,jointset,V,sv,A,add
    real rowvector target,active,keep,groups,control,res,ix,bestmask,profile
    real colvector h,c,allc,sigmas,pens,support,ci,wdir,widx
    real rowvector witness
    real scalar alpha,tail,nout,d,d0,i,o,j,g,v,se,est,cap,cv,rad,ns,
        route,status,logcount,canlist,score,best,winner,limit,fullcount,maxpen
    alpha=1-level/100
    tail=alpha-gamma-mcerr
    if (missing(tail) | tail<=0 | gamma<=0 | mcerr<=0)
        jw_core_error("gamma()+mcerr() must be smaller than 1-level()/100")
    if (scope=="cohort" & rows(joint)>0) jw_core_error("cohort WBIC cannot use jointmatrix()")
    nout=rows(F.A)
    if (!cols(outputs)) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs)))
        jw_core_error("invalid adaptive BB output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.grid=J(nout,12,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,(1,11)]=(o,9)
        R.grid[o,(1,11,12)]=(o,9,9)
    }
    jointset=J(0,M.k,.)
    active=J(1,0,.)
    if (M.k) {
        if (scope=="joint") {
            jointset=jw_joint_masks(M,C,joint,maxjoint)
            active=selectindex(colsum(jointset):>0)
        }
        else active=selectindex(colsum(C.masks):>0)
    }
    d0=cols(active)
    if (d0>0) rad=sqrt(invchi2(d0,1-gamma))
    for (i=1;i<=cols(outputs);i++) {
        o=outputs[i]
        target=F.targets[o,.]
        est=R.table[o,1]
        se=R.table[o,2]
        keep=J(1,0,.)
        if (d0) {
            allc=M.V[M.p:+active,.]*target'
            keep=jw_grid_relevant(M,active,allc,scope)
        }
        d=cols(keep)
        if (!d) {
            if (max(abs(F.A[o,.]-target))>1e-10*(1+max(abs(target))))
                jw_core_error("zero-restriction grid target differs from selected WBIC action")
            /* A truly empty family uses the existing ordinary normal CI.
               If selection still occurs but this target is independent of
               every active restriction, invert its exact law at the BB tail. */
            cv=invnormal(1-(d0 ? tail : alpha)/2)
            R.table[o,.]=(est,se,est-cv*se,est+cv*se,.)
            R.diagnostics[o,.]=(o,0,1,0,0,(d0 ? tail : alpha),cv,cv,0,cv,2,1,1,0,0,cv)
            R.grid[o,.]=(o,0,1,0,0,0,.,cv,cv,0,2,2)
            continue
        }
        if (se<=0 | missing(se)) jw_core_error("adaptive BB requires positive selected standard error")
        ix=M.p:+keep
        h=M.z[ix]
        V=M.V[ix,ix]
        c=M.V[ix,.]*target'
        v=jw_fixed_variance(M.V,target)
        groups=M.hmap[keep,1]'
        canlist=1
        logcount=0
        fullcount=1
        maxpen=0
        if (scope=="joint") masks=jointset
        else {
            restricted=M
            restricted.groups=select(M.groups,rowsum(M.groups*J(1,d,1):==J(rows(M.groups),1,1)*groups):>0)
            for (g=1;g<=rows(restricted.groups);g++) {
                ci=vec(selectindex(C.owners:==restricted.groups[g]))
                ns=rows(ci)
                maxpen=maxpen+max(rowsum(C.masks[ci,keep]))*ln(M.nunit)
                logcount=logcount+ln(ns)
                if (fullcount>floor(maxjoint/ns)) canlist=0
                fullcount=fullcount*ns
            }
            if (canlist) masks=jw_joint_masks(restricted,C,J(0,M.k,.),maxjoint)
        }
        if (!canlist) {
            ns=fullcount
            profile=jw_grid_projection_cap(V,h,c,v,rad,tail,logcount,groups,maxpen)
            cap=profile[2]
            route=4
            res=(profile[1],cap,.,0,0,0,profile[3],4,1)
            wdir=cholsolve(V,c)
        }
        else {
            ns=rows(masks)
            A=J(ns,rows(M.z),.)
            sigmas=J(ns,1,.)
            support=J(ns,1,.)
            sv=V
            if (scope=="cohort") {
                for (j=1;j<=d;j++)
                    sv[j,.]=sv[j,.]:*(groups:==groups[j])
            }
            best=.
            winner=1
            for (j=1;j<=ns;j++) {
                A[j,.]=jw_loading(M.V,target,masks[j,.],M.k)
                sigmas[j]=sqrt(jw_fixed_variance(M.V,A[j,.]))
                if (sigmas[j]<=0 | missing(sigmas[j]))
                    jw_core_error("adaptive BB requires nondegenerate variance for every candidate target")
                bestmask=masks[j,keep]
                ix=selectindex(bestmask)
                score=0
                if (cols(ix)) score=h[ix]'*cholsolve(sv[ix,ix],h[ix])-cols(ix)*ln(M.nunit)
                if (score<best) {
                    best=score
                    winner=j
                }
            }
            add=A[.,M.p:+keep]:-J(ns,1,1)*target[M.p:+keep]
            support=abs(add*h)+rad*sqrt(rowmax((J(ns,1,0),rowsum((add*V):*add))))
            widx=vec(selectindex(support:/sigmas:==max(support:/sigmas)))
            wdir=add[widx[1],.]'
            limit=invnormal(1-tail/(2*ns))
            if (missing(limit)) limit=sqrt(2*ln(2*ns/tail))
            cap=(max(support:/sigmas)+limit)*(1+1e-10)+1e-10
            pens=rowsum(masks)*ln(M.nunit)
            route=1
            /* Predicted memory for the actual number of corners and models,
               not the number of preperiods alone, limits vertex replay. */
            if (d*ln(2)>ln(262144) | d*ln(2)+ln(d+ns)>ln(2e6) |
                max((d,ns))>strtoreal(st_global("c(max_matdim)"))) route=8
            if (max(abs(A[winner,.]-F.A[o,.]))>1e-8*(1+max(abs(F.A[o,.])))) route=7
            profile=(cap,.,0)
            if (route!=7) {
                profile=jw_grid_action_profile(V,h,c,v,add,masks[.,keep],sigmas,pens,
                    groups,scope,rad,tail,cap,cells)
                cap=profile[1]
                if (route==1 & cap-profile[3]<=tol) route=8
            }
            if (route==1) {
                control=(v,rad,tail,mcerr/cols(outputs),tol,cap,mod(seed+o-1,2147483647),cells,reps,pilotreps)
                res=jw_grid_callback(h,V,c,add,masks[.,keep],sigmas,pens,control,sv,route)
            }
            else res=(0,cap,.,0,0,0,.,4,1)
            res[1]=max((res[1],profile[3]))
            if (res[9] & res[5]==0 & profile[2]<.) {
                res[7]=profile[2]
                if (route==1) route=8
            }
        }
        if ((route==4 | route==8) & !res[6] & !hasmissing(wdir) & res[2]-res[1]>tol) {
            if (scope=="cohort") {
                widx=vec(selectindex(rowsum(C.owners*J(1,cols(groups),1):==
                    J(rows(C.owners),1,1)*groups):>0))
                witness=jw_grid_fallback_witness(V,h,c,v,C.masks[widx,keep],C.owners[widx],
                    ln(M.nunit),rad,tail,mcerr/cols(outputs),reps,wdir)
            }
            else witness=jw_grid_fallback_witness(V,h,c,v,masks[.,keep],J(rows(masks),1,1),
                    ln(M.nunit),rad,tail,mcerr/cols(outputs),reps,wdir)
            /* A rare discrepancy between a stochastic lower and a proved
               upper invalidates the witness, never the protected interval. */
            if (witness[1]<=cap) res[1]=max((res[1],witness[1]))
            res[6]=witness[2]
            res[3]=witness[3]
        }
        cv=res[2]
        status=(res[2]-res[1]<=tol ? 0 : 4)
        R.table[o,.]=(est,se,est-cv*se,est+cv*se,.)
        R.diagnostics[o,.]=(o,d,ns,res[5],res[6],tail,res[1],cv,cv-res[1],cv,
            status,1,status==0,res[9],.,res[3])
        R.grid[o,.]=(o,d,ns,res[4],res[5],res[6],res[7],res[1],cv,cv-res[1],status,route)
        if (progress) printf("Adaptive BB output %g: [%g, %g], %g cells, route %g\n",o,res[1],cv,res[4],route)
    }
    return(R)
}

/* ===== jwselect_mc_band.mata ===== */
/* Optional second stage for the verified conditional-score rule. Band
   limits and the nuisance cover depend only on the Gaussian plan, never on
   observed pretrends. All widths and baseline draws below use baseline-SD
   units. The return is a global risk upper, cells visited and node-draws. */
real rowvector jw_mcp_band(real matrix V, real colvector c, real scalar v,
    real matrix masks, real rowvector owner, string scalar scope, real scalar logn,
    real scalar B, real scalar target, real scalar beta, real scalar p, real scalar seed,
    real scalar maxcells, real scalar n, real scalar delta)
{
    struct jw_wbic_score colvector score
    real scalar d,ns,j,k,box,nb,R,tau,lower,top,outer,mass,allow,
        zbeta,guardcv,rad,ub,visited,draws,worst,depth,axis,i,tailrate,minupper
    real matrix L,sv,A,C,Q,eps,H,X,boxes,kids,Z
    real colvector u,q,ix,sd,b,stat,signedb,absb,selected,best,change,
        x,cutoff,foldmean,reject,pens,limits,pass
    real rowvector g,lo,hi,nearest,mid,half,ctr,dist,child
    d=rows(V)
    ns=rows(masks)
    if (maxcells<1 | n<2 | delta<=0 | target<=0 | p>=1 | v<=0) return((1,0,0))
    tau=target/2
    zbeta=invnormal(1-beta/2)
    guardcv=invnormal(1-tau/2)
    top=min((guardcv,invnormal((1+p)/2)))-1e-10
    outer=2*normal(-top)+1e-12
    if (outer>=target | top<=0 | n>2e6) return((1,0,0))
    /* Beyond R, entry into the observation ball has probability <=tau by
       projection onto the nuisance direction. Guard rejection costs <=tau. */
    R=B+invnormal(1-tau)+1e-9
    L=cholesky(V)
    if (hasmissing(L)) return((1,0,0))
    u=lusolve(L,c)/sqrt(v)
    if (hasmissing(u) | sum(u:^2)>=1-1e-12) return((1,0,0))
    C=cholesky(I(d)-u*u')
    if (hasmissing(C)) return((1,0,0))
    Q=I(d)
    if (sum(u:^2)>0) {
        q=Q[.,1]-u/sqrt(sum(u:^2))
        if (sum(q:^2)>1e-24) Q=Q-2*q*q'/sum(q:^2)
    }
    /* The absolute row sums dominate the full metric's quadratic form,
       also with residual off-diagonal terms after rotation. */
    q=Q'*u
    g=J(1,d,1)+sum(abs(q))*abs(q')/(1-sum(u:^2))
    sv=V
    if (scope=="cohort") for (j=1;j<=d;j++) sv[j,.]=sv[j,.]:*(owner:==owner[j])
    if (ns*d*d>5e6 | ns*n>1e8 | n*d>5e6) return((1,0,0))
    score=jw_wbic_score(ns)
    A=J(ns,d,0)
    pens=rowsum(masks)*logn
    for (j=1;j<=ns;j++) {
        ix=vec(selectindex(masks[j,.]))
        score[j].idx=ix'
        if (rows(ix)) {
            A[j,ix]=-cholsolve(V[ix,ix],c[ix])'/sqrt(v)
            score[j].W=luinv(cholesky(sv[ix,ix]))
            if (hasmissing(A[j,ix]) | hasmissing(score[j].W)) return((1,0,0))
        }
    }
    b=sqrt(rowmax((J(ns,1,0),rowsum((A*V):*A))))
    if (max(b:^2)>=1) return((1,0,0))
    sd=sqrt(1:-b:^2)
    limits=zbeta*b:+invnormal(p)*sd
    for (j=1;j<=ns;j++) if (!sum(masks[j,.])) limits[j]=invnormal((1+p)/2)
    lower=max((0,min((guardcv,min(limits)))-1e-10))
    mass=2*(normal(-lower)-normal(-top))+1e-12
    if (2*normal(-lower)+1e-10<=target) return((2*normal(-lower)+1e-10,0,0))
    if (mass<=0) return((1,0,0))
    allow=(target-outer)/mass
    if (allow<=0 | allow>=1) return((1,0,0))
    minupper=jw_cis_kl_upper(0,n,delta/maxcells)
    rad=invnormal(allow)-invnormal(minupper)
    /* Even the tightest possible center bound needs at least (R/rad)^d
       enclosing balls by volume. This is only a budget rejection. */
    if (rad<=0 | missing(rad)) return((1,0,0))
    if (d*ln(R/rad)>ln(maxcells)+1e-8) return((1,0,0))
    rseed(seed)
    boxes=(-J(1,d,R),J(1,d,R))
    visited=draws=worst=0
    for (depth=0;depth<80;depth++) {
        nb=rows(boxes)
        if (visited+nb>maxcells | nb*2*d>5e6) break
        /* Independent from every preceding partition decision. Across
           centers of this stage the same draws may be used. */
        absb=invnormal(normal(lower):+(normal(top)-normal(lower))*runiform(n,1))
        signedb=absb:*(2*(runiform(n,1):>.5):-1)
        eps=rnormal(n,d,0,1)*C':+signedb*u'
        Z=eps*Q
        X=eps*L'
        pass=J(nb,1,0)
        for (box=1;box<=nb;box++) {
            visited++
            lo=boxes[box,1..d]
            hi=boxes[box,d:+(1..d)]
            nearest=colmax((J(1,d,0)\lo\-hi))
            dist=sqrt(colmax((J(1,d,0)\(R^2-sum(nearest:^2):+nearest:^2)))):+1e-9
            lo=colmax((lo\-dist))
            hi=colmin((hi\dist))
            mid=(lo+hi)/2
            half=(hi-lo)/2
            ub=sqrt((half:^2)*g')+1e-8
            if (ub>=rad) continue
            ctr=mid*Q'*L'
            H=X:+J(n,1,1)*ctr
            best=J(n,1,.)
            selected=J(n,1,1)
            for (j=1;j<=ns;j++) {
                ix=score[j].idx'
                stat=J(n,1,-pens[j])
                if (rows(ix)) stat=stat+rowsum((H[.,ix]*score[j].W'):^2)
                change=vec(selectindex(stat:<best))
                if (rows(change)) {
                    best[change]=stat[change]
                    selected[change]=J(rows(change),1,j)
                }
            }
            x=rowsum(A[selected,.]:*H)
            cutoff=abs(signedb+x):/sd[selected]
            foldmean=(abs(x):+zbeta*b[selected]):/sd[selected]
            reject=((normal(foldmean-cutoff)+normal(-foldmean-cutoff)):<=1-p+1e-11):*
                (rowsum((Z:+J(n,1,1)*mid):^2):<=B^2+1e-8)
            draws=draws+n
            tailrate=jw_cis_kl_upper(mean(reject),n,delta/maxcells)
            tailrate=outer+mass*jw_mcr_np_upper(tailrate,ub)+1e-10
            if (tailrate<=target) {
                pass[box]=1
                worst=max((worst,tailrate))
            }
        }
        boxes=select(boxes,!pass)
        if (rows(boxes)==0) return((max((worst,target)),visited,draws))
        if (visited+rows(boxes)>maxcells | rows(boxes)*4*d>5e6) break
        kids=J(2*rows(boxes),2*d,.)
        k=0
        for (box=1;box<=rows(boxes);box++) {
            lo=boxes[box,1..d]
            hi=boxes[box,d:+(1..d)]
            half=(hi-lo)/2; mid=(hi+lo)/2
            ix=vec(selectindex((half:^2):*g:==max((half:^2):*g)))
            axis=ix[1]
            if (mid[axis]<=lo[axis] | mid[axis]>=hi[axis]) return((1,visited,draws))
            for (i=1;i<=2;i++) {
                child=boxes[box,.]
                if (i==1) child[d+axis]=mid[axis]
                else child[axis]=mid[axis]
                nearest=colmax((J(1,d,0)\child[1..d]\-child[d:+(1..d)]))
                if (sum(nearest:^2)<=R^2+1e-8) kids[++k,.]=child
            }
        }
        if (!k) break
        boxes=kids[1..k,.]
    }
    return((1,visited,draws))
}

/* ===== jwselect_mc_profile.mata ===== */
/* Full-space conditional-score calibration of the guarded McCloskey rule.
   All probabilities are marginal for a fixed target and Gaussian input.
   Calibration uses the covariance, family and tuning alone. The observed
   winner is used only afterwards to evaluate the protected width. */

real colvector jw_mcp_foldlo(real colvector mean, real colvector sd, real scalar p)
{
    real colvector lo,hi,mid,more,q1,q2,density
    real scalar j
    lo=rowmax((J(rows(sd),1,0),mean:+invnormal(p)*sd))
    if (p>=.5) {
        /* For q>=mean the folded-normal upper tail is decreasing and
           convex. Its tangent intersects the requested tail to the left
           of the true root. These short Newton iterates are lower bounds
           even without convergence. Pad tail/density outwards. */
        for (j=1;j<=6;j++) {
            q1=(mean-lo):/sd
            q2=(-mean-lo):/sd
            more=normal(q1)+normal(q2):-(1-p+2e-14)
            density=(normalden(q1)+normalden(q2)):*(1+1e-12):+1e-300
            lo=lo:+sd:*rowmax((more,J(rows(sd),1,0))):/density
        }
        return(rowmax((J(rows(sd),1,0),lo:-2e-11:*(1:+abs(lo)))))
    }
    hi=mean:+invnormal((1+p)/2)*sd
    for (j=1;j<=35;j++) {
        mid=(lo+hi)/2
        more=(normal((mean:-mid):/sd)+normal((-mean:-mid):/sd)):>(1-p+2e-14)
        lo=more:*mid+(1:-more):*lo
        hi=(1:-more):*mid+more:*hi
    }
    return(rowmax((J(rows(sd),1,0),lo:-2e-11:*(1:+abs(lo)))))
}

/* Each row (l,r) in curves represents Phi((m-l)/s0)+Phi((-m-r)/s0).
   Dominated rows can be removed exactly. If the Pareto frontier is wide,
   combine neighboring rows by lowering BOTH thresholds: this enlarges
   every affected rejection bound. Monotonicity in m>=0 is preserved. */
real scalar jw_mcp_integral(real matrix curves, real colvector m, real scalar s0)
{
    real colvector keep, val, envelope
    real scalar j,n,minright,step,r,k
    curves=sort(curves,(1,2))
    n=rows(curves)
    keep=J(n,1,0)
    minright=.
    for (j=1;j<=n;j++) {
        if (curves[j,2]<minright) {
            keep[j]=1
            minright=curves[j,2]
        }
    }
    curves=select(curves,keep)
    n=rows(curves)
    step=max((1,ceil(n/256)))
    envelope=J(rows(m),1,0)
    for (j=1;j<=n;j=j+step) {
        k=min((n,j+step-1))
        val=normal((m:-curves[j,1])/s0)+normal((-m:-curves[k,2])/s0)
        envelope=rowmax((envelope,val))
    }
    /* Upper Darboux sum on the half-normal probability scale.
       The last interval is charged one, including its infinite tail. */
    r=(sum(envelope)+1)/(rows(m)+1)
    return(min((1,r+2e-10)))
}

real scalar jw_mcp_upper(real colvector b, real colvector sd,
    real colvector cap, real matrix means, real scalar s0,
    real scalar rho, real scalar p, real scalar zbeta, real scalar guardcv,
    real scalar enumerated, real scalar target)
{
    real colvector lower,zero,w
    real matrix curves
    real scalar k,z,cheap,j
    z=invnormal(p)
    if (enumerated) {
        lower=z*sd:+zbeta*b
        for (j=1;j<=rows(b);j++) if (b[j]==0) lower[j]=invnormal((1+p)/2)*sd[j]
        k=min((guardcv,min(lower)))
    }
    else {
        /* For joint MD b in [0,rho], sd=sqrt(1-b^2). The minimum of
           z*sd+zbeta*b lies at an endpoint for z>=0, and at zero for z<0. */
        k=min((guardcv,z,z*s0+zbeta*rho))
    }
    cheap=min((1,2*normal(-k+1e-10)+2e-10))
    if (cheap<=target) return(cheap)
    if (enumerated) {
        zero=jw_mcp_foldlo(zbeta*b,sd,p)
        w=jw_mcp_foldlo(cap:+zbeta*b,sd,p)
        curves=((w-cap,w+cap) \ (min(zero),min(zero)))
    }
    else {
        /* The maximum over all y>=0 occurs at 0 or infinity:
           the sign of its derivative is sign((m+y)/s0^2-(y+zbeta*b)/sd^2).
           At infinity only the upper tail with threshold k remains.
           Bound the y=0 branch on intervals for b, or by the centered
           folded-normal quantile when p<.5 (sigma monotonicity then need
           not hold). This handles an unlisted Cartesian family safely. */
        if (p>=.5) {
            lower=(0::255)*(rho/256)
            zero=jw_mcp_foldlo(zbeta*lower,sqrt(1:-(lower:+rho/256):^2),p)
        }
        else zero=s0*invnormal((1+p)/2)
        k=min((z,z*s0+zbeta*rho))-1e-10
        curves=((k,1e300) \ (min(zero),min(zero)))
    }
    curves=curves \ (guardcv-1e-10,guardcv-1e-10)
    return(min((cheap,jw_mcp_integral(curves,means,s0))))
}

/* Joint MD loadings use V in both scopes. Only the selection caps differ.
   Estimate all requested targets with a single solve per candidate union. */
void jw_mcp_actions(real matrix V, real matrix c, real rowvector v,
    real matrix masks, real rowvector owners, string scalar scope,
    real scalar logn, real scalar guard, real matrix b,
    real matrix sd, real matrix caps)
{
    real rowvector ix,gain,part,ig,uv
    real matrix coef,S
    real scalar j,g
    b=J(rows(masks),cols(c),0)
    sd=J(rows(masks),cols(c),1)
    caps=J(rows(masks),cols(c),0)
    for (j=1;j<=rows(masks);j++) {
        ix=selectindex(masks[j,.])
        if (!cols(ix)) continue
        S=V[ix,ix]
        coef=cholsolve(S,c[ix,.])
        if (hasmissing(coef)) jw_core_error("singular candidate covariance in McCloskey profile")
        gain=colsum(coef:*(S*coef)):/v
        if (hasmissing(gain) | min(gain)<-1e-10 | max(gain)>=1)
            jw_core_error("degenerate candidate variance in McCloskey profile")
        b[j,.]=sqrt(colmax((J(1,cols(c),0)\gain)))
        sd[j,.]=sqrt(1:-b[j,.]:^2)
        if (scope=="joint") caps[j,.]=sqrt(cols(ix)*logn)*b[j,.]
        else {
            uv=uniqrows(sort(owners[ix]',1))'
            for (g=1;g<=cols(uv);g++) {
                ig=selectindex(owners[ix]:==uv[g])
                part=colsum(coef[ig,.]:*(S[ig,ig]*coef[ig,.])):/v
                caps[j,.]=caps[j,.]+sqrt(cols(ig)*logn):*sqrt(colmax((J(1,cols(c),0)\part)))
            }
            caps[j,.]=colmin((caps[j,.]\guard*b[j,.]))
        }
    }
    caps=caps:*(1+1e-11):+1e-11
}

struct jw_mc_result scalar jw_infer_mc_profile(struct jw_model scalar M,
    struct jw_family scalar C, struct jw_fit scalar F, real matrix joint,
    string scalar scope, real scalar maxjoint, real scalar level,
    real scalar gamma, real scalar mcerr, real rowvector outputs,
    real scalar progress, real scalar steps, real scalar radius, real scalar pmax,
    real scalar mcbeta, string scalar fallback, real scalar reps, real scalar cells,
    real scalar tol, real scalar pilotreps, real scalar seed, real matrix profile,
    real matrix grid, | real scalar refine, real scalar bandcells, real scalar bandreps)
{
    struct jw_mc_result scalar R
    struct jw_grid_result scalar BB
    real matrix masks,V,c,cc,acovs,b,sd,caps,targets
    real rowvector active,v,failed,pgrid,owners,ix,oktargets,band
    real colvector risks,h,localrows,quantiles
    real scalar o,j,g,i,d,nout,ns,maxpen,listed,B,alpha,target,tau,basecv,zbeta,
        rho,s0,guard,bias,move,chosen,width,nq,route,risk,cv,next,bandseed
    if (args()<24) refine=0
    if (args()<25) bandcells=256
    if (args()<26) bandreps=9999
    alpha=1-level/100
    target=alpha-mcerr/2
    if (missing(target) | gamma<=0 | mcerr<=0 | gamma+mcerr>=alpha |
        mcbeta<=0 | mcbeta>=1 | steps<2 | radius<0)
        jw_core_error("invalid McCloskey profile controls")
    if (scope!="joint" & scope!="cohort") jw_core_error("invalid WBIC selection scope")
    if (scope=="cohort" & rows(joint)) jw_core_error("jointmatrix() requires joint WBIC")
    if (fallback!="bb" & fallback!="error") jw_core_error("invalid McCloskey fallback")
    if (!pmax) pmax=1-target/100
    if (pmax<1-target | pmax>=1) jw_core_error("mcpmax() is outside the admissible probability range")
    pgrid=1:-target:+((0..(steps-1))/(steps-1))*(pmax-(1-target))
    nout=rows(F.A)
    if (!cols(outputs)) outputs=1..nout
    if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>nout | any(outputs:!=floor(outputs)))
        jw_core_error("invalid McCloskey output indices")
    outputs=uniqrows(sort(outputs',1))'
    R.table=J(nout,5,.)
    R.diagnostics=J(nout,16,.)
    R.calibration=J(steps,2+3*cols(outputs),.)
    R.calibration[.,1..2]=((1::steps),pgrid')
    R.bb_diagnostics=J(0,0,.)
    R.cis_diagnostics=R.cis_points=grid=J(0,0,.)
    profile=J(nout,11,.)
    for (o=1;o<=nout;o++) {
        R.table[o,1..2]=(F.A[o,.]*M.z,sqrt(jw_fixed_variance(M.V,F.A[o,.])))
        R.diagnostics[o,(1,14)]=(o,1)
        profile[o,(1,8,9)]=(o,1,9)
    }
    listed=1
    ns=1
    maxpen=0
    masks=J(0,M.k,.)
    active=J(1,0,.)
    if (M.k) {
        if (scope=="joint") {
            masks=jw_joint_masks(M,C,joint,maxjoint)
            ns=rows(masks)
            maxpen=max(rowsum(masks))*ln(M.nunit)
            active=selectindex(colsum(masks):>0)
        }
        else {
            active=selectindex(colsum(C.masks):>0)
            for (g=1;g<=rows(M.groups);g++) {
                localrows=vec(selectindex(C.owners:==M.groups[g]))
                ns=ns*rows(localrows)
                maxpen=maxpen+max(rowsum(C.masks[localrows,.]))*ln(M.nunit)
            }
            listed=(ns<=maxjoint & !missing(ns))
            if (listed) masks=jw_joint_masks(M,C,J(0,M.k,.),maxjoint)
        }
    }
    d=cols(active)
    R.domain_radius=R.exterior_radius=0
    if (!d) {
        cv=invnormal(1-alpha/2)
        for (i=1;i<=cols(outputs);i++) {
            o=outputs[i]
            R.table[o,3..4]=R.table[o,1]:+(-1,1)*cv*R.table[o,2]
            R.diagnostics[o,.]=(o,0,1,0,0,0,0,0,0,.,alpha,0,cv,2,1,0)
            profile[o,.]=(o,0,1,.,alpha,alpha,0,2,2,0,0)
        }
        R.calibration=J(0,0,.)
        return(R)
    }
    B=(radius==0 ? sqrt(maxpen)+sqrt(invchi2(d,.99)) : radius)
    R.domain_radius=B
    /* The profile covers all nuisance means with no exterior truncation. */
    R.exterior_radius=.
    tau=target/2
    basecv=invnormal(1-(target-tau)/2)
    zbeta=invnormal(1-mcbeta/2)
    V=M.V[M.p:+active,M.p:+active]
    h=M.z[M.p:+active]
    owners=M.hmap[active,1]'
    targets=F.targets[outputs,.]
    v=rowsum((targets*M.V):*targets)'
    c=M.V[M.p:+active,.]*targets'
    cc=cholsolve(V,(h,c))
    guard=(!hasmissing(cc) ? h'*cc[.,1]>B^2 : .)
    oktargets=J(1,0,.)
    acovs=J(1,cols(v),.)
    if (!hasmissing(cc)) {
        cc=cc[.,2..cols(cc)]
        for (i=1;i<=cols(v);i++) {
            if (v[i]<=0) continue
            acovs[i]=max((0,c[.,i]'*cc[.,i]/v[i]))
            if (!missing(acovs[i]) & 1-acovs[i]>1e-12) oktargets=(oktargets,i)
        }
    }
    b=sd=caps=J(0,cols(v),.)
    if (listed & cols(oktargets)) {
        jw_mcp_actions(V,c[.,oktargets],v[oktargets],masks[.,active],owners,
            scope,ln(M.nunit),B,b,sd,caps)
    }
    failed=J(1,0,.)
    nq=65536
    quantiles=invnormal(.5:+.5*(1::(nq-1))/nq)
    for (i=1;i<=cols(outputs);i++) {
        o=outputs[i]
        route=(listed ? 1 : 3)
        chosen=0
        risks=J(steps,1,.)
        R.diagnostics[o,1..9]=(o,d,ns,0,0,0,0,0,0)
        if (any(oktargets:==i) & all(c[.,i]:==0)) {
            /* With independent baseline noise all MD adjustments vanish.
               The guard and ordinary normal branches each have risk <=
               target even at the first p, allowing equality without a pad. */
            route=2
            chosen=1
            risks[1]=target
        }
        else if (any(oktargets:==i)) {
            ix=selectindex(oktargets:==i)
            j=ix[1]
            rho=sqrt(acovs[i])
            s0=sqrt(1-acovs[i])
            for (g=1;g<=steps;g++) {
                if (listed) risk=jw_mcp_upper(b[.,j],sd[.,j],caps[.,j],rho*quantiles,s0,rho,pgrid[g],zbeta,basecv,1,target)
                else risk=jw_mcp_upper(J(0,1,.),J(0,1,.),J(0,1,.),rho*quantiles,s0,rho,pgrid[g],zbeta,basecv,0,target)
                risks[g]=risk
                if (risk<=target) {
                    chosen=g
                    break
                }
            }
        }
        else route=4
        profile[o,10..11]=(0,0)
        if (refine & route==1) {
            /* Spend only on a narrower p than the retained deterministic
               certificate. Fresh stage seeds and an allocation for every
               possible level allow early stopping at an unresolved level. */
            next=(chosen ? chosen-1 : steps)
            for (g=next;g>=1;g--) {
                bandseed=mod(seed+4099*o+104729*g,2147483646)
                band=jw_mcp_band(V,c[.,i],v[i],masks[.,active],owners,scope,
                    ln(M.nunit),B,target,mcbeta,pgrid[g],bandseed,bandcells,bandreps,
                    mcerr/(2*steps*cols(outputs)))
                if (cols(band)!=3 | hasmissing(band) | band[1]<0 | band[1]>1 |
                    min(band[2..3])<0) jw_core_error("invalid McCloskey band bound")
                profile[o,10..11]=profile[o,10..11]+band[2..3]
                risks[g]=min((risks[g],band[1]))
                if (risks[g]>target) break
                chosen=g
                route=5
            }
        }
        R.calibration[.,1+3*i]=risks
        R.calibration[.,2+3*i]=J(steps,1,0)
        profile[o,(1,2,3,6,7,9)]=(o,d,ns,target,nq,route)
        if (!chosen) {
            failed=(failed,o)
            R.diagnostics[o,11]=risks[steps]
            profile[o,5]=risks[steps]
            continue
        }
        /* Only after choosing p evaluate the actual public joint-MD action. */
        if (max(abs(F.A[o,.]-jw_loading(M.V,F.targets[o,.],F.selected[o,.],M.k)))>
            1e-8*(1+max(abs(F.A[o,.]))))
            jw_core_error("McCloskey profile requires the public union joint-MD action")
        ix=M.p:+active
        move=(F.A[o,ix]-F.targets[o,ix])*h
        bias=sqrt(max((0,(F.A[o,ix]-F.targets[o,ix])*V*(F.A[o,ix]-F.targets[o,ix])')))
        if (guard) width=abs(move)+basecv*sqrt(v[i])
        else width=jw_mcr_foldquantile(abs(move)+zbeta*bias,R.table[o,2],pgrid[chosen])
        R.table[o,3..4]=R.table[o,1]:+(-1,1)*width
        R.diagnostics[o,10..16]=(pgrid[chosen],risks[chosen],0,width/R.table[o,2],0,1,guard)
        profile[o,(4,5,8)]=(pgrid[chosen],risks[chosen],0)
        if (progress) printf("McCloskey profile output %g: p=%g, global risk <= %g, route %g\n",o,pgrid[chosen],risks[chosen],route)
    }
    if (cols(failed)) {
        if (fallback=="error") jw_core_error("McCloskey profile did not certify a level; use mcfallback(bb) or expand the probability grid")
        BB=jw_infer_bb_grid(M,C,F,joint,scope,maxjoint,100*(1-target),gamma,mcerr/2,
            reps,cells,tol,failed,progress,pilotreps,seed)
        R.bb_diagnostics=BB.diagnostics
        grid=BB.grid
        for (i=1;i<=cols(failed);i++) {
            o=failed[i]
            R.table[o,.]=BB.table[o,.]
            R.diagnostics[o,13..15]=(BB.diagnostics[o,10],5,0)
            profile[o,8]=5
        }
    }
    return(R)
}

/* ===== jwselect_jobs.mata ===== */
/* jwselect 0.6.0: portable, resumable objective-by-start jobs.

   The 0.3.0 serialized-format identifiers are retained because the binary
   layouts are unchanged; jw_version() identifies the package build.

   Only numeric/string matrices are serialized with fputmatrix, not opaque
   compiled code or raw observations. Prepare and serial inference consume
   random numbers in the same order. Workers consume no RNG. The Stata
   interface protects immutable files with checksums and acquires job locks.
   A .done marker is created only after its result file has been closed.
   An interrupted unmarked result can be replaced; a marked corrupt result
   is an error. No projection cache is written to disk.
*/

struct jw_job_state {
    real rowvector token
    real rowvector settings
    real rowvector numerical
    real matrix index
    real matrix info
    real matrix outmap
    struct jw_reselect_result scalar result
}

string scalar jw_jobs_num(real scalar x)
{
    return(strtrim(strofreal(x,"%18.0f")))
}

string scalar jw_jobs_output_file(string scalar path, real scalar o)
{
    return(path+"/output_"+jw_jobs_num(o)+".bin")
}

string scalar jw_jobs_result_file(string scalar path, real scalar id)
{
    return(path+"/job_"+jw_jobs_num(id)+".bin")
}

void jw_jobs_write_context(real scalar fh, struct jw_rs_context scalar C)
{
    fputmatrix(fh,C.reportB)
    fputmatrix(fh,C.selectB)
    fputmatrix(fh,C.constant)
    fputmatrix(fh,C.component)
    fputmatrix(fh,C.groups)
    fputmatrix(fh,C.EH)
    fputmatrix(fh,C.baseline)
}

struct jw_rs_context scalar jw_jobs_read_context(real scalar fh)
{
    struct jw_rs_context scalar C
    C.reportB = fgetmatrix(fh)
    C.selectB = fgetmatrix(fh)
    C.constant = fgetmatrix(fh)
    C.component = fgetmatrix(fh)
    C.groups = fgetmatrix(fh)
    C.EH = fgetmatrix(fh)
    C.baseline = fgetmatrix(fh)
    C.projection = J(0,0,.)
    if (rows(C.reportB)!=rows(C.selectB) |
        cols(C.reportB)!=cols(C.selectB) |
        rows(C.constant)!=rows(C.reportB) | cols(C.constant)!=1 |
        rows(C.component)!=rows(C.reportB) | cols(C.component)!=1 |
        cols(C.EH)!=cols(C.reportB) |
        rows(C.baseline)!=rows(C.EH) | cols(C.baseline)!=1 |
        hasmissing(C.reportB) | hasmissing(C.selectB) |
        hasmissing(C.constant) | hasmissing(C.EH) |
        hasmissing(C.baseline)) {
        jw_core_error("invalid or incompatible frozen search context")
    }
    return(C)
}

void jw_jobs_write_output(string scalar filename, struct jw_rs_output scalar P)
{
    real scalar fh
    fh = fopen(filename,"w")
    fputmatrix(fh,"jwselect-output-0.3.0")
    fputmatrix(fh,(P.output,P.estimate))
    fputmatrix(fh,P.active)
    fputmatrix(fh,P.hcols)
    fputmatrix(fh,P.fullmu)
    fputmatrix(fh,P.h)
    fputmatrix(fh,P.B)
    fputmatrix(fh,P.starts)
    jw_jobs_write_context(fh,P.C)
    jw_jobs_write_context(fh,P.CV)
    fputmatrix(fh,"complete")
    fclose(fh)
}

struct jw_rs_output scalar jw_jobs_read_output(string scalar path,
                                               real scalar o)
{
    struct jw_rs_output scalar P
    real scalar fh, d
    real rowvector meta
    string scalar magic
    fh = fopen(jw_jobs_output_file(path,o),"r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-output-0.3.0") jw_core_error("incompatible frozen output format")
    meta = fgetmatrix(fh)
    if (cols(meta)!=2 | rows(meta)!=1) jw_core_error("invalid frozen output header")
    P.output = meta[1]
    P.estimate = meta[2]
    P.active = fgetmatrix(fh)
    P.hcols = fgetmatrix(fh)
    P.fullmu = fgetmatrix(fh)
    P.h = fgetmatrix(fh)
    P.B = fgetmatrix(fh)
    P.starts = fgetmatrix(fh)
    P.C = jw_jobs_read_context(fh)
    P.CV = jw_jobs_read_context(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    d = cols(P.active)
    if (magic!="complete" | P.output!=o | d<1 |
        rows(P.h)!=d | cols(P.h)!=1 | rows(P.B)!=d | cols(P.B)!=d |
        rows(P.starts)!=d | cols(P.starts)<1 | missing(P.estimate) |
        hasmissing(P.h) | hasmissing(P.B) | hasmissing(P.starts) |
        cols(P.C.EH)!=d | cols(P.CV.EH)!=d) {
        jw_core_error("invalid or incomplete frozen output")
    }
    return(P)
}

void jw_jobs_write_state(string scalar path, struct jw_job_state scalar JS)
{
    real scalar fh
    fh = fopen(path+"/jw_state.bin","w")
    fputmatrix(fh,"jwselect-jobs-0.3.0")
    fputmatrix(fh,JS.token)
    fputmatrix(fh,JS.settings)
    fputmatrix(fh,JS.numerical)
    fputmatrix(fh,JS.index)
    fputmatrix(fh,JS.info)
    fputmatrix(fh,JS.outmap)
    fputmatrix(fh,JS.result.table)
    fputmatrix(fh,JS.result.diagnostics)
    fputmatrix(fh,JS.result.extrema)
    fputmatrix(fh,"complete")
    fclose(fh)
}

struct jw_job_state scalar jw_jobs_read_state(string scalar path)
{
    struct jw_job_state scalar JS
    real scalar fh
    string scalar magic
    fh = fopen(path+"/jw_state.bin","r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-jobs-0.3.0") jw_core_error("incompatible jobs state format")
    JS.token = fgetmatrix(fh)
    JS.settings = fgetmatrix(fh)
    JS.numerical = fgetmatrix(fh)
    JS.index = fgetmatrix(fh)
    JS.info = fgetmatrix(fh)
    JS.outmap = fgetmatrix(fh)
    JS.result.table = fgetmatrix(fh)
    JS.result.diagnostics = fgetmatrix(fh)
    JS.result.extrema = fgetmatrix(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    if (magic!="complete" | cols(JS.settings)!=9 | rows(JS.settings)!=1 |
        cols(JS.numerical)!=4 | rows(JS.numerical)!=1 |
        cols(JS.index)!=4 | cols(JS.info)!=3 |
        cols(JS.result.table)!=5 | cols(JS.result.diagnostics)!=13 |
        rows(JS.result.diagnostics)!=4*rows(JS.result.table) |
        rows(JS.result.extrema)!=4*rows(JS.result.table) |
        hasmissing(JS.token) | hasmissing(JS.settings) |
        hasmissing(JS.numerical) | hasmissing(JS.index)) {
        jw_core_error("invalid or incomplete jobs state")
    }
    return(JS)
}

real matrix jw_jobs_index(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.index)
}

real matrix jw_jobs_output_map(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.outmap)
}

real matrix jw_jobs_output_info(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.info)
}

real rowvector jw_jobs_settings(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.settings)
}

string colvector jw_jobs_files(string scalar path)
{
    struct jw_job_state scalar JS
    string colvector files
    real scalar o, n, i
    JS = jw_jobs_read_state(path)
    n = sum(JS.info[,2]:>0)
    files = J(n+1,1,"")
    files[1] = "jw_state.bin"
    i = 1
    for (o=1; o<=rows(JS.info); o++) {
        if (JS.info[o,2]>0) {
            i++
            files[i] = "output_"+jw_jobs_num(JS.info[o,1])+".bin"
        }
    }
    return(files)
}

/* Preparation contains no optimization, and never evaluates an omitted
   output. Empty outputs means all. The caller creates a new directory. */
void jw_jobs_prepare(struct jw_model scalar M, struct jw_fit scalar F,
                    string scalar path, real scalar level,
                    real scalar gamma, real scalar mcerr, real scalar reps,
                    real scalar optstarts, real scalar optmaxeval,
                    real scalar opttol, real scalar optreps,
                    real rowvector outputs, real scalar cachemb,
                    real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    real matrix S, E, EV, index
    real scalar o, j, s, d, initialized, id, nmax, fingerprint
    if (!direxists(path)) jw_core_error("jobs directory does not exist")
    if (fileexists(path+"/jw_state.bin")) {
        jw_core_error("prepared state already exists; use a new directory")
    }
    JS.numerical = jw_rs_validate(M,F,level,gamma,mcerr,reps,optstarts,
                        optmaxeval,opttol,optreps,outputs,cachemb)
    JS.settings = (level,gamma,mcerr,reps,optstarts,optmaxeval,opttol,optreps,cachemb)
    JS.result = jw_rs_result_init(M,F,level)
    JS.outmap = F.outmap
    JS.info = J(rows(F.A),3,0)
    nmax = rows(F.A)*4*optstarts
    index = J(nmax,4,.)
    id = 0
    initialized = 0
    fingerprint = 0
    S = J(0,0,.)
    E = J(0,0,.)
    for (o=1; o<=rows(F.A); o++) {
        JS.info[o,1] = o
        if (cols(outputs)>0) {
            if (!any(outputs:==o)) continue
        }
        P = jw_rs_output_plan(M,F,o,gamma)
        d = cols(P.active)
        if (d==0) {
            jw_rs_zero_output(JS.result,P,jw_fixed_variance(M.V,F.A[o,.]),level)
            continue
        }
        if (JS.numerical[3]<1 | JS.numerical[4]>reps) {
            jw_core_error("reps() too small for finite Monte Carlo quantile bounds; increase reps()")
        }
        if (!initialized) {
            S = jw_infer_sqrt(M.V)
            E = rnormal(optreps,rows(M.z),0,1)*S'
            initialized = 1
            fingerprint = sum(E)
        }
        P.C = jw_rs_prepare_opts(M,F,o,E,P.hcols,0)
        P.starts = jw_rs_starts(P.h,P.B,optstarts)
        EV = rnormal(reps,rows(M.z),0,1)*S'
        P.CV = jw_rs_prepare_opts(M,F,o,EV,P.hcols,0)
        fingerprint = fingerprint+o*sum(EV)+sum(P.starts)
        jw_jobs_write_output(jw_jobs_output_file(path,o),P)
        JS.info[o,.] = (o,d,4*optstarts)
        for (j=1; j<=4; j++) {
            for (s=1; s<=optstarts; s++) {
                id++
                index[id,.] = (id,o,j,s)
            }
        }
        if (progress) {
            printf("Prepared output %g: %g objectives x %g starts, dimension %g\n",
                    o,4,optstarts,d)
            displayflush()
        }
    }
    JS.index = J(0,4,.)
    if (id>0) JS.index = index[|1,1\id,4|]
    /* No random draw is consumed for provenance. Together with the front-
       end's immutable-file checksums this catches copied/mismatched jobs. */
    JS.token = (300,JS.settings[1..8],rows(M.z),M.k,rows(F.A),
                sum(M.z),sum(M.V),sum(F.A),fingerprint)
    jw_jobs_write_state(path,JS)
}

void jw_jobs_save_result(string scalar path, real rowvector token,
                         real rowvector job,
                         struct jw_rs_search_result scalar R)
{
    real scalar fh
    string scalar filename
    filename = jw_jobs_result_file(path,job[1])
    if (fileexists(filename+".done")) {
        jw_core_error("refusing to overwrite a completed search job")
    }
    /* fopen(...,"w") requires a new file. Only an unmarked partial file
       may be removed here, after the front-end acquires this job's lock. */
    unlink(filename)
    fh = fopen(filename,"w")
    fputmatrix(fh,"jwselect-search-0.3.0")
    fputmatrix(fh,token)
    fputmatrix(fh,job)
    fputmatrix(fh,R.x)
    fputmatrix(fh,(R.value,R.evaluations,R.starts,R.converged,R.exhausted))
    fputmatrix(fh,"complete")
    fclose(fh)
    /* The marker is deliberately written last. */
    fh = fopen(filename+".done","w")
    fputmatrix(fh,token)
    fputmatrix(fh,job)
    fclose(fh)
}

struct jw_rs_search_result scalar jw_jobs_load_result(string scalar path,
                    real rowvector token, real rowvector job,
                    real scalar d, real scalar budget)
{
    struct jw_rs_search_result scalar R
    real scalar fh
    real rowvector savedtoken, savedjob, meta
    string scalar filename, magic
    filename = jw_jobs_result_file(path,job[1])
    if (!fileexists(filename+".done")) jw_core_error("required search job is unfinished")
    fh = fopen(filename+".done","r")
    savedtoken = fgetmatrix(fh)
    savedjob = fgetmatrix(fh)
    fclose(fh)
    if (rows(savedtoken)!=1 | cols(savedtoken)!=cols(token) |
        rows(savedjob)!=1 | cols(savedjob)!=4) jw_core_error("invalid job completion marker")
    if (any(savedtoken:!=token) | any(savedjob:!=job)) {
        jw_core_error("job completion marker belongs to a different preparation or job")
    }
    fh = fopen(filename,"r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-search-0.3.0") jw_core_error("incompatible search result format")
    savedtoken = fgetmatrix(fh)
    savedjob = fgetmatrix(fh)
    R.x = fgetmatrix(fh)
    meta = fgetmatrix(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    if (magic!="complete" | rows(savedtoken)!=1 | cols(savedtoken)!=cols(token) |
        rows(savedjob)!=1 | cols(savedjob)!=4 | rows(R.x)!=d | cols(R.x)!=1 |
        rows(meta)!=1 | cols(meta)!=5 | hasmissing(R.x) | hasmissing(meta)) {
        jw_core_error("invalid or incomplete completed search result")
    }
    if (any(savedtoken:!=token) | any(savedjob:!=job) |
        sqrt(sum(R.x:^2))>1+1e-10 | meta[2]<1 | meta[2]>budget |
        meta[2]!=floor(meta[2]) | meta[3]!=1 |
        (meta[4]!=0 & meta[4]!=1) | (meta[5]!=0 & meta[5]!=1) |
        meta[4]+meta[5]!=1 | (meta[5]==1 & meta[2]!=budget)) {
        jw_core_error("completed search result is incompatible with its job")
    }
    R.value = meta[1]
    R.evaluations = meta[2]
    R.starts = meta[3]
    R.converged = meta[4]
    R.exhausted = meta[5]
    return(R)
}

void jw_jobs_check_filter(real rowvector filter, real scalar upper,
                           string scalar label)
{
    if (cols(filter)>0) {
        if (hasmissing(filter) | min(filter)<1 | max(filter)>upper |
            any(filter:!=floor(filter))) {
            jw_core_error("invalid "+label+" filter")
        }
    }
}

/* Filters are rowvectors; empty means all. Modulo assignment always uses
   the original global job id, independently of filters. Root owns locks. */
void jw_jobs_worker_ids(string scalar path, real scalar worker,
                    real scalar nworkers, real rowvector outputfilter,
                    real rowvector objectivefilter, real rowvector startfilter,
                    real rowvector jobfilter, real scalar cachemb,
                    real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    struct jw_rs_search_result scalar R
    real scalar i, o, objective, start, lastoutput, budget
    real rowvector job
    JS = jw_jobs_read_state(path)
    if (worker<1 | worker>nworkers | nworkers<1 |
        worker!=floor(worker) | nworkers!=floor(nworkers) |
        missing(worker) | missing(nworkers) | missing(cachemb) | cachemb<0) {
        jw_core_error("invalid worker index, worker count, or cache budget")
    }
    jw_jobs_check_filter(outputfilter,rows(JS.info),"output")
    jw_jobs_check_filter(objectivefilter,4,"objective")
    jw_jobs_check_filter(startfilter,JS.settings[5],"starting-point")
    jw_jobs_check_filter(jobfilter,rows(JS.index),"job")
    lastoutput = 0
    for (i=1; i<=rows(JS.index); i++) {
        job = JS.index[i,.]
        if (mod(job[1]-1,nworkers)!=worker-1) continue
        if (cols(jobfilter)>0) {
            if (!any(jobfilter:==job[1])) continue
        }
        o = job[2]
        objective = job[3]
        start = job[4]
        if (cols(outputfilter)>0) {
            if (!any(outputfilter:==o)) continue
        }
        if (cols(objectivefilter)>0) {
            if (!any(objectivefilter:==objective)) continue
        }
        if (cols(startfilter)>0) {
            if (!any(startfilter:==start)) continue
        }
        budget = jw_rs_start_budget(JS.settings[6],JS.settings[5],start)
        if (fileexists(jw_jobs_result_file(path,job[1])+".done")) {
            R = jw_jobs_load_result(path,JS.token,job,JS.info[o,2],budget)
            if (progress) printf("Search job %g already complete; retained.\n",job[1])
            continue
        }
        if (o!=lastoutput) {
            P = jw_jobs_read_output(path,o)
            jw_rs_cache(P.C,cachemb)
            lastoutput = o
        }
        if (progress) {
            printf("Search job %g/%g: output %g, objective %g, start %g, budget %g\n",
                    job[1],rows(JS.index),o,objective,start,budget)
            displayflush()
        }
        R = jw_rs_search_start(P.C,P.h,P.B,P.starts[,start],JS.numerical[1],
                P.estimate,objective,budget,JS.settings[7],progress)
        jw_jobs_save_result(path,JS.token,job,R)
        if (progress) {
            printf("Search job %g complete: %g evaluations, converged=%g, budget_exhausted=%g\n",
                    job[1],R.evaluations,R.converged,R.exhausted)
            displayflush()
        }
    }
}

/* Compatibility entry point: no explicit job-id filter. */
void jw_jobs_worker(string scalar path, real scalar worker,
                    real scalar nworkers, real rowvector outputfilter,
                    real rowvector objectivefilter, real rowvector startfilter,
                    real scalar cachemb, real scalar progress)
{
    jw_jobs_worker_ids(path,worker,nworkers,outputfilter,objectivefilter,
                       startfilter,J(1,0,.),cachemb,progress)
}

struct jw_reselect_result scalar jw_jobs_collect(string scalar path,
                                                 real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    struct jw_rs_search_result scalar R, S
    real matrix located
    real rowvector searchvalues, job
    real scalar o, objective, start, id, budget, d
    JS = jw_jobs_read_state(path)
    id = 0
    /* Inspect every required completion before doing costly validation. */
    for (o=1; o<=rows(JS.index); o++) {
        if (!fileexists(jw_jobs_result_file(path,JS.index[o,1])+".done")) {
            jw_core_error("cannot collect: one or more required search jobs are unfinished")
        }
    }
    for (o=1; o<=rows(JS.info); o++) {
        d = JS.info[o,2]
        if (d==0) continue
        P = jw_jobs_read_output(path,o)
        located = J(d,4,.)
        searchvalues = J(1,4,.)
        for (objective=1; objective<=4; objective++) {
            R.value = .
            R.x = J(d,1,.)
            R.evaluations = 0
            R.starts = 0
            R.converged = 0
            R.exhausted = 0
            for (start=1; start<=JS.settings[5]; start++) {
                id++
                job = JS.index[id,.]
                if (any(job:!=(id,o,objective,start))) jw_core_error("invalid search job order")
                budget = jw_rs_start_budget(JS.settings[6],JS.settings[5],start)
                S = jw_jobs_load_result(path,JS.token,job,d,budget)
                jw_rs_merge(R,S,objective)
            }
            located[,objective] = P.h+P.B*R.x
            searchvalues[objective] = R.value
            jw_rs_record_search(JS.result,P,R,objective,JS.settings[8],JS.settings[4])
        }
        /* CV is deliberately uncached; only four validation points. */
        jw_rs_finish_output(JS.result,P,located,searchvalues,JS.numerical,
                            JS.settings[2],JS.settings[3],progress)
    }
    return(JS.result)
}

/* ===== jwselect_bridge.mata ===== */
/* Stata transport only. Numerical modules never depend on st_* state. */
real matrix jw_optional_matrix(string scalar localname)
{
    string scalar name
    name=st_local(localname)
    if (name=="") return(J(0,0,.))
    return(st_matrix(name))
}

void jw_put(string scalar localname, real matrix value)
{
    if (rows(value)>0 & cols(value)>0) st_matrix(st_local(localname),value)
}

void jw_names(string scalar localname, string colvector names, real scalar rownames)
{
    string matrix stripe
    stripe=J(rows(names),1,"") , names
    if (rownames) st_matrixrowstripe(st_local(localname),stripe)
    else st_matrixcolstripe(st_local(localname),stripe)
}

void jw_stata()
{
    struct jw_model scalar M
    struct jw_family scalar C
    struct jw_fit scalar F
    struct jw_reselect_result scalar RI
    struct jw_mc_result scalar MC
    struct jw_pw_result scalar PW
    real matrix table, Vout, Vcells, joint, pw, blocks,profile,profilegrid
    real rowvector outputs
    real colvector estimate
    real scalar i, j, G, m, indexblocks
    string colvector names, cnames, hnames, znames, familycols, diagcols
    string rowvector scorecols
    string scalar selector, focus, inference

    selector=st_local("select")
    focus=st_local("focus")
    inference=st_local("inference")
    pw=J(1,0,.)
    if (st_local("postwindow")!="") pw=strtoreal(tokens(st_local("postwindow")))
    indexblocks = st_local("restriction_scale")=="conditional_index"
    if (indexblocks) {
        M=jw_design_index(st_matrix(st_local("rawb")),st_matrix(st_local("rawV")),
            st_matrix(st_local("rawcells")),jw_optional_matrix("rawh"),
            jw_optional_matrix("rawVh"),jw_optional_matrix("rawVdh"),
            jw_optional_matrix("rawhmap"),strtoreal(st_local("nunits")),
            strtoreal(st_local("gap")),strtoreal(st_local("anticipation")),
            st_local("baseline"),pw,st_local("aggregate"),jw_optional_matrix("aggweights"))
        C=jw_index_candidates(M,selector,st_local("candidates"),jw_optional_matrix("prematrix"),
            jw_optional_matrix("candmatrix"),strtoreal(st_local("maxcand")))
    }
    else {
        M=jw_design(st_matrix(st_local("rawb")),st_matrix(st_local("rawV")),
            st_matrix(st_local("rawcells")),strtoreal(st_local("nunits")),
            strtoreal(st_local("gap")),strtoreal(st_local("anticipation")),
            st_local("baseline"),pw,st_local("aggregate"),jw_optional_matrix("aggweights"))
        C=jw_candidates(M,selector,st_local("candidates"),jw_optional_matrix("prematrix"),
            jw_optional_matrix("candmatrix"),strtoreal(st_local("maxcand")))
    }
    joint=jw_optional_matrix("jointmatrix")
    if (indexblocks & rows(joint)>0) jw_block_masks_check(M.hmap,joint)
    F=jw_run_scope(M,C,joint,selector,focus,strtoreal(st_local("maxjoint")),st_local("wbicscope"))
    outputs=J(1,0,.)
    if (st_local("inferoutputs")!="") outputs=strtoreal(tokens(st_local("inferoutputs")))
    if (cols(outputs)>0) {
        if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>rows(F.A) |
            any(outputs:!=floor(outputs))) jw_core_error("invalid inferoutputs() row index")
    }
    if (st_local("prepare")!="") {
        jw_jobs_prepare(M,F,st_local("prepare"),strtoreal(st_local("level")),
            strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),
            strtoreal(st_local("reps")),strtoreal(st_local("optstarts")),
            strtoreal(st_local("optmaxeval")),strtoreal(st_local("opttol")),
            strtoreal(st_local("optreps")),outputs,strtoreal(st_local("cachemb")),
            st_local("progress")!="")
    }
    estimate=F.A*M.z
    Vout=F.A*M.V*F.A'
    Vout=(Vout+Vout')/2
    for (i=1; i<=rows(F.A); i++) Vout[i,i]=jw_fixed_variance(M.V,F.A[i,.])
    Vcells=F.cellA*M.V*F.cellA'
    Vcells=(Vcells+Vcells')/2
    for (i=1; i<=M.p; i++) Vcells[i,i]=jw_fixed_variance(M.V,F.cellA[i,.])
    if (inference=="none") {
        table=J(rows(F.A),5,.)
        table[.,1]=estimate
        for (i=1; i<=rows(F.A); i++) table[i,2]=sqrt(jw_fixed_variance(M.V,F.A[i,.]))
    }
    else if (inference=="wbic") {
        if (st_local("wbicmethod")=="bootstrap") {
            PW=jw_infer_wbic_pointwise(M,C,F,joint,st_local("wbicscope"),
                strtoreal(st_local("maxjoint")),strtoreal(st_local("level")),
                strtoreal(st_local("reps")),outputs,st_local("wbicengine"))
            table=PW.table
            jw_put("outpwdiag",PW.diagnostics)
            jw_names("outpwdiag",("output"\"nactive"\"candidates"\"reps"\
                "q_lo"\"q_hi"\"adjustment"\"status"),0)
        }
        else if (st_local("wbicmethod")=="mc") {
            MC=jw_infer_mc_profile(M,C,F,joint,st_local("wbicscope"),
                strtoreal(st_local("maxjoint")),strtoreal(st_local("level")),
                strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),outputs,
                st_local("progress")!="",strtoreal(st_local("mcsteps")),
                strtoreal(st_local("mcradius")),strtoreal(st_local("mcpmax")),
                strtoreal(st_local("mcbeta")),st_local("mcfallback"),
                strtoreal(st_local("reps")),strtoreal(st_local("bbcells")),
                strtoreal(st_local("bbtol")),strtoreal(st_local("isreps")),
                strtoreal(st_local("seed")),profile,profilegrid,
                st_local("mcrefine")!="",strtoreal(st_local("mccells")),
                strtoreal(st_local("mcvalidate")))
            jw_put("outprofile",profile)
            jw_names("outprofile",("output"\ "nactive"\ "candidates"\ "pstar"\
                "risk_upper"\ "risk_target"\ "quadrature"\ "status"\ "route"\
                "band_cells"\ "band_draws"),0)
            if (rows(profilegrid)>0) {
                jw_put("outgrid",profilegrid)
                jw_names("outgrid",("output"\ "nactive"\ "candidates"\ "cells"\ "vertices"\
                    "draws"\ "risk_upper"\ "critical_lo"\ "critical_hi"\ "gap"\ "status"\ "route"),0)
            }
            table=MC.table
            if (rows(MC.bb_diagnostics)>0) {
                jw_put("outbbdiag",MC.bb_diagnostics)
                jw_names("outbbdiag",("output"\"nactive"\"candidates"\"nodes"\
                    "draws"\"target_tail"\"critical_lo"\"critical_hi"\"gap"\
                    "critical"\"status"\"certified"\"accuracy_certified"\
                    "cap_used"\"max_radius"\"fixed_point_cv"),0)
            }
            jw_put("outmcdiag",MC.diagnostics)
            jw_names("outmcdiag",("output"\"nactive"\"candidates"\"inner_nodes"\
                "inner_reps"\"search_nodes"\"outer_reps"\"validation_nodes"\
                "validation_reps"\"pstar"\"risk_upper"\"shift_radius"\
                "critical"\"status"\"certified"\"rule_fallback"),0)
            if (rows(MC.calibration)>0) jw_put("outmccal",MC.calibration)
            st_local("_mc_domain_radius",strofreal(MC.domain_radius,"%21.17g"))
            st_local("_mc_exterior_radius",strofreal(MC.exterior_radius,"%21.17g"))
        }
        else jw_core_error("unknown adjusted WBIC method")
    }
    else if (selector=="focused" & (inference=="auto" | inference=="reselect")) {
        RI=jw_infer_reselect_opts(M,F,strtoreal(st_local("level")),
            strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),
            strtoreal(st_local("reps")),strtoreal(st_local("optstarts")),
            strtoreal(st_local("optmaxeval")),strtoreal(st_local("opttol")),
            strtoreal(st_local("optreps")),outputs,strtoreal(st_local("cachemb")),
            st_local("progress")!="")
        table=RI.table
        jw_put("outdiagnostics",RI.diagnostics)
        jw_put("outextrema",RI.extrema)
    }
    else table=jw_infer(M,F,selector,strtoreal(st_local("level")),
        strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),strtoreal(st_local("reps")))
    G=rows(M.groups)
    m=rows(M.z)
    names=J(rows(F.A),1,"")
    for (i=1; i<=rows(F.A); i++) {
        if (F.outmap[i,1]==1 | F.outmap[i,1]==5) names[i]="cohort"+strtrim(strofreal(i))
        else if (F.outmap[i,1]==2) names[i]="agg_cohort"
        else if (F.outmap[i,1]==3) names[i]="agg_focused"
        else names[i]="aggregate"
    }
    cnames=J(M.p,1,"")
    for (i=1; i<=M.p; i++) cnames[i]="cell"+strtrim(strofreal(i))
    hnames=J(M.k,1,"")
    for (i=1; i<=M.k; i++) hnames[i]="h"+strtrim(strofreal(i))
    if (selector=="focused" & (inference=="auto" | inference=="reselect")) {
        diagcols=("output"\"objective"\"nactive"\"search_reps"\"validation_reps"\
            "evaluations"\"starts_used"\"converged_starts"\"budget_exhausted"\
            "search_value"\"validated_value"\"validation_gap"\"radius_fraction")
        jw_names("outdiagnostics",diagcols,0)
        if (rows(RI.extrema)>0) jw_names("outextrema",("output"\"objective"\hnames),0)
    }
    znames=cnames\hnames
    jw_put("outb",estimate')
    jw_put("outV",Vout)
    jw_put("outtable",table)
    jw_put("outmap",F.outmap)
    jw_put("outselection",F.selected)
    jw_put("outbases",(M.groups,M.bases))
    jw_put("outpost",M.postmap)
    jw_put("outh",M.hmap)
    blocks=jw_block_map(M.hmap)
    if (indexblocks & rows(blocks)>0) {
        jw_put("outblockmap",blocks)
        jw_names("outblockmap",("cohort"\"time"\"rank"),0)
    }
    jw_put("outbcells",(F.cellA*M.z)')
    jw_put("outVcells",Vcells)
    jw_put("outcellsel",F.cell_selected)
    jw_put("outweights",M.w)
    jw_put("outscores",F.scores)
    jw_put("outfamily",(C.owners,C.ids,C.masks))
    if ((selector=="focused" & focus!="cohort") | (selector=="wbic" & st_local("wbicscope")=="joint")) {
        joint=jw_joint_masks(M,C,joint,strtoreal(st_local("maxjoint")))
        jw_put("outjoint",joint)
        if (M.k>0) jw_names("outjoint",hnames,0)
    }
    jw_put("outz",M.z')
    jw_put("outzV",M.V)
    jw_put("outA",F.A)
    jw_put("outT",F.targets)
    jw_names("outb",names,0)
    jw_names("outV",names,0)
    jw_names("outV",names,1)
    jw_names("outtable",names,1)
    if (inference=="wbic" & st_local("wbicmethod")=="bootstrap") {
        jw_names("outpwdiag",names,1)
    }
    jw_names("outtable",("estimate"\"se_fixed"\"ll"\"ul"\"pvalue"),0)
    jw_names("outmap",names,1)
    jw_names("outmap",("kind"\"cohort"),0)
    jw_names("outbases",("cohort"\"baseline"),0)
    jw_names("outpost",("cohort"\"time"),0)
    jw_names("outpost",cnames,1)
    if (M.k>0) {
        jw_names("outh",("cohort"\"time"),0)
        jw_names("outh",hnames,1)
        jw_names("outselection",names,1)
        jw_names("outselection",hnames,0)
        jw_names("outcellsel",cnames,1)
        jw_names("outcellsel",hnames,0)
    }
    jw_names("outbcells",cnames,0)
    jw_names("outVcells",cnames,0)
    jw_names("outVcells",cnames,1)
    jw_names("outweights",cnames,0)
    scorecols=("scope","cohort","candidate","nrestr","Wald","WBIC","estimate","se_fixed","adjustment","gain","FMSC")
    jw_names("outscores",scorecols',0)
    familycols=("cohort"\"candidate"\hnames)
    jw_names("outfamily",familycols,0)
    jw_names("outz",znames,0)
    jw_names("outzV",znames,0)
    jw_names("outzV",znames,1)
    jw_names("outA",names,1)
    jw_names("outA",znames,0)
    jw_names("outT",names,1)
    jw_names("outT",znames,0)
    st_local("_G",strofreal(G))
    st_local("_K",strofreal(M.k))
    st_local("_K_blocks",strofreal(rows(blocks)))
    st_local("_P",strofreal(M.p))
}

/* r-class transport for completed divided searches; no e() mutation. */
void jw_worker_index_stata()
{
    real matrix index
    index=jw_jobs_index(st_local("directory"))
    st_local("njobs",strtrim(strofreal(rows(index))))
    if (rows(index)>0) st_matrix(st_local("jobs"),index)
}

void jw_collect_stata()
{
    struct jw_reselect_result scalar R
    string scalar path
    real matrix index
    path=st_local("directory")
    R=jw_jobs_collect(path,st_local("progress")!="")
    st_matrix(st_local("table"),R.table)
    st_matrix(st_local("diagnostics"),R.diagnostics)
    st_matrix(st_local("extrema"),R.extrema)
    st_matrix(st_local("output_map"),jw_jobs_output_map(path))
    index=jw_jobs_index(path)
    st_local("njobs",strtrim(strofreal(rows(index))))
    if (rows(index)==0) index=J(1,4,.)
    st_matrix(st_local("jobs"),index)
    st_matrix(st_local("settings"),jw_jobs_settings(path))
}

/* Joint testing uses the saved canonical contrasts without refitting,
   changing e(), or drawing random numbers. Calendar-date subsets map only
   to eligible optional contrasts; maintained baselines are never tested. */
void jw_pretrend_stata()
{
    real colvector z, ii
    real matrix V, hm, requested, selections
    real rowvector mask, result
    real scalar k, oi, i, h, indexblocks
    string scalar testset

    z = st_matrix("e(z)")'
    V = st_matrix("e(V_z)")
    k = st_numscalar("e(K)")
    indexblocks = st_global("e(restriction_scale)")=="conditional_index"
    mask = J(1,k,1)
    hm = J(0,2,.)
    if (k>0) hm = st_matrix("e(hmap)")
    if (rows(hm)!=k | cols(hm)!=2) {
        jw_core_error("saved optional-preperiod map is incompatible")
    }
    testset = st_local("testset")
    if (testset=="selected") {
        oi = strtoreal(st_local("output"))
        if (k>0) {
            selections = st_matrix("e(selection)")
            if (oi<1 | oi>rows(selections) | cols(selections)!=k) {
                jw_core_error("selected output has no compatible restriction mask")
            }
            mask = selections[oi,.]
        }
    }
    else if (testset=="custom") {
        requested = st_matrix(st_local("prematrix"))
        if (cols(requested)!=2 | hasmissing(requested)) {
            jw_core_error("prematrix() must contain finite cohort and calendar-time pairs")
        }
        mask = J(1,k,0)
        for (i=1; i<=rows(requested); i++) {
            if (indexblocks) {
                ii = vec(selectindex((hm[,1]:==requested[i,1]):&
                    (abs(hm[,2]:-requested[i,2]):<=1e-8*st_numscalar("e(gap)"))))
                if (rows(ii)==0)
                    jw_core_error("prematrix() contains an ineligible date or maintained baseline")
                if (any(mask[ii]:==1)) jw_core_error("prematrix() contains duplicate cohort-time pairs")
                mask[ii] = J(1,rows(ii),1)
            }
            else {
                h = jw_findcell(hm,requested[i,1],requested[i,2],st_numscalar("e(gap)"))
                if (h==0) {
                    jw_core_error("prematrix() contains an ineligible date or maintained baseline")
                }
                if (mask[h]==1) jw_core_error("prematrix() contains duplicate cohort-time pairs")
                mask[h] = 1
            }
        }
    }
    if (indexblocks) jw_block_masks_check(hm,mask)
    result = jw_pretrend_test(z,V,mask,k)
    st_matrix(st_local("testtable"),result)
    if (k>0) st_matrix(st_local("testmask"),mask)
    if (sum(mask)>0) {
        ii = vec(selectindex(mask'))
        if (indexblocks) st_matrix(st_local("testperiods"),uniqrows(sort(hm[ii,.],(1,2))))
        else st_matrix(st_local("testperiods"),hm[ii,.])
    }
}

/* Completion marker must be LAST: a partially compiled bridge is not ready. */
string scalar jw_version()
{
    return("0.8.8")
}

end
