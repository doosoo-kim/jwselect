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
