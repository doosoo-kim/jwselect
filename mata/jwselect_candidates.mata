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
