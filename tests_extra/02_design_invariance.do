version 17.0
set more off

/* jwselect 0.5.0: design, units, normalization, and fixed-target checks.
   Run from the add-on test runner, or with jwselect already on adopath.
   Synthetic B and V describe 16 raw cohort/date estimates. V is their
   SAMPLING covariance. N_units=500 is the independent-unit count for BIC;
   do not multiply or divide V by that count.
   This file uses only public commands, not jwselect's Mata functions.
   No log is opened here. The runner collects output. */

capture mata: mata drop xdi_fixture()
capture mata: mata drop xdi_close()
capture mata: mata drop xdi_expected()
capture mata: mata drop xdi_transform()
capture mata: mata drop xdi_joint_scores()
capture mata: mata drop xdi_aggregate()
capture mata: mata drop xdi_custom_test()

mata:
void xdi_fixture()
{
    real scalar block, date, i, j, cohort
    real matrix cellmap, covmat
    real colvector coef, slope

    cellmap = J(16,2,.)
    coef = J(16,1,.)
    slope = 1 :+ (1::16):/10
    covmat = J(16,16,.)
    for (block=1; block<=2; block++)
    {
        cohort = 2+2*block
        for (date=1; date<=8; date++)
        {
            i = 8*(block-1)+date
            cellmap[i,.] = (cohort,date)
            coef[i] = .14*date + .08*sin((block+1)*date) + .3*block
            if (date>=cohort)
            {
                coef[i] = coef[i]+.8+.1*(date-cohort)
            }
        }
    }
    for (i=1; i<=16; i++)
    {
        for (j=1; j<=16; j++)
        {
            covmat[i,j] = .015*(.55^abs(i-j)) + .005*(i==j) +
                         .003*slope[i]*slope[j]
        }
    }
    /* One distinctive preperiod gives a mixed selected/not-selected mask. */
    coef[1] = coef[1]+.65
    st_matrix("XDI_B",coef')
    st_matrix("XDI_V",covmat)
    st_matrix("XDI_C",cellmap)
}

void xdi_close(string scalar first, string scalar second, real scalar tol)
{
    real matrix left, right
    real scalar scale
    left = st_matrix(first)
    right = st_matrix(second)
    assert(rows(left)==rows(right))
    assert(cols(left)==cols(right))
    assert(!hasmissing(left))
    assert(!hasmissing(right))
    scale = max((1,max(abs(left)),max(abs(right))))
    assert(max(abs(left-right))<=tol*scale)
}

/* Construct L directly from unit rows of the RAW 16-cell vector.
   Original calendar dates are 1,...,8; group dates are 4 and 6.
   Changing gap just expresses all dates in another unit. */
void xdi_expected(real scalar unitgap, real scalar delta, string scalar basetype)
{
    real scalar block, cohort, base, date, rawindex, baseindex
    real scalar pindex, hindex, npost, npre
    real matrix lpost, lpre, postmap, premap, basemap, linmap, cv
    real colvector zval

    npost = 8
    npre = 6-2*delta
    lpost = J(npost,16,0)
    lpre = J(npre,16,0)
    postmap = J(npost,2,.)
    premap = J(npre,2,.)
    basemap = J(2,2,.)
    pindex = 0
    hindex = 0
    for (block=1; block<=2; block++)
    {
        cohort = 2+2*block
        base = cohort-delta-1
        if (basetype=="common")
        {
            base = 4-delta-1
        }
        baseindex = 8*(block-1)+base
        basemap[block,.] = (cohort,base)*unitgap
        for (date=cohort; date<=8; date++)
        {
            pindex = pindex+1
            rawindex = 8*(block-1)+date
            lpost[pindex,rawindex] = 1
            lpost[pindex,baseindex] = -1
            postmap[pindex,.] = (cohort,date)*unitgap
        }
        for (date=1; date<cohort-delta; date++)
        {
            if (date!=base)
            {
                hindex = hindex+1
                rawindex = 8*(block-1)+date
                lpre[hindex,rawindex] = 1
                lpre[hindex,baseindex] = -1
                premap[hindex,.] = (cohort,date)*unitgap
            }
        }
    }
    assert(pindex==npost & hindex==npre)
    linmap = lpost \ lpre
    zval = linmap*st_matrix("XDI_B")'
    cv = linmap*st_matrix("XDI_V")*linmap'
    st_matrix("XDI_CG",st_matrix("XDI_C")*unitgap)
    st_matrix("XDI_EBASE",basemap)
    st_matrix("XDI_EP",postmap)
    st_matrix("XDI_EH",premap)
    st_matrix("XDI_EZ",zval')
    st_matrix("XDI_EV",cv)
    st_matrix("XDI_EBC",zval[1..npost]')
    st_matrix("XDI_EVC",cv[1..npost,1..npost])
}

/* Alternative representations of exactly the same canonical contrasts. */
void xdi_transform(string scalar method)
{
    real colvector coef, shift, perm
    real matrix covmat, cellmap, rebase
    real scalar block, i, refindex

    coef = st_matrix("XDI_B")'
    covmat = st_matrix("XDI_V")
    cellmap = st_matrix("XDI_C")
    if (method=="permutation")
    {
        perm = (16\1\9\4\12\6\14\3\11\8\2\15\5\13\7\10)
        coef = coef[perm]
        covmat = covmat[perm,perm]
        cellmap = cellmap[perm,.]
    }
    if (method=="shift")
    {
        shift = J(8,1,2.25) \ J(8,1,-1.75)
        coef = coef+shift
    }
    if (method=="rebase")
    {
        rebase = I(16)
        for (block=1; block<=2; block++)
        {
            refindex = 8*(block-1)+1
            for (i=refindex; i<=8*block; i++)
            {
                rebase[i,refindex] = rebase[i,refindex]-1
            }
        }
        coef = rebase*coef
        covmat = rebase*covmat*rebase'
    }
    if (method=="scale")
    {
        coef = 7*coef
        covmat = 49*covmat
    }
    st_matrix("XDI_BT",coef')
    st_matrix("XDI_VT",covmat)
    st_matrix("XDI_CT",cellmap)
}

/* Only restriction count, Wald, and WBIC are target-independent. Other
   score columns contain estimates/variances for the aggregation target. */
void xdi_joint_scores(string scalar destination)
{
    real matrix scores
    scores = st_matrix("e(scores)")
    scores = select(scores,scores[,1]:==2)
    st_matrix(destination,scores[,4..6])
}

/* Compute all reported cohort/aggregate targets independently from cells.
   This checks off-diagonal covariance, not only the aggregate's SE. */
void xdi_aggregate(real rowvector weights)
{
    real matrix postmap, targets, cellcov, expectedcov
    real rowvector cellb, expectedb
    real scalar i, g, countg

    postmap = st_matrix("e(postmap)")
    cellb = st_matrix("e(b_cells)")
    cellcov = st_matrix("e(V_cells)")
    assert(cols(weights)==rows(postmap))
    targets = J(3,rows(postmap),0)
    for (g=1; g<=2; g++)
    {
        countg = sum(postmap[,1]:==(2+2*g))
        for (i=1; i<=rows(postmap); i++)
        {
            if (postmap[i,1]==2+2*g)
            {
                targets[g,i] = 1/countg
            }
        }
    }
    targets[3,.] = weights
    expectedb = cellb*targets'
    expectedcov = targets*cellcov*targets'
    st_matrix("XDI_EOUTB",expectedb)
    st_matrix("XDI_EOUTV",expectedcov)
    st_matrix("XDI_EW",weights)
    xdi_close("e(b)","XDI_EOUTB",1e-9)
    xdi_close("e(V)","XDI_EOUTV",1e-9)
    xdi_close("e(weights)","XDI_EW",1e-12)
}

/* For common baseline 3, (6,5) is an eligible optional date even though
   it is later than the baseline. Its position in hmap is six. */
void xdi_custom_test()
{
    real matrix cv
    real rowvector zval
    real scalar waldval
    zval = st_matrix("XDI_EZ")
    cv = st_matrix("XDI_EV")
    waldval = zval[14]^2/cv[14,14]
    st_numscalar("XDI_WTEST",waldval)
    st_matrix("XDI_EMASK",(0,0,0,0,0,1))
}
end

mata: xdi_fixture()
clear
quietly set obs 3
generate byte xdi_marker = _n

/* Cases 1-18: independent canonical maps and covariance in every design.
   No optional restrictions are imposed here, so cell estimates must be
   exactly the independently constructed post-minus-baseline contrasts. */
local designcases 0
foreach gap in 1 2 .25 {
    foreach base in cohort common {
        foreach delta in 0 1 2 {
            mata: xdi_expected(strtoreal(st_local("gap")),strtoreal(st_local("delta")),st_local("base"))
            quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_CG) ///
                nunits(500) gap(`gap') baseline(`base') anticipation(`delta') ///
                select(manual)
            assert e(P)==8
            assert e(K)==6-2*`delta'
            assert e(G)==2
            assert e(gap)==`gap'
            assert e(N_units)==500
            assert e(anticipation)==`delta'
            mata: xdi_close("e(baselines)","XDI_EBASE",1e-12)
            mata: xdi_close("e(postmap)","XDI_EP",1e-12)
            mata: xdi_close("e(hmap)","XDI_EH",1e-12)
            mata: xdi_close("e(z)","XDI_EZ",1e-10)
            mata: xdi_close("e(V_z)","XDI_EV",1e-10)
            mata: xdi_close("e(b_cells)","XDI_EBC",1e-10)
            mata: xdi_close("e(V_cells)","XDI_EVC",1e-10)
            local designcases = `designcases'+1
        }
    }
}
assert `designcases'==18
display as text "  18 baseline/anticipation/gap designs checked"

/* Cases 19-20: same FULL preperiod set, with a different maintained base.
   The full set is {1,2,3} for cohort 4 and {1,2,3,4,5} for cohort 6.
   prematrix() lists that set excluding the relevant maintained baseline.
   These are fixed manual restrictions. This test makes NO claim that
   nested Wald-BIC selection must be invariant to a baseline change. */
mata: xdi_expected(1,0,"cohort")
matrix XDI_PRE = XDI_EH
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(cohort) select(manual) prematrix(XDI_PRE)
matrix XDI_FULLB = e(b)
matrix XDI_FULLV = e(V)
matrix XDI_FULLBC = e(b_cells)
matrix XDI_FULLVC = e(V_cells)
mata: xdi_expected(1,0,"common")
matrix XDI_PRE = XDI_EH
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(common) select(manual) prematrix(XDI_PRE)
mata: xdi_close("e(b)","XDI_FULLB",1e-9)
mata: xdi_close("e(V)","XDI_FULLV",1e-9)
mata: xdi_close("e(b_cells)","XDI_FULLBC",1e-9)
mata: xdi_close("e(V_cells)","XDI_FULLVC",1e-9)

/* Case 21: explicit calendar dates after a common baseline. */
matrix XDI_CUSTOMP = (6,5)
matrix XDI_SAVEB = e(b)
matrix XDI_SAVEV = e(V)
mata: xdi_custom_test()
quietly estat pretrend, prematrix(XDI_CUSTOMP)
assert r(df)==1
assert abs(r(chi2)-scalar(XDI_WTEST))<1e-9
assert abs(r(J)-scalar(XDI_WTEST))<1e-9
assert abs(r(p)-chi2tail(1,scalar(XDI_WTEST)))<1e-10
assert r(postselection)==0
assert "`r(test_set)'"=="custom"
local returned_matrices : r(matrices)
assert strpos(" `returned_matrices' "," mask ")>0
assert strpos(" `returned_matrices' "," preperiods ")>0
matrix XDI_GOTMASK = r(mask)
matrix XDI_GOTPRE = r(preperiods)
mata: xdi_close("XDI_GOTMASK","XDI_EMASK",1e-12)
mata: xdi_close("XDI_GOTPRE","XDI_CUSTOMP",1e-12)
mata: xdi_close("e(b)","XDI_SAVEB",1e-12)
mata: xdi_close("e(V)","XDI_SAVEV",1e-12)

/* Case 22: joint WBIC reference. Check the fixture is not at a score tie,
   because invariance assertions must not hinge on rounding at a boundary.
   The documented tie rule is the first row in the candidate family. */
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(cohort) select(wbic) wbicscope(joint) candidates(all)
matrix XDI_REFB = e(b)
matrix XDI_REFV = e(V)
matrix XDI_REFBC = e(b_cells)
matrix XDI_REFVC = e(V_cells)
matrix XDI_REFSEL = e(selection)
matrix XDI_REFZ = e(z)
matrix XDI_REFZV = e(V_z)
matrix XDI_REFFAM = e(joint_candidates)
mata: xdi_joint_scores("XDI_REFSCORE")
mata: st_matrix("XDI_SORTSCORE",sort(st_matrix("XDI_REFSCORE"),(3)))
assert XDI_SORTSCORE[2,3]-XDI_SORTSCORE[1,3]>1e-5
assert rowsof(XDI_REFSCORE)==64
mata: assert(sum(st_matrix("XDI_REFSEL")[1,.])==5)
mata: xdi_aggregate((J(1,5,1/10),J(1,3,1/6)))
quietly estat pretrend
scalar XDI_REFWALD = r(chi2)
scalar XDI_REFJ = r(J)

/* Cases 23-25: permutation, cohort shifts, and a covariance-aware raw
   normalization change. Rebased raw V is positive semidefinite; the
   identifiable contrast covariance stays positive definite. */
foreach transform in permutation shift rebase {
    mata: xdi_transform(st_local("transform"))
    quietly jwselect, bmat(XDI_BT) vmat(XDI_VT) cells(XDI_CT) nunits(500) ///
        baseline(cohort) select(wbic) wbicscope(joint) candidates(all)
    mata: xdi_close("e(b)","XDI_REFB",1e-9)
    mata: xdi_close("e(V)","XDI_REFV",1e-9)
    mata: xdi_close("e(b_cells)","XDI_REFBC",1e-9)
    mata: xdi_close("e(V_cells)","XDI_REFVC",1e-9)
    mata: xdi_close("e(z)","XDI_REFZ",1e-9)
    mata: xdi_close("e(V_z)","XDI_REFZV",1e-9)
    mata: xdi_close("e(selection)","XDI_REFSEL",1e-12)
    mata: xdi_close("e(joint_candidates)","XDI_REFFAM",1e-12)
    mata: xdi_joint_scores("XDI_NEWSCORE")
    mata: xdi_close("XDI_NEWSCORE","XDI_REFSCORE",1e-9)
    quietly estat pretrend
    assert abs(r(chi2)-scalar(XDI_REFWALD))<1e-8
    assert abs(r(J)-scalar(XDI_REFJ))<1e-8
}

/* Case 26: change the outcome's unit by a factor of seven. */
mata: xdi_transform("scale")
quietly jwselect, bmat(XDI_BT) vmat(XDI_VT) cells(XDI_CT) nunits(500) ///
    baseline(cohort) select(wbic) wbicscope(joint) candidates(all)
matrix XDI_SCALEB = 7*XDI_REFB
matrix XDI_SCALEV = 49*XDI_REFV
matrix XDI_SCALEBC = 7*XDI_REFBC
matrix XDI_SCALEVC = 49*XDI_REFVC
mata: xdi_close("e(b)","XDI_SCALEB",1e-9)
mata: xdi_close("e(V)","XDI_SCALEV",1e-9)
mata: xdi_close("e(b_cells)","XDI_SCALEBC",1e-9)
mata: xdi_close("e(V_cells)","XDI_SCALEVC",1e-9)
mata: xdi_close("e(selection)","XDI_REFSEL",1e-12)
mata: xdi_joint_scores("XDI_NEWSCORE")
mata: xdi_close("XDI_NEWSCORE","XDI_REFSCORE",1e-9)
quietly estat pretrend
assert abs(r(chi2)-scalar(XDI_REFWALD))<1e-8
assert abs(r(J)-scalar(XDI_REFJ))<1e-8

/* Cases 27-29: target weights and post horizon do not change h, so joint
   Wald/BIC scores and selection must stay fixed. The reporting targets do
   change. Check aggregate covariance, including every cross-cohort term. */
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(cohort) select(wbic) wbicscope(joint) candidates(all) ///
    aggregate(equalcell)
mata: xdi_aggregate(J(1,8,1/8))
mata: xdi_close("e(selection)","XDI_REFSEL",1e-12)
mata: xdi_joint_scores("XDI_NEWSCORE")
mata: xdi_close("XDI_NEWSCORE","XDI_REFSCORE",1e-9)

matrix XDI_W = (4,4,.05 \ 4,5,.10 \ 4,6,.15 \ 4,7,.10 \ ///
               4,8,.10 \ 6,6,.20 \ 6,7,.25 \ 6,8,.05)
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(cohort) select(wbic) wbicscope(joint) candidates(all) ///
    aggweights(XDI_W)
mata: xdi_aggregate((.05,.10,.15,.10,.10,.20,.25,.05))
mata: xdi_close("e(b_cells)","XDI_REFBC",1e-9)
mata: xdi_close("e(V_cells)","XDI_REFVC",1e-9)
mata: xdi_close("e(selection)","XDI_REFSEL",1e-12)
mata: xdi_joint_scores("XDI_NEWSCORE")
mata: xdi_close("XDI_NEWSCORE","XDI_REFSCORE",1e-9)

matrix XDI_SHORTMAP = (4,4 \ 4,5 \ 6,6 \ 6,7)
matrix XDI_SHORTW = (4,4,.10 \ 4,5,.20 \ 6,6,.30 \ 6,7,.40)
quietly jwselect, bmat(XDI_B) vmat(XDI_V) cells(XDI_C) nunits(500) ///
    baseline(cohort) select(wbic) wbicscope(joint) candidates(all) ///
    postwindow(0 1) aggweights(XDI_SHORTW)
assert e(P)==4
assert e(K)==6
mata: xdi_close("e(postmap)","XDI_SHORTMAP",1e-12)
mata: st_matrix("XDI_SHORTBC",st_matrix("XDI_REFBC")[.,(1,2,6,7)])
mata: st_matrix("XDI_SHORTVC",st_matrix("XDI_REFVC")[(1,2,6,7),(1,2,6,7)])
mata: xdi_close("e(b_cells)","XDI_SHORTBC",1e-9)
mata: xdi_close("e(V_cells)","XDI_SHORTVC",1e-9)
mata: xdi_aggregate((.10,.20,.30,.40))
mata: xdi_close("e(selection)","XDI_REFSEL",1e-12)
mata: xdi_joint_scores("XDI_NEWSCORE")
mata: xdi_close("XDI_NEWSCORE","XDI_REFSCORE",1e-9)

/* Cases 30-35: invalid designs fail and preserve active e(), RNG, data.
   Test the known validation error (3498), not merely any error code. */
matrix XDI_HELDB = e(b)
matrix XDI_HELDV = e(V)
matrix XDI_HELDSEL = e(selection)
local rng_before `c(rngstate)'
local failures 0
foreach invalid in grid missingbase postwindow nonsymmetric indefinite duplicate {
    matrix XDI_BADV = XDI_V
    matrix XDI_BADC = XDI_C
    local badoptions
    if "`invalid'"=="grid" matrix XDI_BADC[1,2] = 1.5
    if "`invalid'"=="missingbase" matrix XDI_BADC[3,2] = 0
    if "`invalid'"=="postwindow" local badoptions postwindow(0 3)
    if "`invalid'"=="nonsymmetric" matrix XDI_BADV[1,2] = XDI_BADV[1,2]+.05
    if "`invalid'"=="indefinite" matrix XDI_BADV[1,1] = -1
    if "`invalid'"=="duplicate" matrix XDI_BADC[2,2] = 1
    capture noisily jwselect, bmat(XDI_B) vmat(XDI_BADV) cells(XDI_BADC) ///
        nunits(500) baseline(cohort) select(wbic) wbicscope(joint) ///
        candidates(all) `badoptions'
    local failed_rc = _rc
    assert `failed_rc'==3498
    assert "`c(rngstate)'"=="`rng_before'"
    assert "`e(cmd)'"=="jwselect"
    assert "`e(wbic_scope)'"=="joint"
    assert e(P)==4
    assert e(K)==6
    assert e(N_units)==500
    mata: xdi_close("e(b)","XDI_HELDB",1e-12)
    mata: xdi_close("e(V)","XDI_HELDV",1e-12)
    mata: xdi_close("e(selection)","XDI_HELDSEL",1e-12)
    assert _N==3
    assert xdi_marker==_n
    local failures = `failures'+1
}
assert `failures'==6
display as text "  35 design and invariance cases completed"
display as result "PASS 02_design_invariance"
