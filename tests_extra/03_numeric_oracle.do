*! Extra numerical tests for jwselect 0.5.0, Stata 17 or newer.
*! Run from the folder containing tests_extra/: do tests_extra/03_numeric_oracle.do
*! Uses only public jwselect/estat commands and precomputed independent fixtures.
*! No Python is required. Does not clear data or open/close a log.
version 17.0
set more off

capture mata: mata drop xo_near()
capture mata: mata drop xo_equal()
capture mata: mata drop xo_checkfit()
capture mata: mata drop xo_scores()
capture mata: mata drop xo_selfcheck()

mata:
real scalar xo_near(real matrix A, real matrix B,
                   real scalar atol, real scalar rtol)
{
    real colvector idx, av, bv, delta, bound
    if (rows(A)!=rows(B) | cols(A)!=cols(B)) return(0)
    if (missing(atol) | missing(rtol) | atol<0 | rtol<0) return(0)
    if (rows(A)*cols(A)==0) return(1)
    if (!all(vec(A:>=.):==vec(B:>=.))) return(0)
    idx = vec(selectindex(vec(A:<.)))
    if (rows(idx)*cols(idx)==0) return(1)
    av = vec(A)
    bv = vec(B)
    av = av[idx]
    bv = bv[idx]
    delta = abs(av-bv)
    bound = atol :+ rtol:*(abs(av)+abs(bv))
    if (hasmissing(delta) | hasmissing(bound)) return(0)
    return(all(delta:<=bound))
}

void xo_equal(string scalar actual, string scalar expected)
{
    real matrix A, B
    real scalar okay
    A = st_matrix(actual)
    B = st_matrix(expected)
    okay = xo_near(A,B,2e-9,2e-9)
    if (!okay) {
        printf("Oracle mismatch: %s versus %s\n",actual,expected)
        printf("Actual (%g x %g):\n",rows(A),cols(A))
        A
        printf("Expected (%g x %g):\n",rows(B),cols(B))
        B
    }
    assert(okay)
}

void xo_checkfit(string scalar prefix)
{
    real matrix info, chosen, got
    real scalar g, p, k
    info = st_matrix("XO_INFO")
    g = info[1,2]
    p = info[1,3]
    k = info[1,4]
    chosen = st_matrix(prefix+"MASK")
    got = st_matrix("e(selection)")
    assert(rows(got)==g+1 & cols(got)==k)
    assert(xo_near(got,J(g+1,1,1)*chosen,0,0))
    assert(st_numscalar("e(P)")==p)
    assert(st_numscalar("e(K)")==k)
    assert(st_numscalar("e(G)")==g)
    assert(st_numscalar("e(N_units)")==info[1,5])
    xo_equal("e(b)",prefix+"B")
    xo_equal("e(V)",prefix+"V")
    xo_equal("e(b_cells)",prefix+"BC")
    xo_equal("e(V_cells)",prefix+"VC")
    xo_equal("e(loading)",prefix+"LOAD")
    xo_equal("e(z)","XO_Z")
    xo_equal("e(V_z)","XO_VZ")
    xo_equal("e(hmap)","XO_HMAP")
    xo_equal("e(postmap)","XO_POSTMAP")
    xo_equal("e(targets)","XO_TARGETS")
}

void xo_scores()
{
    real matrix E, F, S, expected, actual
    real colvector score_rows, found, used
    real scalar k, r, candidate, oracle_row
    E = st_matrix("XO_SCORES")
    F = st_matrix("e(joint_candidates)")
    S = st_matrix("e(scores)")
    k = cols(F)
    assert(cols(E)==k+4)
    assert(rows(F)==rows(E))
    assert(all(vec((F:==0):|(F:==1))))
    score_rows = vec(selectindex(S[,1]:==2))
    assert(rows(score_rows)==rows(E))
    used = J(rows(E),1,0)
    for (r=1; r<=rows(score_rows); r++) {
        candidate = S[score_rows[r],3]
        assert(candidate>=1 & candidate<=rows(F))
        assert(candidate==floor(candidate))
        assert(S[score_rows[r],2]==0)
        found = vec(selectindex(rowsum(E[,1..k]:!=
                    J(rows(E),1,1)*F[candidate,.]):==0))
        assert(rows(found)==1)
        oracle_row = found[1]
        used[oracle_row] = used[oracle_row]+1
        actual = S[score_rows[r],4..6]
        expected = (E[oracle_row,k+1],E[oracle_row,k+2],
                    E[oracle_row,k+4])
        if (!xo_near(actual,expected,2e-8,2e-9)) {
            printf("Joint score mismatch at oracle row %g\n",oracle_row)
            actual
            expected
        }
        assert(xo_near(actual,expected,2e-8,2e-9))
    }
    // A duplicate row must not substitute for an omitted candidate.
    assert(all(used:==1))
}

void xo_selfcheck()
{
    assert(xo_near((0,.,2),(0,.,2),2e-9,2e-9))
    assert(!xo_near((0,.,2),(0,3,2),2e-9,2e-9))
    assert(!xo_near((0,3,2),(0,.,2),2e-9,2e-9))
    assert(!xo_near((0,.,2),(0,.,3),2e-9,2e-9))
    assert(!xo_near((1,2),(1\2),2e-9,2e-9))
    assert(xo_near(0,1e-10,2e-9,2e-9))
    assert(!xo_near(0,1e-5,2e-9,2e-9))
    assert(xo_near(1e6,1e6+.001,2e-9,2e-9))
    assert(!xo_near(1e6,1e6+1,2e-9,2e-9))
}
xo_selfcheck()
end

forvalues xo_case = 1/24 {
    display as text "Numerical oracle case `xo_case'/24: loading fixture"
    quietly do tests_extra/numeric_oracle_fixture.do `xo_case'
    local xo_n = XO_INFO[1,5]
    display as text "  Fixture loaded (N = `xo_n'); fitting joint WBIC"

    quietly jwselect, bmat(XO_B) vmat(XO_V) cells(XO_C) ///
        nunits(`xo_n') select(wbic) wbicscope(joint) candidates(all)
    assert "`e(version)'" == "0.5.0"
    assert "`e(wbic_scope)'" == "joint"
    display as text "  Checking joint estimates and candidate scores"
    mata: xo_checkfit("XO_J")
    mata: xo_scores()

    display as text "  Checking joint fit: all-period pretrend test"
    quietly estat pretrend
    matrix XO_GOTTEST = (r(df),r(chi2),r(J),r(p))
    matrix XO_GOTTABLE = r(table)
    matrix XO_GOTMASK = r(mask)
    assert "`r(test_set)'" == "all"
    assert r(postselection) == 0
    mata: xo_equal("XO_GOTTEST","XO_ALLTEST")
    mata: xo_equal("XO_GOTTABLE","XO_ALLTEST")
    mata: assert(xo_near(st_matrix("XO_GOTMASK"),J(1,cols(st_matrix("XO_JMASK")),1),0,0))

    display as text "  Checking joint fit: selected-period pretrend test"
    quietly estat pretrend, selected
    matrix XO_GOTTEST = (r(df),r(chi2),r(J),r(p))
    matrix XO_GOTTABLE = r(table)
    matrix XO_GOTMASK = r(mask)
    assert "`r(test_set)'" == "selected"
    assert r(postselection) == 1
    mata: xo_equal("XO_GOTTEST","XO_JTEST")
    mata: xo_equal("XO_GOTTABLE","XO_JTEST")
    mata: xo_equal("XO_GOTMASK","XO_JMASK")
    if XO_JTEST[1,1] == 0 {
        assert r(df) == 0
        assert r(chi2) == 0
        assert r(J) == 0
        assert missing(r(p))
    }

    // Separately minimize each cohort's four local candidates, then verify
    // that the public estimates jointly impose their selected union.
    display as text "  Fitting cohort WBIC"
    quietly jwselect, bmat(XO_B) vmat(XO_V) cells(XO_C) ///
        nunits(`xo_n') select(wbic) wbicscope(cohort) candidates(all)
    assert "`e(wbic_scope)'" == "cohort"
    display as text "  Checking cohort estimates"
    mata: xo_checkfit("XO_C")

    display as text "  Checking cohort fit: all-period pretrend test"
    quietly estat pretrend
    matrix XO_GOTTEST = (r(df),r(chi2),r(J),r(p))
    matrix XO_GOTTABLE = r(table)
    matrix XO_GOTMASK = r(mask)
    assert "`r(test_set)'" == "all"
    assert r(postselection) == 0
    mata: xo_equal("XO_GOTTEST","XO_ALLTEST")
    mata: xo_equal("XO_GOTTABLE","XO_ALLTEST")
    mata: assert(xo_near(st_matrix("XO_GOTMASK"),J(1,cols(st_matrix("XO_CMASK")),1),0,0))

    display as text "  Checking cohort fit: selected-period pretrend test"
    quietly estat pretrend, selected
    matrix XO_GOTTEST = (r(df),r(chi2),r(J),r(p))
    matrix XO_GOTTABLE = r(table)
    matrix XO_GOTMASK = r(mask)
    assert "`r(test_set)'" == "selected"
    assert r(postselection) == 1
    mata: xo_equal("XO_GOTTEST","XO_CTEST")
    mata: xo_equal("XO_GOTTABLE","XO_CTEST")
    mata: xo_equal("XO_GOTMASK","XO_CMASK")
    if XO_CTEST[1,1] == 0 {
        assert r(df) == 0
        assert r(chi2) == 0
        assert r(J) == 0
        assert missing(r(p))
    }
    display as result "  PASS numerical oracle case `xo_case'/24"
}

capture mata: mata drop xo_near()
capture mata: mata drop xo_equal()
capture mata: mata drop xo_checkfit()
capture mata: mata drop xo_scores()
capture mata: mata drop xo_selfcheck()
display as result "PASS 03_numeric_oracle"
