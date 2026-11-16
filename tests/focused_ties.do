*! Focused tie rule: fewer restrictions, then original family order.
*! Run from the package root; leaves data, estimates, and RNG unchanged.
version 16.1
capture mata: assert(findexternal("jw_focused_order()") != NULL)
if _rc quietly do _jwselect_mata.do
capture mata: mata drop jw_test_focused_ties()
mata:
void jw_test_focused_ties()
{
    struct jw_model scalar M
    struct jw_family scalar C
    struct jw_fit scalar F
    struct jw_rs_context scalar RS
    real matrix supplied, E
    real rowvector expected
    real colvector mu, err, ii
    real scalar s, o, cache, equalrank
    string rowvector focus

    M.z = (4 \ 1 \ 0 \ -1)
    M.V = (4,1,0,1 \ 1,1,0,0 \ 0,0,1,0 \ 1,0,0,1)
    M.p = 1
    M.k = 3
    M.postmap = (5,5)
    M.hmap = (5,1 \ 5,2 \ 5,3)
    M.groups = 5
    M.bases = 4
    M.T = (1,0,0,0)
    M.w = 1
    M.nunit = 100
    M.gap = 1
    C.owners = (5 \ 5 \ 5)
    C.ids = (1 \ 2 \ 3)
    focus = ("cohort","aggregate","both")
    E = J(3,4,0)
    E[2,1] = 2
    E[3,1] = -3
    mu = (1 \ 0 \ -1)

    for (equalrank=0; equalrank<=1; equalrank++) {
        /* Baseline score=4. Both alternatives score=3 exactly, but their
           estimates are 3 and 5, so a wrong tie rule changes the answer. */
        C.masks = (0,0,0 \ 1,1,0 \ 0,equalrank,1)
        supplied = C.masks
        expected = (0,0,1)
        if (equalrank) expected = (1,1,0)
        for (s=1; s<=3; s++) {
            F = jw_run(M,C,supplied,"focused",focus[s],100)
            assert(F.scores[2,11]==F.scores[3,11])
            for (o=1; o<=rows(F.A); o++) {
                assert(sum(abs(F.selected[o,.]-expected))==0)
                assert(F.A[o,.]*M.z==5-2*equalrank)
                /* Simulated first-minimum must use the same rule. Test
                   streaming, cached, and the public deterministic helper. */
                err = E[,1]:+(1-2*equalrank)
                assert(max(abs(jw_rs_errors(M,F,o,E,mu)-err))==0)
                for (cache=0; cache<=1; cache++) {
                    RS = jw_rs_prepare_opts(M,F,o,E,(2,3,4),cache)
                    assert(max(abs(jw_rs_draw_errors(RS,mu)-err))==0)
                }
            }
        }
    }

    /* Near ties remain strict comparisons; there is no numerical tolerance
       that could replace a slightly smaller score by a lower-rank model. */
    C.masks = (0,0,0 \ 1,1,0 \ 0,0,1)
    M.z[4] = -1-2^(-25)
    F = jw_run(M,C,C.masks,"focused","both",100)
    for (o=1; o<=rows(F.A); o++) {
        assert(sum(abs(F.selected[o,.]-(1,1,0)))==0)
    }

    /* WBIC did not acquire the focused tie rule. Its preexisting
       first-in-family convention is unchanged at an exact negative tie. */
    M.z = (4 \ 1 \ 0 \ 0)
    M.nunit = exp(1)
    F = jw_run_scope(M,C,C.masks,"wbic","both",100,"cohort")
    assert(F.scores[2,6]==F.scores[3,6] & F.scores[2,6]<0)
    assert(sum(abs(F.selected[1,.]-(1,1,0)))==0)
    F = jw_run_scope(M,C,C.masks,"wbic","both",100,"joint")
    assert(sum(abs(F.selected[1,.]-(1,1,0)))==0)

    /* Baseline first among zero-restriction ties; empty models supported. */
    ii = jw_focused_order(J(0,0,.))
    assert(rows(ii)==0 & cols(ii)==1)
    ii = jw_focused_order((1,1,0 \ 0,0,1 \ 1,0,0 \ 0,0,0))
    assert(sum(abs(ii-(4 \ 2 \ 3 \ 1)))==0)
}
jw_test_focused_ties()
end
display as result "Focused exact ties, all scopes and cached/streamed reselection: PASS"
