*! Actual Mata check: partial projection-cache prefixes and mixed components.
*! Requires loaded jwselect 0.3.0; no RNG draws or data changes.
version 17.0
capture mata: mata drop jw030_partial_cache()
mata:
void jw030_partial_cache()
{
    struct jw_rs_context scalar C
    real colvector mu, reference, actual
    real rowvector prefixes
    real scalar j, prefix, mb
    // colshape's second argument sets columns; these loadings are 340 x 3.
    C.selectB = colshape(sin((1::1020):/7),3)
    C.reportB = colshape(cos((1::1020):/11),3)
    C.constant = sin((1::340):/3)
    C.component = mod((0::339),2):+1
    C.groups = (1 \ 2)
    C.EH = colshape(cos((1::123):/13),3)
    C.baseline = sin((1::41):/5)
    C.projection = J(0,0,.)
    mu = (.12 \ -.43 \ .71)
    assert(rows(C.selectB)==340 & cols(C.selectB)==3)
    assert(rows(C.reportB)==340 & cols(C.reportB)==3)
    assert(rows(C.EH)==41 & cols(C.EH)==3)
    assert(rows(C.constant)==340 & rows(C.component)==340)
    assert(rows(C.baseline)==41 & rows(mu)==3)
    reference = jw_rs_draw_errors(C,mu)
    prefixes = (0,1,63,64,65,170,339,340)
    for (j=1; j<=cols(prefixes); j++) {
        prefix = prefixes[j]
        // Fractional MB, with a quarter-column margin against rounding down.
        mb = (8*rows(C.EH)*(prefix+.25))/1000000
        if (prefix==0) mb = 0
        jw_rs_cache(C,mb)
        assert(cols(C.projection)==prefix)
        actual = jw_rs_draw_errors(C,mu)
        assert(max(abs(actual-reference))<1e-10)
    }
}
jw030_partial_cache()
mata drop jw030_partial_cache()
end
display as result "Partial cache prefixes 0/1/63/64/65/170/339/340, interleaved components: PASS"
