*! Shared missing-aware matrix comparisons for release tests.
*! Mata missing(X) counts missing entries; use colon comparisons for masks.
version 17.0
capture mata: mata drop jw030_match()
capture mata: mata drop jw030_equal()
capture mata: mata drop jw030_assertion_checks()
mata:
real scalar jw030_match(real matrix A, real matrix B, real scalar tol)
{
    real colvector good, difference
    if (rows(A)!=rows(B) | cols(A)!=cols(B)) return(0)
    if (missing(tol) | tol<=0) return(0)
    if (rows(A)*cols(A)==0) return(1)
    if (!all(vec(A:>=.):==vec(B:>=.))) return(0)
    good = vec(selectindex(vec(A:<.)))
    if (rows(good)*cols(good)==0) return(1)
    difference = vec(A)-vec(B)
    difference = difference[good]
    if (hasmissing(difference)) return(0)
    return(max(abs(difference))<tol)
}

void jw030_equal(string scalar left, string scalar right, real scalar tol)
{
    real matrix A, B
    real scalar matches
    A = st_matrix(left)
    B = st_matrix(right)
    matches = jw030_match(A,B,tol)
    if (!matches) {
        printf("Matrix comparison failed: %s versus %s, tolerance %g\n",
               left,right,tol)
        printf("Left matrix (%g x %g):\n",rows(A),cols(A))
        A
        printf("Right matrix (%g x %g):\n",rows(B),cols(B))
        B
    }
    assert(matches)
}

void jw030_assertion_checks()
{
    real matrix A, B
    A = (1,2 \ 3,4)
    assert(jw030_match(A,A,1e-10))
    B = A
    B[2,2] = 5
    assert(!jw030_match(A,B,1e-10))
    A = (1,. \ 3,4)
    assert(jw030_match(A,A,1e-10))
    B = A
    B[2,2] = 5
    assert(!jw030_match(A,B,1e-10))
    B = (.,2 \ 3,4)
    assert(!jw030_match(A,B,1e-10))
    assert(jw030_match(J(2,2,.),J(2,2,.),1e-10))
    assert(!jw030_match(J(2,2,.),J(2,2,0),1e-10))
    assert(jw030_match(3,3,1e-10))
    assert(!jw030_match(3,4,1e-10))
    assert(jw030_match(.,.,1e-10))
    assert(!jw030_match(.,3,1e-10))
    assert(jw030_match(J(0,2,.),J(0,2,.),1e-10))
    assert(jw030_match(J(2,0,.),J(2,0,.),1e-10))
    assert(!jw030_match(J(0,2,.),J(0,3,.),1e-10))
    assert(!jw030_match((1,2),(1 \ 2),1e-10))
    assert(jw030_match(1,1+1e-6,1e-5))
    assert(!jw030_match(1,1+1e-6,1e-8))
    assert(!jw030_match(1,1,0))
    assert(!jw030_match(1,1,.))
}
jw030_assertion_checks()
end
display as result "Matrix comparison helper: PASS"
