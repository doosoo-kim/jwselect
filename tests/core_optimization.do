*! Compare optimized Mata routines directly to frozen0.2 source paths.
*! Run from package root. Does not clear the dataset or existing estimates.
version 16.1
set more off
capture mata: assert(findexternal("jw_run()") != NULL)
if _rc quietly do _jwselect_mata.do
do tests/core_optimization_reference.do
capture mata: mata drop jw_test_core_equal()
capture mata: mata drop jw_test_core_fit()
capture mata: mata drop jw_test_core_modes()
capture mata: mata drop jw_test_core_zero_groups()
mata:
void jw_test_core_equal(real matrix a, real matrix b)
{
    assert(rows(a)==rows(b) & cols(a)==cols(b))
    if (rows(a)*cols(a)>0) assert(max(abs(a-b))<1e-10)
}

void jw_test_core_fit(struct jw_fit scalar a, struct jw_fit scalar b)
{
    jw_test_core_equal(a.A,b.A)
    jw_test_core_equal(a.targets,b.targets)
    jw_test_core_equal(a.selected,b.selected)
    jw_test_core_equal(a.scores,b.scores)
    jw_test_core_equal(a.cellA,b.cellA)
    jw_test_core_equal(a.cell_selected,b.cell_selected)
    jw_test_core_equal(a.allA,b.allA)
    jw_test_core_equal(a.allT,b.allT)
    jw_test_core_equal(a.allmeta,b.allmeta)
    jw_test_core_equal(a.selA,b.selA)
    jw_test_core_equal(a.selT,b.selT)
    jw_test_core_equal(a.outmap,b.outmap)
}

void jw_test_core_modes()
{
    struct jw_model scalar M
    struct jw_family scalar f, oldf
    struct jw_fit scalar a, b
    real matrix q, prem, supplied, empty
    real colvector rev
    real scalar G, K, gi, ki, p, m, i, h, s, j, c, trial, fam, method, scope, calls
    string rowvector families, methods, scopes
    real rowvector kvals, A, stat, target, mask

    rseed(63004)
    kvals = (0,1,4)
    families = ("nested","all")
    methods = ("manual","wbic","focused")
    scopes = ("cohort","aggregate","both")
    calls = 0
    for (gi=1; gi<=3; gi++) {
        G = gi
        for (ki=1; ki<=cols(kvals); ki++) {
            K = kvals[ki]
            p = 2*G
            m = p+K
            M.p = p
            M.k = K
            M.nunit = 500
            M.gap = 1
            M.groups = (10::(9+G))
            M.bases = M.groups:-1
            M.postmap = J(p,2,.)
            M.hmap = J(K,2,.)
            M.T = J(G,m,0)
            M.w = (1..p)/sum(1..p)
            for (i=1; i<=G; i++) {
                M.postmap[2*i-1,.] = (M.groups[i],M.groups[i])
                M.postmap[2*i,.] = (M.groups[i],M.groups[i]+1)
                M.T[i,(2*i-1)..(2*i)] = (.5,.5)
            }
            for (h=1; h<=K; h++) M.hmap[h,.] = (M.groups[mod(h-1,G)+1],h)
            q = rnormal(m,m,0,1)
            M.V = q*q'+I(m)
            M.z = rnormal(m,1,0,1)
            empty = J(0,0,.)
            for (fam=1; fam<=2; fam++) {
                for (method=1; method<=3; method++) {
                    prem = J(0,0,.)
                    if (method==1 & K>0) {
                        prem = M.hmap[selectindex(mod((1::K),2):==1),.]
                    }
                    f = jw_candidates(M,methods[method],families[fam],prem,empty,256)
                    oldf = jw020_candidates(M,methods[method],families[fam],prem,empty,256)
                    jw_test_core_equal(f.masks,oldf.masks)
                    jw_test_core_equal(f.owners,oldf.owners)
                    jw_test_core_equal(f.ids,oldf.ids)
                    /* Reverse supplied row order to test stable tie/order behavior.
                       This changes IDs/order deliberately but identically in both paths. */
                    for (trial=1; trial<=2; trial++) {
                        if (trial==2) {
                            rev = rows(f.masks)::1
                            f.masks = f.masks[rev,.]
                            f.owners = f.owners[rev,.]
                            f.ids = f.ids[rev,.]
                        }
                        for (scope=1; scope<=3; scope++) {
                            for (c=1; c<=2; c++) {
                                supplied = J(0,K,.)
                                if (c==2) supplied = f.masks \ J(1,K,0) \ f.masks
                                a = jw_run(M,f,supplied,methods[method],scopes[scope],4096)
                                b = jw020_run(M,f,supplied,methods[method],scopes[scope],4096)
                                jw_test_core_fit(a,b)
                                jw_test_core_equal(a.A*M.V*a.A',b.A*M.V*b.A')
                                calls++
                            }
                        }
                    }
                }
            }
            /* Compare supplied candidate handling, including repeated baseline. */
            f = jw_candidates(M,"focused","all",empty,empty,256)
            supplied = (f.owners,f.masks) \ (M.groups,J(G,K,0))
            oldf = jw020_candidates(M,"focused","all",empty,supplied,256)
            f = jw_candidates(M,"focused","all",empty,supplied,256)
            jw_test_core_equal(f.masks,oldf.masks)
            jw_test_core_equal(f.owners,oldf.owners)
            jw_test_core_equal(f.ids,oldf.ids)
        }
    }
    printf("Core optimization source regression: %g complete-fit comparisons PASS\n",calls)
}

void jw_test_core_zero_groups()
{
    struct jw_model scalar M
    struct jw_family scalar F
    struct jw_fit scalar R
    M.groups = J(0,1,.)
    M.z = J(0,1,.)
    M.p = 0
    M.k = 0
    R = jw_run(M,F,J(0,0,.),"focused","cohort",4096)
}

jw_test_core_modes()
end
capture noisily mata: jw_test_core_zero_groups()
assert _rc == 3498
display as result "Core optimization: zero-cohort model correctly rejected."
