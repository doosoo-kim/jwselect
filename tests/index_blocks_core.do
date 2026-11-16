*! jwselect 0.7.0 conditional-index core regression checks.
*! Run from the extracted package root in a fresh Stata session.
version 16.1
set more off
adopath ++ "."
quietly _jwselect_load

mata:
void jw_test_index_blocks_core()
{
    struct jw_model scalar M, Z, MD
    struct jw_family scalar C, N, U, C0
    struct jw_fit scalar F, FM, FD
    struct jw_rs_context scalar RS
    real matrix H, A, S, V, L, cells, hm, b, Vraw, cross, blocks, joint, expectedA, expectedV, table, cm, E, errors
    real colvector z, h, ix, mu, shift
    real rowvector target, mask, stat, test, adjusted, perm
    real scalar i, o, cache, W, gain, vd, delta

    /* Four outcome ATT cells followed by two two-coefficient lead blocks.
       A nonzero cross-cohort covariance makes the joint system substantive. */
    H = (.01,0,.002,0 \ 0,.02,0,0 \ .002,0,.015,0 \ 0,0,0,.025)
    A = (.5,.2,.1,-.1 \ .3,-.4,.2,.1 \ -.2,.1,.6,.3 \ .1,.1,.4,-.2)
    S = diag((.02,.03,.025,.035))
    V = (S+A*H*A',A*H \ H*A',H)
    z = (.3 \ .5 \ .2 \ .4 \ .03 \ .02 \ .45 \ -.3)
    h = z[5..8]
    cells = (3,1 \ 3,2 \ 3,3 \ 3,4 \ 4,2 \ 4,3 \ 4,4 \ 4,5)
    hm = (3,1 \ 3,1 \ 4,2 \ 4,2)
    L = J(8,8,0)
    L[1,5..6] = (1,.2)
    L[3,1] = 1
    L[4,2] = 1
    L[5,7..8] = (1,-.1)
    L[7,3] = 1
    L[8,4] = 1
    b = (L*z)'
    Vraw = L*V*L'
    cross = L*V[,5..8]
    M = jw_design_index(b,Vraw,cells,h,H,cross,hm,100,1,0,"cohort",(0,1),"equalcohort",J(0,0,.))
    assert(M.k==4 & M.p==4 & rows(M.groups)==2)
    assert(max(abs(M.z-z))<1e-14)
    assert(max(abs(M.V-V))<1e-14)
    assert(max(abs(M.hmap-hm))==0)
    blocks = jw_block_map(M.hmap)
    assert(max(abs(blocks-(3,1,2 \ 4,2,2)))==0)
    assert(max(abs(M.w-J(1,4,.25)))==0)
    C = jw_index_candidates(M,"wbic","all",J(0,0,.),J(0,0,.),2)
    assert(rows(C.masks)==4)
    assert(max(abs(C.masks-(0,0,0,0 \ 1,1,0,0 \ 0,0,0,0 \ 0,0,1,1)))==0)
    N = jw_index_candidates(M,"wbic","nested",J(0,0,.),J(0,0,.),2)
    assert(max(abs(C.masks-N.masks))==0)
    joint = jw_joint_masks(M,C,J(0,0,.),4)
    assert(rows(joint)==4 & cols(joint)==4)
    jw_block_masks_check(M.hmap,joint)

    /* Independent block-rank score and target-MD calculations, using an
       explicit inverse rather than the production Cholesky loading solver. */
    target = (.25,.25,.25,.25,0,0,0,0)
    table = J(4,8,.)
    for (i=1; i<=4; i++) {
        mask = joint[i,.]
        adjusted = target
        W = 0
        gain = 0
        if (sum(mask)>0) {
            ix = vec(selectindex(mask')):+4
            W = z[ix]'*invsym(V[ix,ix])*z[ix]
            adjusted[ix] = adjusted[ix] - (target*V[,ix]*invsym(V[ix,ix]))
            gain = target*V*(target-adjusted)'
        }
        delta = (target-adjusted)*z
        vd = adjusted*V*adjusted'
        table[i,.] = (sum(mask),W,W-sum(mask)*ln(100),adjusted*z,sqrt(vd),delta,gain,delta^2+target*V*target'-2*gain)
        stat = jw_candidate_stats(M,target,mask)
        assert(max(abs(stat-table[i,.]))<1e-11)
    }
    F = jw_run_scope(M,C,J(0,0,.),"wbic","both",4,"joint")
    ix = vec(selectindex(table[,3]:==min(table[,3])))
    assert(max(abs(F.selected[3,.]-joint[ix[1],.]))==0)
    assert(abs(F.A[3,.]*z-table[ix[1],4])<1e-12)
    assert(any(F.scores[,4]:==2) & any(F.scores[,4]:==4))
    assert(!any(F.scores[,4]:==1) & !any(F.scores[,4]:==3))
    U = jw_index_candidates(M,"manual","nested",(3,1 \ 4,2),J(0,0,.),2)
    FM = jw_run_scope(M,U,J(0,0,.),"manual","both",4,"cohort")
    expectedA = (I(4),-A)
    expectedV = S
    assert(max(abs(FM.cellA-expectedA))<1e-12)
    assert(max(abs(FM.cellA*V*FM.cellA'-expectedV))<1e-12)
    assert(abs(FM.A[3,.]*z-.206)<1e-12)
    /* Coefficient order is an adapter contract, not a date-sort assumption.
       Interleaving two date blocks must preserve the final MD estimate. */
    perm = (1,3,2,4)
    MD = jw_design_index(b,Vraw,cells,h[perm],H[perm,perm],cross[,perm],hm[perm,.],100,1,0,"cohort",(0,1),"equalcohort",J(0,0,.))
    U = jw_index_candidates(MD,"wbic","all",J(0,0,.),J(0,0,.),2)
    assert(max(abs(U.masks-(0,0,0,0 \ 1,0,1,0 \ 0,0,0,0 \ 0,1,0,1)))==0)
    U = jw_index_candidates(MD,"manual","nested",(3,1 \ 4,2),J(0,0,.),2)
    FD = jw_run_scope(MD,U,J(0,0,.),"manual","both",4,"cohort")
    assert(max(abs(FD.A*MD.z-FM.A*M.z))<1e-12)
    assert(max(abs(FD.A*MD.V*FD.A'-FM.A*M.V*FM.A'))<1e-12)
    test = jw_pretrend_test(M.z,M.V,(1,1,0,0),4)
    assert(test[1]==2)
    assert(abs(test[2]-h[1..2]'*invsym(H[1..2,1..2])*h[1..2])<1e-12)

    /* A custom family uses coefficient bits, but never splits a date. */
    cm = (3,1,1,0,0 \ 4,0,0,1,1)
    U = jw_index_candidates(M,"focused","all",J(0,0,.),cm,2)
    assert(max(abs(U.masks-C.masks))==0)
    F = jw_run_scope(M,U,joint,"focused","both",4,"cohort")
    jw_block_masks_check(M.hmap,F.selected)
    assert(rows(F.A)==4)
    assert(max(abs(F.A*M.V*F.A'-(F.A*M.V*F.A')'))<1e-12)

    /* Reselection must reselect entire date blocks on every draw. Compare
       cached and streaming simulation with full model refits at fixed
       nuisance means, keeping all cross-cohort covariance unchanged. */
    E = (0,0,0,0,0,0,0,0 \ ///
        .1,-.2,.3,-.1,.02,.3,-.1,.2 \ ///
        -.3,.2,.1,.4,-.25,.1,.4,-.2 \ ///
        .2,.1,-.4,.1,.3,-.1,-.2,.1 \ ///
        .5,-.4,.2,-.3,-.1,-.4,.2,.3)
    mu = (.05 \ -.1 \ .2 \ -.15)
    shift = J(M.p,1,0) \ mu
    errors = J(rows(E),rows(F.A),.)
    MD = M
    for (i=1; i<=rows(E); i++) {
        MD.z = E[i,.]' + shift
        FD = jw_run_scope(MD,U,joint,"focused","both",4,"cohort")
        jw_block_masks_check(M.hmap,FD.selected)
        errors[i,.] = (FD.A*MD.z-FD.targets*shift)'
    }
    for (o=1; o<=rows(F.A); o++) {
        assert(max(abs(jw_rs_errors(M,F,o,E,mu)-errors[,o]))<1e-12)
        for (cache=0; cache<=1; cache++) {
            RS = jw_rs_prepare_opts(M,F,o,E,(5,6,7,8),cache)
            assert(max(abs(jw_rs_draw_errors(RS,mu)-errors[,o]))<1e-12)
        }
    }

    /* Zero optional blocks remain a supported analytical design. */
    Z = jw_design_index((0,.4),(0,0 \ 0,.04),(3,2 \ 3,3),
        J(0,0,.),J(0,0,.),J(0,0,.),J(0,0,.),100,1,0,"cohort",(0,0),"equalcohort",J(0,0,.))
    assert(Z.k==0 & Z.p==1 & Z.z[1]==.4)
    C0 = jw_index_candidates(Z,"focused","all",J(0,0,.),J(0,0,.),1)
    assert(rows(C0.masks)==1 & cols(C0.masks)==0)
    F = jw_run_scope(Z,C0,J(0,0,.),"focused","aggregate",1,"cohort")
    assert(max(abs(F.A*Z.z:-.4))<1e-14)

    /* Export fixture data to exercise invalid inputs and Stata transport. */
    st_matrix("IB_B",b)
    st_matrix("IB_V",Vraw)
    st_matrix("IB_C",cells)
    st_matrix("IB_H",h')
    st_matrix("IB_VH",H)
    st_matrix("IB_VDH",cross)
    st_matrix("IB_HMAP",hm)
    st_matrix("IB_Z",z')
    st_matrix("IB_ZV",V)
    st_matrix("IB_ALLB",(FM.A*z)')
    st_matrix("IB_ALLV",FM.A*V*FM.A')
}
jw_test_index_blocks_core()
mata drop jw_test_index_blocks_core()
end
display as result "PASS index blocks: complete-date families, rank penalties, independent MD and zero-block case"

capture mata: jw_block_masks_check(st_matrix("IB_HMAP"),(1,0,0,0))
assert _rc==3498
capture mata: jw_block_masks_check(st_matrix("IB_HMAP"),(0,0,0,.))
assert _rc==3498
capture mata: jw_design_index(st_matrix("IB_B"),st_matrix("IB_V"),st_matrix("IB_C"),st_matrix("IB_H"),J(4,4,0),st_matrix("IB_VDH"),st_matrix("IB_HMAP"),100,1,0,"cohort",(0,1),"equalcohort",J(0,0,.))
assert _rc==3498
capture mata: jw_design_index(st_matrix("IB_B"),st_matrix("IB_V"),st_matrix("IB_C"),st_matrix("IB_H"),st_matrix("IB_VH"),J(8,4,10),st_matrix("IB_HMAP"),100,1,0,"cohort",(0,1),"equalcohort",J(0,0,.))
assert _rc==3498
capture mata: jw_design_index(st_matrix("IB_B"),st_matrix("IB_V"),st_matrix("IB_C"),st_matrix("IB_H"),st_matrix("IB_VH"),st_matrix("IB_VDH"),st_matrix("IB_HMAP"),100,1,0,"common",(0,1),"equalcohort",J(0,0,.))
assert _rc==3498
display as result "PASS index blocks: split blocks, invalid covariance, and unsupported rebasing rejected"
display as result "JWSELECT_INDEX_BLOCKS_CORE_PASS"
