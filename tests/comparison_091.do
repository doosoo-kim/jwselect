*! jwselect 0.9.1 matched-estimator comparison regression checks.
*! Run from the package root in a fresh Stata session.
version 16.1
clear all
set more off
adopath ++ "."
log using tests/comparison_091.log, text replace

// Analytic two-cohort example: V(h) has correlation .2 and both targets
// load on both leads. Own MD has a nonzero cross-cohort covariance.
matrix B = (.2,0,1,-.1,0,2)
matrix C = (3,1 \ 3,2 \ 3,3 \ 4,2 \ 4,3 \ 4,4)
matrix P = (3,1 \ 4,2)
mata:
H=(1,.2\ .2,1)
A=(.5,.3\ -.2,.4)
S=diag((.5,.6))
ZV=(S+A*H*A',A*H\ H*A',H)
L=J(6,4,0)
L[1,3]=1
L[3,1]=1
L[4,4]=1
L[6,2]=1
st_matrix("V",L*ZV*L')
end
jwselect, bmat(B) vmat(V) cells(C) nunits(500) select(manual) prematrix(P)
matrix before_b=e(b)
matrix before_V=e(V_fixed)
matrix before_sel=e(selection)
estat compare
matrix result=r(table)
matrix ownV=r(V_own)
matrix jointV=r(V_joint)
matrix ownA=r(loading_own)
matrix jointA=r(loading_joint)
mata:
R=st_matrix("result")
assert(max(abs(R[,1]-(1\2\1.5)))<1e-12)
assert(max(abs(R[,3]-(.888\2.036\1.462)))<1e-12)
assert(max(abs(R[,5]-(.93\2.08\1.505)))<1e-12)
assert(max(abs(R[,2]:^2-(.9\.768\.441)))<1e-12)
assert(max(abs(R[,4]:^2-(.5864\.6384\.31196)))<1e-12)
assert(max(abs(R[,6]:^2-(.5\.6\.275)))<1e-12)
assert(abs(st_matrix("ownV")[1,2]-.01152)<1e-12)
assert(max(abs(st_matrix("jointA")-st_matrix("e(loading)")))<1e-12)
assert(max(abs(st_matrix("jointV")-st_matrix("e(V_fixed)")))<1e-12)
assert(max(abs(st_matrix("before_b")-st_matrix("e(b)")))==0)
assert(max(abs(st_matrix("before_V")-st_matrix("e(V_fixed)")))==0)
assert(max(abs(st_matrix("before_sel")-st_matrix("e(selection)")))==0)
end
estat compare, output(3)
assert rowsof(r(table))==1
assert abs(r(table)[1,5]-1.505)<1e-12
estat compare, cells
assert rowsof(r(table))==2
assert abs(r(table)[1,3]-.888)<1e-12
assert abs(r(table)[2,5]-2.08)<1e-12
capture noisily estat compare, cells output(1)
assert _rc==198
capture noisily estat compare, output(100)
assert _rc==198

// Baseline-only selection must make all three estimates identical even
// when optional restrictions exist. An empty eligible lead set also works.
jwselect, bmat(B) vmat(V) cells(C) nunits(500) select(manual)
estat compare
mata: assert(max(abs(st_matrix("r(table)")[.,1]-st_matrix("r(table)")[.,5]))==0)
matrix B0=(0,1,0,2)
matrix C0=(3,2\3,3\4,3\4,4)
matrix V0=(0,0,0,0\0,.9,0,.048\0,0,0,0\0,.048,0,.768)
jwselect, bmat(B0) vmat(V0) cells(C0) nunits(500) select(wbic)
assert e(K)==0
estat compare
mata: assert(max(abs(st_matrix("r(table)")[.,1]-st_matrix("r(table)")[.,5]))==0)

// Cohort and joint WBIC must reproduce the fitted joint estimator without
// rerunning selection. Custom weights are retained in every aggregate.
matrix W=(3,3,.25\4,4,.75)
foreach scope in cohort joint {
    jwselect, bmat(B) vmat(V) cells(C) nunits(500) wbicscope(`scope') aggweights(W)
    estat compare
    mata: assert(max(abs(st_matrix("r(table)")[.,5]-st_matrix("e(b)")'))<1e-12)
    mata: assert(max(abs(st_matrix("r(V_joint)")-st_matrix("e(V_fixed)")))<1e-12)
    assert abs(r(table)[3,1]-1.75)<1e-12
}

// Focused reporting retains its target-specific masks. With focus(both),
// primary rows reproduce own MD; the separately selected aggregate row
// reproduces joint MD. Inference is deliberately disabled for this test.
// Large first lead makes cohort selection differ from aggregate selection.
matrix BF=B
matrix BF[1,1]=2
matrix BF[1,4]=.1
mata:
void jw_test_focused_comparison()
{
    real matrix R, OM
    real colvector EB
    real scalar i
    R=st_matrix("r(table)")
    EB=st_matrix("e(b)")'
    OM=st_matrix("e(output_map)")
    for (i=1;i<=rows(R);i++) {
        if (OM[i,1]==3 | OM[i,1]==5) assert(abs(R[i,5]-EB[i])<1e-12)
        else assert(abs(R[i,3]-EB[i])<1e-12)
    }
}
end
foreach scope in cohort aggregate both {
    jwselect, bmat(BF) vmat(V) cells(C) nunits(500) select(focused) focus(`scope') inference(none) aggweights(W)
    estat compare
    mata: jw_test_focused_comparison()
    if "`scope'"=="both" mata: assert(max(abs(st_matrix("r(selection)")[3,.]-st_matrix("r(selection)")[4,.]))==1)
    estat compare, cells
    if "`scope'"=="aggregate" mata: assert(max(abs(st_matrix("r(table)")[.,5]-st_matrix("e(b_cells)")'))<1e-12)
    else mata: assert(max(abs(st_matrix("r(table)")[.,3]-st_matrix("e(b_cells)")'))<1e-12)
}

// Full nonlinear index blocks: two independent coefficients per date.
// Each own-cohort projection must retain its entire two-column block.
mata:
H=I(4)
A=(.5,.2,.1,-.1\-.2,.3,.4,.2)
ZV=(diag((.5,.6))+A*A',A\A',H)
T=(I(2),J(2,4,0))
U=J(2,4,1)
own=joint=J(0,0,.)
jw_compare_loadings(ZV,T,U,(3,3\4,4),(3,1\3,1\4,2\4,2),own,joint)
assert(max(abs(joint-(I(2),-A)))<1e-12)
assert(max(abs(own-(I(2),(-.5,-.2,0,0\0,0,-.4,-.2))))<1e-12)
assert(max(abs(joint*ZV*joint'-diag((.5,.6))))<1e-12)
end

// Several cells per cohort and unequal within-cohort target weights.
// The aggregate must be formed from the cell loadings, not a weighted
// average of equally weighted cohort summaries.
do tests/raw_fixture.do
matrix WM=(3,3,.05\3,4,.10\3,5,.15\3,6,.20\4,4,.10\4,5,.10\4,6,.30)
jwselect, bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500) aggweights(WM)
estat compare
matrix multi=r(table)
matrix multi_own_V=r(V_own)
matrix multi_joint_V=r(V_joint)
estat compare, cells
mata:
w=st_matrix("e(weights)")
R=st_matrix("r(table)")
M=st_matrix("multi")
assert(abs(w*R[,1]-M[3,1])<1e-12)
assert(abs(w*R[,3]-M[3,3])<1e-12)
assert(abs(w*R[,5]-M[3,5])<1e-12)
assert(abs(w*st_matrix("r(V_own)")*w'-st_matrix("multi_own_V")[3,3])<1e-12)
assert(abs(w*st_matrix("r(V_joint)")*w'-st_matrix("multi_joint_V")[3,3])<1e-12)
end

// Scalar population sources do not expose a cohort-cell decomposition.
program define jw_test_scalar_flag, eclass
    ereturn scalar scalar_target=1
end
jw_test_scalar_flag
capture noisily estat compare
assert _rc==198
display "PASS 0.9.1 matched baseline/own/joint comparison, covariance, selections, and state preservation"
log close
