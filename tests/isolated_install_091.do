*! jwselect 0.9.1 isolated installation check. Does not alter the user's ado installation.
*! do tests/isolated_install_091.do "absolute/package/folder" "new/empty/test/folder"
version 17
clear all
set more off
args source destination
if `"`source'"'=="" | `"`destination'"'=="" error 198
cd `"`destination'"'
mkdir plus
mkdir personal
sysdir set PLUS `"`destination'/plus"'
sysdir set PERSONAL `"`destination'/personal"'
capture adopath - `"`source'"'
net install jwselect, from(`"`source'"')
foreach f in jwselect.ado jwselect_estat.ado _jwselect_compare.mata _jwselect_mata.do _jwselect_rcs.mata _jwselect_sa_iw_mata.do _jwselect_bb_w64.plugin {
    findfile `f'
    local installed=subinstr(`"`r(fn)'"',char(92),"/",.)
    assert strpos(lower(`"`installed'"'),lower(`"`destination'/plus/"'))==1
}
foreach h in jwselect jwselect_rcs jwselect_sa jwselect_sa_iw jwselect_worker jwselect_collect {
    findfile `h'.sthlp
}
do `"`source'/examples/comparison_demo.do"'
assert "`e(version)'"=="0.9.1"
assert "`e(build)'"=="0.9.1-20261116"
estat compare
assert colsof(r(table))==9
mata: assert(max(abs(st_matrix("r(table)")[.,5]-st_matrix("e(b)")'))<1e-12)
do `"`source'/examples/matrix_demo.do"'
do `"`source'/examples/rcs_composition_demo.do"'
assert "`e(build)'"=="0.9.1-20261116"
capture noisily estat compare
assert _rc==198
display "ISOLATED_INSTALL_091_PASS"
