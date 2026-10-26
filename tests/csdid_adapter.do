*! jwselect 0.5.0: CS-DR transport, validation, and state-preservation tests
*! Run from package root in a fresh test session. Does not need csdid installed.
version 16.1
set more off
adopath ++ "."

capture mata: mata drop jw_csf_fixture()
mata:
void jw_csf_fixture(string scalar bname, string scalar vname, string scalar aname,
    real scalar anti, real scalar permuted)
{
    real matrix csA, csV
    real rowvector csB, csOrd
    real scalar csJ, csI, csG, csT, csBase
    string matrix csNames
    csA = J(16,10,0)
    csB = J(1,14,0)
    csV = I(14)*.02 + J(14,14,.003)
    csNames = J(14,2,"")
    csJ = 0
    for (csI=1; csI<=16; csI++) {
        csG = csI<=8 ? 4 : 6
        csT = mod(csI-1,8)+1
        csBase = csG-anti-1
        csA[csI,.] = (csG,csT,csT-csG,0,.,20,20,20,20,csBase)
        if (csT!=csBase) {
            csJ++
            csB[csJ] = csT>=csG ? .4+.03*csG+.02*csT : .015*(csT-csBase)
            csA[csI,4] = csB[csJ]
            csA[csI,5] = sqrt(csV[csJ,csJ])
            csNames[csJ,2] = sprintf("g%g___%g_%g",csG,csT,csBase)
        }
    }
    if (permuted) {
        csOrd = (14..1)
        csB = csB[,csOrd]
        csV = csV[csOrd,csOrd]
        csNames = csNames[csOrd',.]
        csA = csA[(16..1)',.]
    }
    st_matrix(bname,csB)
    st_matrix(vname,csV)
    st_matrix(aname,csA)
    st_matrixcolstripe(bname,csNames)
    st_matrixcolstripe(vname,csNames)
    st_matrixrowstripe(vname,csNames)
}
end

capture program drop _jw_cs_fixture
program define _jw_cs_fixture, eclass
    version 16.1
    syntax [, PERMuted ANTI(integer 0)]
    tempname B V A
    tempvar sample
    quietly generate byte `sample' = 1
    mata: jw_csf_fixture("`B'", "`V'", "`A'", `anti', "`permuted'"!="")
    matrix colnames `A' = group time event_time att se n_treat_t n_treat_pre n_control_t n_control_pre base_time
    ereturn post `B' `V', esample(`sample') obs(480)
    ereturn matrix attgt = `A'
    ereturn local cmdline "csdid y x, ivar(id) time(time) gvar(group) method(dr) nevertreated base_period(universal) analytical pointwise anticipation(`anti')"
    ereturn local version "2.0.1"
    ereturn local method "dr"
    ereturn local control_group "nevertreated"
    ereturn local base_period "universal"
    ereturn local panel_mode "panel"
    ereturn local idvar "id"
    ereturn local timevar "time"
    ereturn local gvar "group"
    ereturn local depvar "y"
    ereturn local vce "analytical"
    ereturn scalar bstrap = 0
    ereturn scalar N_units = 60
    ereturn scalar N_time = 8
    ereturn scalar N_groups = 2
    ereturn scalar anticipation = `anti'
    ereturn scalar pscoretrim = .995
    ereturn scalar level = 95
    quietly signestimationsample id time group y x
    ereturn local cmd "csdid"
end

capture program drop _jw_cs_local
program define _jw_cs_local, eclass
    args key val
    ereturn local `key' "`val'"
end
capture program drop _jw_cs_scalar
program define _jw_cs_scalar, eclass
    args key val
    ereturn scalar `key' = `val'
end
capture program drop _jw_cs_badtable
program define _jw_cs_badtable, eclass
    args bad
    tempname A
    matrix `A' = e(attgt)
    if "`bad'" == "normalization" matrix `A'[3,4] = .1
    if "`bad'" == "duplicate" matrix `A'[2,1] = 4
    if "`bad'" == "duplicate" matrix `A'[2,2] = 1
    if "`bad'" == "duplicate" matrix `A'[2,3] = -3
    if "`bad'" == "missing" matrix `A'[1,4] = .
    if "`bad'" == "se" matrix `A'[1,5] = 0
    if "`bad'" == "base" matrix `A'[1,10] = 2
    ereturn matrix attgt = `A'
end

clear
set seed 260526
set obs 480
generate long id = ceil(_n/8)
generate byte time = mod(_n-1,8)+1
generate byte group = cond(mod(id,3)==0,0,cond(mod(id,3)==1,4,6))
generate double x = sin(id/3)
generate double y = .05*id + .2*time + x + rnormal()
generate double ordering = runiform()
sort ordering
quietly _jw_cs_fixture
tempname SRC B0 V0 A0 BC VC CC BT
quietly estimates store `SRC'
matrix `B0' = e(b)
matrix `V0' = e(V)
matrix `A0' = e(attgt)
generate byte saved_esample = e(sample)
quietly datasignature
local signature `"`r(datasignature)'"'
local rng `c(rngstate)'
quietly _jwselect_csdid
matrix `BC' = r(bgt)
matrix `VC' = r(Vgt)
matrix `CC' = r(cells)
matrix `BT' = r(upstream_table)
assert r(N_units)==60 & r(N_obs)==480 & r(N_clusters)==60
assert r(N_covariates)==1 & r(gap)==1 & r(anticipation)==0
assert "`r(model)'" == "csdr" & "`r(method)'" == "dr"
assert "`r(vce)'" == "supplied" & "`r(sourceversion)'" == "2.0.1"
assert "`r(restriction_scale)'" == "outcome_placebo"
assert "`e(cmd)'" == "csdid"
assert e(sample)==saved_esample
mata: assert(max(abs(st_matrix("e(b)")-st_matrix("`B0'")))==0)
mata: assert(max(abs(st_matrix("e(V)")-st_matrix("`V0'")))==0)
mata: assert(max(abs(st_matrix("`BC'")'-st_matrix("`A0'")[.,4]))<1e-12)
mata: assert(max(abs(st_matrix("`CC'")-st_matrix("`A0'")[.,1..2]))==0)
mata: assert(max(abs(st_matrix("`BT'")[.,3]-st_matrix("`BC'")'))==0)
mata: assert(max(abs(st_matrix("`BT'")[.,4]:^2-diagonal(st_matrix("`VC'"))))<1e-12)
mata: assert(st_matrix("`VC'")[1,9]==.003)
mata: assert(max(abs(st_matrix("`VC'")[3,.]))==0)
mata: assert(max(abs(st_matrix("`VC'")[13,.]))==0)
quietly datasignature
assert `"`r(datasignature)'"' == `"`signature'"'
assert "`c(rngstate)'" == "`rng'"

// Independent order permutations in the coefficient and cell tables must
// preserve coordinates and the full covariance, including across cohorts.
quietly _jw_cs_fixture, permuted
quietly _jwselect_csdid
matrix CS_PERM_B = r(bgt)
matrix CS_PERM_V = r(Vgt)
matrix CS_PERM_C = r(cells)
mata: assert(max(abs(st_matrix("CS_PERM_B")[,(16..1)]-st_matrix("`BC'")))<1e-12)
mata: assert(max(abs(st_matrix("CS_PERM_V")[(16..1)',(16..1)]-st_matrix("`VC'")))<1e-12)
mata: assert(max(abs(st_matrix("CS_PERM_C")[(16..1)',.]-st_matrix("`CC'")))==0)

// A nonzero anticipation bound moves the structural zero to g-a-1.
quietly _jw_cs_fixture, anti(1)
quietly _jwselect_csdid
assert r(anticipation)==1
matrix CS_ANTI_V = r(Vgt)
mata: assert(max(abs(st_matrix("CS_ANTI_V")[2,.]))==0)
mata: assert(max(abs(st_matrix("CS_ANTI_V")[12,.]))==0)

foreach specification in "version 1.82" "version 2.0.0" "version 3.0.0" ///
    "method reg" "control_group notyettreated" "base_period varying" ///
    "panel_mode pair-balanced" "vce bootstrap" "wtype iweight" ///
    "clustervar group" {
    quietly estimates restore `SRC'
    quietly _jw_cs_local `specification'
    capture noisily _jwselect_csdid
    assert _rc==198
    assert "`e(cmd)'" == "csdid"
}
quietly estimates restore `SRC'
quietly _jw_cs_scalar bstrap 1
capture noisily _jwselect_csdid
assert _rc==198
quietly estimates restore `SRC'
capture noisily _jwselect_csdid, vce(unconditional)
assert _rc==198
foreach defect in normalization duplicate missing se base {
    quietly estimates restore `SRC'
    quietly _jw_cs_badtable `defect'
    capture noisily _jwselect_csdid
    assert _rc==498
    assert "`e(cmd)'" == "csdid"
}

// Signed inputs detect changed data even if e(sample), row counts, and cells
// remain identical; the adapter must not alter the changed data on refusal.
quietly estimates restore `SRC'
quietly replace y = y+1 in 1
capture noisily _jwselect_csdid
assert _rc==459
quietly replace y = y-1 in 1
// Use a fresh signature after round-trip arithmetic rather than assume bits.
quietly _jw_cs_fixture
quietly replace x = x+.1 if time==2
capture noisily _jwselect_csdid
assert _rc==459
quietly replace x = sin(id/3)
quietly replace x = x+.1 if time==2
quietly _jw_cs_fixture
capture noisily _jwselect_csdid
assert _rc==198
quietly replace x = sin(id/3)
quietly _jw_cs_fixture

// Mislabelled aggregation estimates have incompatible names even if their
// dimensions happen to match; no positional fallback may accept them.
capture program drop _jw_cs_badnames
program define _jw_cs_badnames, eclass
    tempname B V
    matrix `B' = e(b)
    matrix `V' = e(V)
    local names
    forvalues i=1/14 {
        local names `names' agg`i'
    }
    matrix colnames `B' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'
    ereturn repost b=`B' V=`V', rename
end
quietly _jw_cs_badnames
local badnames : colnames e(b)
local firstname : word 1 of `badnames'
assert "`firstname'" == "agg1"
capture noisily _jwselect_csdid
assert _rc==498

quietly estimates drop `SRC'
program drop _jw_cs_fixture _jw_cs_local _jw_cs_scalar _jw_cs_badtable _jw_cs_badnames
mata: mata drop jw_csf_fixture()
display as result "jwselect 0.5.0 CS-DR adapter contract checks: PASS"
