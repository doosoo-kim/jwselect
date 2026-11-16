*! jwselect 0.9.1 16nov2026: composition-robust repeated-cross-section DR inputs
program define jwselect_rcs, eclass
    version 16.1
    if replay() {
        if "`e(cmd)'" != "jwselect_rcs" error 301
        ereturn display
        exit
    }
    local cmdline `"jwselect_rcs `0'"'
    syntax varlist(numeric min=1) [if] [in], TVAR(varname numeric) ///
        GVAR(varname numeric) COHORT(integer) TARGET(integer) BASELINE(integer) ///
        [PREPeriods(numlist integer sort) GAP(integer 1) PSVars(varlist numeric) ///
         ORConstant PSConstant SAMPling(string) Level(cilevel) ITERate(integer 100)]
    if "`sampling'"=="" local sampling pooled
    local sampling=lower("`sampling'")
    if !inlist("`sampling'","pooled","stratified") {
        di as err "sampling() must be pooled or stratified"
        exit 198
    }
    gettoken y xvars : varlist
    if `gap'<=0 | `cohort'<=0 | `target'<`cohort' | `baseline'>=`cohort' | ///
        `baseline'<=0 | mod(`cohort'-`baseline',`gap')!=0 | ///
        mod(`target'-`cohort',`gap')!=0 | `iterate'<1 {
        di as err "Require positive integer dates/gap, baseline < cohort <= target, compatible gap, and iterate()>0"
        exit 198
    }
    local anti=(`cohort'-`baseline')/`gap'-1
    local pre : list uniq preperiods
    foreach tt of local pre {
        if `tt'<=0 | `tt'>`baseline' | mod(`tt'-`baseline',`gap')!=0 {
            di as err "preperiods() must contain untreated dates no later than baseline(), on the specified gap"
            exit 198
        }
    }
    local dates `pre' `baseline' `target'
    local dates : list uniq dates
    numlist "`dates'", sort
    local dates `r(numlist)'
    if "`psvars'"=="" local psvars `xvars'
    local orvars `xvars'
    if "`orconstant'"!="" local orvars
    if "`psconstant'"!="" local psvars
    local covars `orvars' `psvars'
    local covars : list uniq covars
    foreach vv of local covars {
        if inlist("`vv'","`y'","`tvar'","`gvar'") {
            di as err "Covariates must differ from outcome, date, and cohort variables"
            exit 198
        }
    }
    marksample touse, novarlist
    markout `touse' `y' `tvar' `gvar' `covars'
    tempvar eligible
    quietly gen byte `eligible'=0
    foreach tt of local dates {
        quietly replace `eligible'=1 if `tvar'==`tt'
    }
    quietly replace `touse'=0 if !`eligible' | !inlist(`gvar',0,`cohort')
    quietly count if `touse'
    local n=r(N)
    if `n'<=1 error 2000
    quietly findfile _jwselect_rcs.mata
    local mfile `"`r(fn)'"'
    capture mata: assert(jwrcs_version()=="0.9.1")
    if _rc {
        capture mata: mata drop jwrcs_*()
        quietly do `"`mfile'"'
    }
    tempname B V CELLS MEANS MV DIAG JD JH JV JHM JB JTM
    mata: jwrcs_fit("`y'","`orvars'","`psvars'","`tvar'","`gvar'", ///
        "`touse'",strtoreal(tokens("`dates'")),`cohort',`target',`baseline', ///
        `iterate',"`sampling'","`B'","`V'","`CELLS'","`MEANS'","`MV'","`DIAG'")
    local names
    foreach tt of local dates {
        local names `names' g`cohort'_t`tt'
    }
    matrix colnames `B' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'
    matrix colnames `CELLS' = cohort time
    matrix colnames `MEANS' = cohort time estimate
    matrix colnames `DIAG' = cohort time N_source min_pair_p max_pair_p max_ratio max_score
    ereturn post `B' `V', esample(`touse') obs(`n') depname(`y')
    matrix `B'=e(b)
    matrix `V'=e(V)
    mata: st_matrix("`JD'",st_matrix("`B'")[,selectindex(st_matrix("`CELLS'")[,2]:==`target')])
    mata: st_matrix("`JV'",st_matrix("`V'")[(selectindex(st_matrix("`CELLS'")[,2]:==`target')\selectindex(st_matrix("`CELLS'")[,2]:<`baseline')),(selectindex(st_matrix("`CELLS'")[,2]:==`target')\selectindex(st_matrix("`CELLS'")[,2]:<`baseline'))])
    local kh=colsof(`B')-2
    if `kh'>0 {
        mata: st_matrix("`JH'",st_matrix("`B'")[,selectindex(st_matrix("`CELLS'")[,2]:<`baseline')])
        mata: st_matrix("`JHM'",st_matrix("`CELLS'")[selectindex(st_matrix("`CELLS'")[,2]:<`baseline'),])
        ereturn matrix jw_h=`JH'
        ereturn matrix jw_hmap=`JHM'
    }
    matrix `JB'=(`cohort',`baseline')
    matrix `JTM'=(`cohort',`target')
    ereturn matrix jw_d=`JD'
    ereturn matrix jw_V=`JV'
    ereturn matrix jw_bases=`JB'
    ereturn matrix jw_targetmap=`JTM'
    ereturn scalar jw_nunits=`n'
    ereturn scalar jw_K=`kh'
    ereturn local jw_target_label "ATT for cohort `cohort' at date `target'; post-treated population"
    ereturn local jw_schema "scalar_v1"
    ereturn local jw_restriction_scale "target_standardized_average"
    ereturn matrix bgt=`B'
    ereturn matrix Vgt=`V'
    ereturn matrix cells=`CELLS'
    ereturn matrix transported_means=`MEANS'
    ereturn matrix transported_V=`MV'
    ereturn matrix rcs_diagnostics=`DIAG'
    ereturn scalar N_units=`n'
    ereturn scalar N_clusters=`n'
    ereturn scalar N_covariates=`: word count `covars''
    ereturn scalar gap=`gap'
    ereturn scalar anticipation=`anti'
    ereturn scalar rcs_cohort=`cohort'
    ereturn scalar rcs_target=`target'
    ereturn scalar rcs_baseline=`baseline'
    ereturn scalar level=`level'
    ereturn scalar pscoretrim=1
    ereturn local cmd jwselect_rcs
    ereturn local cmdline `"`cmdline'"'
    ereturn local version "0.9.1"
    ereturn local rcs_contract "0.9.0-composition1"
    ereturn local rcs_tvar "`tvar'"
    ereturn local rcs_gvar "`gvar'"
    ereturn local rcs_outcome "`y'"
    ereturn local rcs_orvars "`orvars'"
    ereturn local rcs_psvars "`psvars'"
    ereturn local rcs_dates "`dates'"
    ereturn local vce robust
    ereturn local vcetype Robust
    ereturn local method "transport_dr"
    ereturn local composition "changing"
    ereturn local rcs_sampling "`sampling'"
    signestimationsample `y' `tvar' `gvar' `covars'
    di as text "Composition-robust DR repeated cross sections; target distribution: cohort `cohort', date `target'"
    di as text "Untrimmed transported scores; full nuisance and target-distribution covariance; N = `n'"
    if "`sampling'"=="stratified" di as text "Fixed date sample sizes: complete influences centered within each date"
    ereturn display, level(`level')
end
