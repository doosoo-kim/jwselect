*! jwselect estat 0.9.1 16nov2026
program define jwselect_estat, rclass
    version 16.1
    if "`e(cmd)'" != "jwselect" error 301
    gettoken subcmd rest : 0, parse(" ,")
    if "`subcmd'" == "" local subcmd design
    if "`subcmd'" == "pretrend" {
        _jwselect_pretrend `rest'
        return add
        exit
    }
    if "`subcmd'" == "compare" {
        _jwselect_compare `rest'
        return add
        exit
    }
    if trim(`"`rest'"') != "" error 198
    if "`subcmd'" == "design" {
        matlist e(baselines), format(%12.0g)
        if e(scalar_target)==1 {
            di as txt "One prespecified population target: `e(target_label)'"
            di as txt "The postmap row is a representative label, not a disaggregated target cell."
        }
        matlist e(postmap), format(%12.0g)
        if e(K)>0 matlist e(hmap), format(%12.0g)
        local copyfrom baselines postmap
        if e(K)>0 local copyfrom `copyfrom' hmap
        if "`e(restriction_scale)'" == "conditional_index" & e(K)>0 {
            di as txt "Selectable date blocks: cohort, time, number of independent index coefficients."
            matlist e(blockmap), format(%12.0g)
            di as txt "Coefficient components use the raw covariate basis; component 0 is the intercept."
            di as txt "Covariate components 1,... follow: `e(index_covariates)'"
            matlist e(index_terms), format(%12.0g)
            local copyfrom `copyfrom' blockmap index_terms
        }
        local copyto `copyfrom' 
    }
    else if "`subcmd'" == "selection" {
        if "`e(restriction_scale)'" == "conditional_index" {
            di as txt "Mask columns index coefficients in e(hmap); every date block moves together."
        }
        matlist e(candidates), format(%10.0g)
        if e(K)>0 matlist e(selection), format(%10.0g)
        local copyfrom candidates
        if e(K)>0 local copyfrom `copyfrom' selection
        local matrices : e(matrices)
        if strpos(" `matrices' ", " joint_candidates ") {
            matlist e(joint_candidates), format(%10.0g)
            local copyfrom `copyfrom' joint_candidates
        }
        local copyto `copyfrom' 
    }
    else if "`subcmd'" == "scores" {
        matlist e(scores), format(%10.5g)
        local copyfrom scores
        local copyto scores
    }
    else if "`subcmd'" == "covariance" {
        di as txt "Fixed-selection joint covariance; for selection-adjusted inference use e(ci) and e(pvalue)."
        matlist e(V_fixed), format(%11.5g)
        local copyfrom V_fixed
        local copyto V_fixed
    }
    else if "`subcmd'" == "inference" {
        di as txt "Inference method: `e(inference)'"
        tempname shownci shownpv
        matrix `shownci' = e(ci)
        matrix `shownpv' = e(pvalue)
        matlist `shownci', format(%12.5g)
        matlist `shownpv', format(%12.5g)
        if "`e(wbicmethod)'"=="mc" {
            forvalues i=1/`=rowsof(e(mc_diagnostics))' {
                if el(e(mc_diagnostics),`i',14)==8 ///
                    di as txt "Target row `i': CI [-inf,+inf] (whole real line; arithmetic or regularity exception)."
            }
        }
        local copyfrom ci pvalue
        local copyto ci pvalue
        if "`e(inference)'" == "reselect" {
            di as txt "Numerical search diagnostics; global extrema are not certified."
            matlist e(inference_diagnostics), format(%12.5g)
            local copyfrom `copyfrom' inference_diagnostics
            local copyto `copyto' inference_diagnostics
            local matrices : e(matrices)
            if strpos(" `matrices' ", " inference_extrema ") {
                local copyfrom `copyfrom' inference_extrema
                local copyto `copyto' inference_extrema
            }
        }
        if "`e(inference)'" == "wbic" {
            local matrices : e(matrices)
            if strpos(" `matrices' ", " wbic_bb_grid_diagnostics ") {
                di as txt "Berger-Boos fallback. Route 1 evaluated the grid. Route 2 is an exact normal calculation."
                di as txt "Route 8 uses the dimension-independent Mata bound. See help jwselect for the route codes."
                matlist e(wbic_bb_grid_diagnostics), format(%12.5g)
                local copyfrom `copyfrom' wbic_bb_grid_diagnostics
                local copyto `copyto' wbic_bb_grid_diagnostics
            }
            di as txt "WBIC calibration: `e(wbicmethod)'"
            if "`e(wbicmethod)'"=="mc" {
                di as txt "Mode: `e(mcmode)'; certified=1 is a continuum-calibration diagnostic or exact no-restriction normal inference."
                di as txt "The label is not a formal outward-rounded floating-point certificate; see the documented numerical limits."
                di as txt "Status: 0 calibrated grid; 1 skipped; 2 normal; 7 analytic terminal; 8 whole real line."
                di as txt "rule_fallback=1 means the baseline-containing safeguard of the calibrated rule."
                matlist e(mc_diagnostics), format(%12.5g)
                local copyfrom `copyfrom' mc_diagnostics
                local copyto `copyto' mc_diagnostics
                local matrices : e(matrices)
                if strpos(" `matrices' ", " mc_terminal_diagnostics ") {
                    di as txt "Terminal normal score k_safe remains finite even when its probability rounds to one."
                    di as txt "unbounded=1 explicitly denotes [-inf,+inf]; its numeric CI endpoints are stored missing."
                    matlist e(mc_terminal_diagnostics), format(%12.5g)
                    local copyfrom `copyfrom' mc_terminal_diagnostics
                    local copyto `copyto' mc_terminal_diagnostics
                }
                if strpos(" `matrices' ", " mc_profile_diagnostics ") {
                    di as txt "Conditional-score calibration. Inspect pstar, risk_upper and risk_target for each requested output."
                    matlist e(mc_profile_diagnostics), format(%12.5g)
                    local copyfrom `copyfrom' mc_profile_diagnostics
                    local copyto `copyto' mc_profile_diagnostics
                }
                if strpos(" `matrices' ", " mc_calibration ") {
                    local copyfrom `copyfrom' mc_calibration
                    local copyto `copyto' mc_calibration
                }
                if strpos(" `matrices' ", " bb_diagnostics ") {
                    di as txt "Fallback BB bracket: accuracy_certified=0 means its requested precision was not reached."
                    matlist e(bb_diagnostics), format(%12.5g)
                    local copyfrom `copyfrom' bb_diagnostics
                    local copyto `copyto' bb_diagnostics
                }
            }
            else if "`e(wbicmethod)'"=="bootstrap" {
                di as txt "Centered reselection bootstrap-t: signed quantiles are in standard-error units."
                di as txt "Gaussian draws: `e(wbicengine)' engine. CIs are equal-tailed."
                di as txt "Finite draws approximate pointwise asymptotic quantiles; no higher-order refinement is established."
                di as txt "Adjustment is the observed selected estimate minus its baseline estimate."
                di as txt "Status: 0 computed, 2 exact normal with no restrictions, 9 skipped."
                matlist e(wbic_pointwise_diagnostics), format(%12.5g)
                local copyfrom `copyfrom' wbic_pointwise_diagnostics
                local copyto `copyto' wbic_pointwise_diagnostics
            }
            di as txt "CIs are marginal, not simultaneous; estimated inputs require asymptotic conditions."
            di as txt "e(V) and se_fixed are fixed-selection diagnostics, not adjusted inference."
        }
    }
    else if "`subcmd'" == "cells" {
        if e(scalar_target)==1 {
            di as err "This source supplies one prespecified population target, with no disaggregated selected cells"
            exit 198
        }
        tempname celltab
        matrix `celltab' = e(postmap)
        mata: st_matrix("`celltab'",(st_matrix("`celltab'"),st_matrix("e(b_cells)")',sqrt(diagonal(st_matrix("e(V_cells)")))))
        matrix colnames `celltab' = cohort time estimate se_fixed
        di as txt "Cell estimates use the reported selection rule; SEs condition on selection."
        matlist `celltab', format(%12.5g)
        return matrix table = `celltab'
        local copyfrom b_cells V_cells
        local copyto b V_fixed
    }
    else {
        di as err "estat subcommands: design, selection, scores, covariance, cells, inference, pretrend, compare"
        exit 198
    }
    // return matrix requires an ordinary matrix name, not an e() expression.
    // Copy first so returning diagnostics leaves the active estimation intact.
    local ncopy : word count `copyfrom'
    forvalues i=1/`ncopy' {
        local source : word `i' of `copyfrom'
        local destination : word `i' of `copyto'
        tempname result
        matrix `result' = e(`source')
        return matrix `destination' = `result'
    }
end

program define _jwselect_compare, rclass
    version 16.1
    syntax [, CELLS OUTput(numlist integer min=1 max=1 >0)]
    if "`e(cmd)'" != "jwselect" error 301
    if e(scalar_target)==1 {
        di as err "estat compare requires disaggregated cohort cells; this source supplies one population target"
        exit 198
    }
    if "`cells'"!="" & "`output'"!="" {
        di as err "cells and output() are mutually exclusive"
        exit 198
    }
    if "`output'"=="" local output 0
    if `output'>colsof(e(b)) {
        di as err "output() is outside e(output_map)"
        exit 198
    }
    quietly _jwselect_load
    capture mata: assert(findexternal("jw_compare_stata()")!=NULL)
    if _rc {
        quietly findfile _jwselect_compare.mata
        quietly do `"`r(fn)'"'
    }
    tempname table Vbaseline Vown Vjoint Abaseline Aown Ajoint mask map
    mata: jw_compare_stata()
    di as txt "Matched baseline, own-cohort MD, and joint MD estimates"
    di as txt "Own and joint MD use the same selected configuration; targets and weights are held fixed."
    di as txt "SEs are fixed-selection diagnostics; this table supplies no selection-adjusted inference."
    if "`e(selection_method)'"=="focused" {
        if "`e(focus)'"=="aggregate" di as txt "All rows retain the configuration selected for the aggregate target."
        else {
            di as txt "Primary joint comparisons use the union of the primary cohort selections."
            if "`e(focus)'"=="both" di as txt "The independently aggregate-focused row retains its own selected configuration."
        }
        di as txt "The comparison does not extend focused inference guarantees to a different target or estimator."
    }
    matlist `table', format(%12.5g)
    return matrix table = `table'
    return matrix V_baseline = `Vbaseline'
    return matrix V_own = `Vown'
    return matrix V_joint = `Vjoint'
    return matrix loading_baseline = `Abaseline'
    return matrix loading_own = `Aown'
    return matrix loading_joint = `Ajoint'
    if e(K)>0 return matrix selection = `mask'
    if "`cells'"!="" return matrix postmap = `map'
    else return matrix output_map = `map'
    return scalar output = `output'
    return scalar cells = ("`cells'"!="")
    return local inference "fixed_selection_diagnostics"
end

program define _jwselect_pretrend, rclass
    version 16.1
    syntax [, SELECTED PREMatrix(name) OUTput(numlist integer min=1 max=1 >0)]
    if "`e(cmd)'" != "jwselect" error 301
    if "`selected'" != "" & "`prematrix'" != "" {
        di as err "selected and prematrix() are mutually exclusive"
        exit 198
    }
    if "`output'" != "" {
        if `output' < 1 | `output' >= . | "`selected'" == "" {
            di as err "output() requires selected and a positive output index"
            exit 198
        }
    }
    else local output 0
    local testset all
    if "`prematrix'" != "" {
        confirm matrix `prematrix'
        local testset custom
    }
    if "`selected'" != "" {
        local testset selected
        if `output' == 0 & "`e(selection_method)'" == "focused" {
            di as err "selected after focused selection requires output(); see e(output_map)"
            exit 198
        }
        if `output' == 0 local output 1
        if `output' > colsof(e(b)) {
            di as err "output() is outside e(output_map)"
            exit 198
        }
    }
    quietly _jwselect_load
    tempname testtable testmask testperiods
    mata: jw_pretrend_stata()
    matrix colnames `testtable' = df Wald J pvalue
    matrix rownames `testtable' = joint
    di as txt "Joint pretrend test: `testset' optional restrictions"
    di as txt "Wald and J use the same compatible covariance-weighted MD system."
    matlist `testtable', format(%12.6g)
    if `testtable'[1,1] == 0 di as txt "No optional restrictions; there is no test and the p-value is missing."
    if "`selected'" != "" {
        di as txt "Selected restrictions: the chi-squared p-value is a postselection diagnostic."
    }
    else di as txt "Chi-squared inference requires a prespecified test set and the maintained model."
    di as txt "This test does not test the maintained baseline identifying restrictions."
    return scalar chi2 = `testtable'[1,2]
    return scalar J = `testtable'[1,3]
    return scalar df = `testtable'[1,1]
    return scalar p = `testtable'[1,4]
    return scalar postselection = ("`selected'" != "")
    return scalar output = `output'
    return local test_set `testset'
    if e(K)>0 return matrix mask = `testmask'
    if `testtable'[1,1]>0 {
        matrix colnames `testperiods' = cohort time
        return matrix preperiods = `testperiods'
    }
    return matrix table = `testtable'
end
