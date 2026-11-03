*! jwselect estat 0.8.8 03nov2026
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
    if trim(`"`rest'"') != "" error 198
    if "`subcmd'" == "design" {
        matlist e(baselines), format(%12.0g)
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
                di as txt "Mode: `e(mcmode)'; certified=1 means a continuum certificate or exact no-restriction normal inference."
                di as txt "Status: 0 calibrated; 1 skipped; 2 normal; 5/6 protected BB fallback."
                di as txt "rule_fallback=1 means the baseline-containing safeguard of the calibrated rule."
                matlist e(mc_diagnostics), format(%12.5g)
                local copyfrom `copyfrom' mc_diagnostics
                local copyto `copyto' mc_diagnostics
                local matrices : e(matrices)
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
        di as err "estat subcommands: design, selection, scores, covariance, cells, inference, pretrend"
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
