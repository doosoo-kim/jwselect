*! jwselect estat 0.5.0 23oct2026
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
        local copyto `copyfrom' 
    }
    else if "`subcmd'" == "selection" {
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
        di as txt "Fixed-selection joint covariance; for focused inference use e(ci) and e(pvalue)."
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
