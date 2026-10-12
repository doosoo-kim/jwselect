*! jwselect estat 0.2.0
program define jwselect_estat, rclass
    version 16.1
    if "`e(cmd)'" != "jwselect" error 301
    gettoken subcmd rest : 0, parse(" ,")
    if "`subcmd'" == "" local subcmd design
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
        di as txt "Cell estimates use each cohort's selected model; SEs condition on selection."
        matlist `celltab', format(%12.5g)
        return matrix table = `celltab'
        local copyfrom b_cells V_cells
        local copyto b V_fixed
    }
    else {
        di as err "estat subcommands: design, selection, scores, covariance, cells, inference"
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
