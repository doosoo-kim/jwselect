*! jwselect 0.9.1 16nov2026: composition-robust RCS input adapter
program define _jwselect_rcs, rclass
    version 16.1
    syntax [, VCE(string)]
    if "`vce'"=="" local vce supplied
    if "`vce'"!="supplied" {
        di as err "The RCS adapter requires vce(supplied), including nuisance and target-distribution uncertainty"
        exit 198
    }
    if "`e(cmd)'"!="jwselect_rcs" | "`e(rcs_contract)'"!="0.9.0-composition1" {
        di as err "Run jwselect_rcs first"
        exit 301
    }
    capture checkestimationsample
    if _rc {
        di as err "The repeated-cross-section source data changed; restore them or refit jwselect_rcs"
        exit 459
    }
    if e(pscoretrim)!=1 {
        di as err "Active trimming is not supported by the DR RCS construction"
        exit 198
    }
    tempname B V C T
    matrix `B'=e(bgt)
    matrix `V'=e(Vgt)
    matrix `C'=e(cells)
    local k=colsof(`B')
    matrix `T'=J(`k',8,.)
    matrix `T'[1,1]=`C'
    local z=invnormal(1-(100-e(level))/200)
    forvalues j=1/`k' {
        matrix `T'[`j',3]=`B'[1,`j']
        matrix `T'[`j',4]=sqrt(`V'[`j',`j'])
        if `V'[`j',`j']>0 {
            matrix `T'[`j',5]=`B'[1,`j']/sqrt(`V'[`j',`j'])
            matrix `T'[`j',6]=2*normal(-abs(`T'[`j',5]))
            matrix `T'[`j',7]=`B'[1,`j']-`z'*sqrt(`V'[`j',`j'])
            matrix `T'[`j',8]=`B'[1,`j']+`z'*sqrt(`V'[`j',`j'])
        }
    }
    matrix colnames `T'=cohort time estimate se statistic pvalue ll ul
    return matrix bgt=`B'
    return matrix Vgt=`V'
    return matrix cells=`C'
    return matrix upstream_table=`T'
    return scalar N_units=e(N)
    return scalar N_obs=e(N)
    return scalar N_clusters=e(N)
    return scalar N_covariates=e(N_covariates)
    return scalar gap=e(gap)
    return scalar anticipation=e(anticipation)
    return scalar df_r=.
    return scalar pscoretrim=1
    return local ivar ""
    return local tvar "`e(rcs_tvar)'"
    return local gvar "`e(rcs_gvar)'"
    return local sourcecmd `"`e(cmdline)'"'
    return local sourcevce "robust"
    return local sourceversion "0.9.1"
    return local sourceestimator "jwselect_rcs"
    return local vce "supplied"
    return local model "rcs"
    return local method "transport_dr"
    return local restriction_scale "target_standardized_average"
    return local upstream_stat "z"
end
