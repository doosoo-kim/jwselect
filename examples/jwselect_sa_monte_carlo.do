*! Sun-Abraham Monte Carlo for jwselect 0.7.0, draft v46_14, 30oct2026
/*
  RUN IN A FRESH STATA SESSION: this do-file clears the current data/results.
  Requires Stata 16.1+ and already installed eventstudyinteract, reghdfe,
  ftools, and avar. Installs nothing. See the companion README for details.

  do "<package>/examples/jwselect_sa_monte_carlo.do" [mode] [outer_reps] ///
       [N] [package_directory] [output_parent]

  Example from the package root:
  do examples/jwselect_sa_monte_carlo.do smoke 2 600 "." "."
  do examples/jwselect_sa_monte_carlo.do pilot 100 600 "." "."
  do examples/jwselect_sa_monte_carlo.do full 1000 1000 "." "."

  With no arguments: full, 1000 replications per DGP, N=1000. This is slow.
  Output is placed in a NEW timestamped subdirectory of output_parent.
  The smoke profile tests execution; it cannot establish nominal coverage.
  Edit profile budgets below for numerical sensitivity runs. Keep outer seeds
  unchanged to compare numerical settings using the same simulated panels.

  Target: equal-cohort ATT averaged over event times 0,1,2, tau=.655.
  Five methods: baseline; all_pre; wbic_joint_normal; focused_normal;
  focused_reselect. The last two have exactly the same point estimator.
  WBIC normal inference is pointwise-oracle inference, not uniform inference
  under local misspecification. jwselect0.7.0 has NO WBIC reselect interface.
  focused_normal is a deliberately naive fixed-selection benchmark.
  focused_reselect uses e(ci), NOT the fixed-selection covariance.
*/
version 16.1
args mode outer_reps nunits package_directory output_parent
if "`mode'"=="" local mode full
local mode=lower("`mode'")
if !inlist("`mode'","smoke","pilot","full") {
    display as error "mode must be smoke, pilot, or full"
    exit 198
}
if "`mode'"=="smoke" {
    local default_reps 2
    local default_n 600
    local search_starts 2
    local search_evals 200
    local search_tol .05
    local search_draws 199
    local check_draws 999
}
if "`mode'"=="pilot" {
    local default_reps 100
    local default_n 600
    local search_starts 4
    local search_evals 800
    local search_tol .02
    local search_draws 999
    local check_draws 9999
}
if "`mode'"=="full" {
    local default_reps 1000
    local default_n 1000
    local search_starts 8
    local search_evals 6400
    local search_tol .01
    local search_draws 9999
    local check_draws 49999
}
if "`outer_reps'"=="" local outer_reps `default_reps'
if "`nunits'"=="" local nunits `default_n'
capture confirm integer number `outer_reps'
if _rc exit 198
capture confirm integer number `nunits'
if _rc exit 198
if `outer_reps'<1 | `outer_reps'>100000 | `nunits'<100 | mod(`nunits',10)!=0 {
    display as error "Use 1<=outer_reps<=100000 and N>=100 divisible by 10."
    exit 198
}
if `"`package_directory'"'=="" local package_directory "."
if `"`output_parent'"'=="" local output_parent "."
local initial_directory `"`c(pwd)'"'
capture cd `"`package_directory'"'
if _rc {
    display as error "package_directory must be an existing jwselect0.7.0 folder."
    exit 170
}
local package_directory `"`c(pwd)'"'
quietly cd `"`initial_directory'"'
capture cd `"`output_parent'"'
if _rc {
    display as error "output_parent must already exist."
    exit 170
}
local output_parent `"`c(pwd)'"'
quietly cd `"`initial_directory'"'
capture confirm file `"`package_directory'/jwselect.ado"'
if _rc {
    display as error "Cannot find jwselect.ado in package_directory."
    exit 601
}
tempname headerfile
file open `headerfile' using `"`package_directory'/jwselect.ado"', read text
file read `headerfile' package_header
file close `headerfile'
if !strpos(`"`package_header'"',"jwselect 0.7.0") {
    display as error "This experiment requires the supplied jwselect0.7.0 package."
    exit 459
}
clear all
set more off
set rng mt64
adopath ++ `"`package_directory'"'
foreach command in jwselect jwselect_sa eventstudyinteract reghdfe ftools avar {
    capture which `command'
    if _rc {
        display as error "Missing dependency: `command'. See the package installation guide."
        exit 111
    }
}
local stamp=subinstr("`c(current_date)'"," ","",.)+"_"+subinstr("`c(current_time)'",":","",.)
local created 0
forvalues suffix=1/999 {
    if !`created' {
        local outdir `"`output_parent'/sa_mc_`stamp'_`suffix'"'
        capture mkdir `"`outdir'"'
        if !_rc local created 1
    }
}
if !`created' {
    display as error "Could not create a unique simulation directory in output_parent."
    exit 603
}
capture log close SA_MC_LOG
log using `"`outdir'/run.log"', text name(SA_MC_LOG)
display as result "Sun-Abraham Monte Carlo output: `outdir'"
display as text "mode=`mode'; outer reps per DGP=`outer_reps'; N=`nunits'"
foreach command in jwselect jwselect_sa eventstudyinteract reghdfe ftools avar {
    which `command'
}
local level 95
local gamma .01
local mcerr .001
local truth .655
local data_seed_base 614729
local inference_seed_base 371829
local normalcrit=invnormal(1-(1-`level'/100)/2)
tempname settings
file open `settings' using `"`outdir'/settings.txt"', write text
file write `settings' "jwselect 0.7.0 Sun-Abraham Monte Carlo; draft v46_14" _n
file write `settings' "Created: `c(current_date)' `c(current_time)'" _n
file write `settings' "Stata: `c(stata_version)' `c(flavor)'" _n
file write `settings' "Package directory: `package_directory'" _n
file write `settings' "mode=`mode'; outer_reps_per_DGP=`outer_reps'; N=`nunits'; T=8" _n
file write `settings' "target=.655; event window=0:2; equalcohort weights=1/2,1/2" _n
file write `settings' "DGPs: valid(delta=0), local(delta=3/sqrt(N)), fixed(delta=.5)" _n
file write `settings' "candidates=all; K=6; joint candidates=64; focus=aggregate" _n
file write `settings' "level=`level'; gamma=`gamma'; mcerr=`mcerr'" _n
file write `settings' "optstarts=`search_starts'; optmaxeval=`search_evals'; opttol=`search_tol'" _n
file write `settings' "optreps=`search_draws'; reps=`check_draws'" _n
file write `settings' "rng=mt64; data_seed=`data_seed_base'+replication (same innovations across DGPs)" _n
file write `settings' "inference_seed=`inference_seed_base'+replication (separate from DGP seed)" _n
file write `settings' "All attempted replications are retained; no retries or silent exclusions." _n
file write `settings' "Coverage conditional on a returned valid interval is reported with its denominator." _n
file write `settings' "covered_per_attempt treats failed/unavailable intervals as failures to deliver coverage." _n
file write `settings' "Numerical convergence is not certification of global nuisance extrema." _n
file write `settings' "Fixed-invalid is a stress test outside the local-sequence coverage theorem." _n
file close `settings'

// The fitted source is validated before using it. These assertion failures
// are recorded as failed replications rather than terminating the experiment.
capture program drop _jwsa_mc_source
program define _jwsa_mc_source, eclass
    quietly jwselect_sa y, ivar(id) tvar(time) gvar(group)
    assert "`e(cmd)'"=="jwselect_sa"
    assert "`e(version)'"=="0.7.0"
    assert "`e(sa_contract)'"=="0.6.0-sa1"
end

capture program drop _jwsa_mc_point
program define _jwsa_mc_point, eclass
    syntax, Method(string) Level(real)
    local options "select(manual)"
    if "`method'"=="all_pre" local options "select(manual) prematrix(SA_MC_PRE)"
    if "`method'"=="wbic_joint_normal" local options "select(wbic) wbicscope(joint) candidates(all)"
    if "`method'"=="focused" local options "select(focused) focus(aggregate) candidates(all) inference(none)"
    quietly jwselect, `options' postwindow(0 2) aggregate(equalcohort) level(`level')
    assert "`e(version)'"=="0.7.0"
    assert "`e(build)'"=="0.7.0-index1"
    assert "`e(source_estimator)'"=="eventstudyinteract"
    assert e(G)==2 & e(K)==6 & e(P)==6
    assert colsof(e(b))==3 & rowsof(e(output_map))==3
    local expected_kind=cond("`method'"=="focused",3,4)
    assert e(output_map)[3,1]==`expected_kind' & e(output_map)[3,2]==0
    assert !missing(e(b)[1,3],e(V_fixed)[3,3]) & e(V_fixed)[3,3]>0
end

tempname RESULTS DIAGS SOURCE
postfile `RESULTS' str12 design long replication str24 method ///
    long N long data_seed long inference_seed double drift double truth ///
    double estimate double se_fixed double ci_lo double ci_hi ///
    int source_rc int point_rc int inference_rc byte selected_n ///
    double selected_mask double nactive double evaluations double starts_used ///
    double converged_starts double budget_exhausted double max_validation_gap ///
    double max_radius_fraction double elapsed_seconds ///
    using `"`outdir'/sa_mc_draws.dta"', every(1)
postfile `DIAGS' str12 design long replication byte output byte objective ///
    double nactive double search_reps double validation_reps double evaluations ///
    double starts_used double converged_starts double budget_exhausted ///
    double search_value double validated_value double validation_gap ///
    double radius_fraction using `"`outdir'/inference_diagnostics.dta"', every(1)
matrix SA_MC_PRE=(4,1 \ 4,2 \ 6,1 \ 6,2 \ 6,3 \ 6,4)
// Mask bit order is (4,1),(4,2),(6,1),(6,2),(6,3),(6,4).
timer clear 91
timer on 91
foreach design in valid local fixed {
    local drift 0
    if "`design'"=="local" local drift=3/sqrt(`nunits')
    if "`design'"=="fixed" local drift .5
    forvalues replication=1/`outer_reps' {
        clear
        local data_seed=`data_seed_base'+`replication'
        local inference_seed=`inference_seed_base'+`replication'
        quietly set seed `data_seed'
        quietly set obs `nunits'
        generate long id=_n
        generate byte group=cond(mod(id,10)<4,0,cond(mod(id,10)<8,4,6))
        generate double alpha=.4*rnormal()+.25*(group==4)+.5*(group==6)
        quietly expand 8
        bysort id: generate byte time=_n
        generate double innovation=rnormal()
        bysort id (time): generate double error=innovation if _n==1
        by id (time): replace error=.5*error[_n-1]+sqrt(.75)*innovation if _n>1
        generate double y=alpha+.12*time+(.5+.15*(group==6)+.08*(time-group))* ///
            (group>0 & time>=group)+.8*error+`drift'*(group>0 & time<=2)
        capture quietly _jwsa_mc_source
        local source_rc=_rc
        if !`source_rc' quietly estimates store `SOURCE'
        if `source_rc' display as error "Source failed: `design' replication `replication', rc=`source_rc'"
        foreach method in baseline all_pre wbic_joint_normal focused {
            timer clear 92
            timer on 92
            local point_rc `source_rc'
            local inference_rc `source_rc'
            foreach item in estimate se_fixed ci_lo ci_hi selected_n selected_mask nactive evaluations starts_used converged_starts budget_exhausted max_validation_gap max_radius_fraction {
                local `item' .
            }
            if !`source_rc' {
                quietly estimates restore `SOURCE'
                capture quietly _jwsa_mc_point, method(`method') level(`level')
                local point_rc=_rc
                local inference_rc=`point_rc'
                if !`point_rc' {
                    local estimate=e(b)[1,3]
                    local se_fixed=sqrt(e(V_fixed)[3,3])
                    local selected_n 0
                    local selected_mask 0
                    forvalues k=1/6 {
                        local selected_n=`selected_n'+e(selection)[3,`k']
                        local selected_mask=`selected_mask'+e(selection)[3,`k']*2^(`k'-1)
                    }
                    if "`method'"!="focused" {
                        local ci_lo=e(ci)[3,1]
                        local ci_hi=e(ci)[3,2]
                    }
                    else {
                        local ci_lo=`estimate'-`normalcrit'*`se_fixed'
                        local ci_hi=`estimate'+`normalcrit'*`se_fixed'
                    }
                }
            }
            timer off 92
            quietly timer list 92
            local elapsed_seconds=r(t92)
            local report_method "`method'"
            if "`method'"=="focused" local report_method focused_normal
            post `RESULTS' ("`design'") (`replication') ("`report_method'") ///
                (`nunits') (`data_seed') (`inference_seed') (`drift') (`truth') ///
                (`estimate') (`se_fixed') (`ci_lo') (`ci_hi') ///
                (`source_rc') (`point_rc') (`inference_rc') (`selected_n') ///
                (`selected_mask') (`nactive') (`evaluations') (`starts_used') ///
                (`converged_starts') (`budget_exhausted') (`max_validation_gap') ///
                (`max_radius_fraction') (`elapsed_seconds')
            if "`method'"=="focused" {
                local ci_lo .
                local ci_hi .
                timer clear 92
                timer on 92
                if !`point_rc' {
                    quietly estimates restore `SOURCE'
                    capture quietly jwselect, select(focused) focus(aggregate) candidates(all) ///
                        postwindow(0 2) aggregate(equalcohort) inference(reselect) inferoutputs(3) ///
                        level(`level') gamma(`gamma') mcerr(`mcerr') seed(`inference_seed') ///
                        optstarts(`search_starts') optmaxeval(`search_evals') opttol(`search_tol') ///
                        optreps(`search_draws') reps(`check_draws')
                    local inference_rc=_rc
                    if !`inference_rc' {
                        // Point estimate must be identical; only inference changed.
                        if abs(e(b)[1,3]-`estimate')>1e-10 local inference_rc 459
                        if !`inference_rc' {
                            local ci_lo=e(ci)[3,1]
                            local ci_hi=e(ci)[3,2]
                            matrix SA_MC_D=e(inference_diagnostics)
                            local nactive 0
                            local evaluations 0
                            local starts_used 0
                            local converged_starts 0
                            local budget_exhausted 0
                            local max_validation_gap 0
                            local max_radius_fraction 0
                            forvalues d=1/`=rowsof(SA_MC_D)' {
                                if SA_MC_D[`d',1]==3 {
                                    local nactive=max(`nactive',SA_MC_D[`d',3])
                                    local evaluations=`evaluations'+SA_MC_D[`d',6]
                                    local starts_used=`starts_used'+SA_MC_D[`d',7]
                                    local converged_starts=`converged_starts'+SA_MC_D[`d',8]
                                    local budget_exhausted=`budget_exhausted'+SA_MC_D[`d',9]
                                    local max_validation_gap=max(`max_validation_gap',SA_MC_D[`d',12])
                                    local max_radius_fraction=max(`max_radius_fraction',SA_MC_D[`d',13])
                                    post `DIAGS' ("`design'") (`replication') ///
                                        (SA_MC_D[`d',1]) (SA_MC_D[`d',2]) (SA_MC_D[`d',3]) ///
                                        (SA_MC_D[`d',4]) (SA_MC_D[`d',5]) (SA_MC_D[`d',6]) ///
                                        (SA_MC_D[`d',7]) (SA_MC_D[`d',8]) (SA_MC_D[`d',9]) ///
                                        (SA_MC_D[`d',10]) (SA_MC_D[`d',11]) ///
                                        (SA_MC_D[`d',12]) (SA_MC_D[`d',13])
                                }
                            }
                        }
                    }
                }
                timer off 92
                quietly timer list 92
                local elapsed_seconds=`elapsed_seconds'+r(t92)
                post `RESULTS' ("`design'") (`replication') ("focused_reselect") ///
                    (`nunits') (`data_seed') (`inference_seed') (`drift') (`truth') ///
                    (`estimate') (`se_fixed') (`ci_lo') (`ci_hi') ///
                    (`source_rc') (`point_rc') (`inference_rc') (`selected_n') ///
                    (`selected_mask') (`nactive') (`evaluations') (`starts_used') ///
                    (`converged_starts') (`budget_exhausted') (`max_validation_gap') ///
                    (`max_radius_fraction') (`elapsed_seconds')
                if `inference_rc' display as error "Focused inference failed: `design' replication `replication', rc=`inference_rc'"
            }
            if `point_rc' display as error "Point fit failed: `design' replication `replication', `method', rc=`point_rc'"
        }
        capture estimates drop `SOURCE'
        if `replication'==1 | mod(`replication',10)==0 | `replication'==`outer_reps' {
            display as result "Completed `design': `replication' / `outer_reps'"
        }
    }
}
timer off 91
postclose `RESULTS'
postclose `DIAGS'

use `"`outdir'/inference_diagnostics.dta"', clear
export delimited using `"`outdir'/inference_diagnostics.csv"', replace
use `"`outdir'/sa_mc_draws.dta"', clear
isid design replication method
assert _N==3*`outer_reps'*5
generate byte point_ok=point_rc==0 & !missing(estimate)
generate byte interval_ok=point_ok & inference_rc==0 & !missing(ci_lo,ci_hi) & ci_lo<=ci_hi
generate double error=estimate-truth if point_ok
generate double squared_error=error^2 if point_ok
generate byte covered=(ci_lo<=truth & truth<=ci_hi) if interval_ok
generate double interval_width=ci_hi-ci_lo if interval_ok
generate byte covered_per_attempt=cond(interval_ok,covered,0)
generate byte point_failure=!point_ok
generate byte interval_failure=!interval_ok
generate byte budget_flag=budget_exhausted>0 if method=="focused_reselect" & interval_ok
generate byte mesh_flag=(converged_starts<starts_used) if method=="focused_reselect" & interval_ok & nactive>0
replace mesh_flag=0 if method=="focused_reselect" & interval_ok & nactive==0
label variable covered "Coverage among valid returned intervals"
label variable covered_per_attempt "Covered / all attempts (unavailable intervals count zero)"
label variable selected_mask "Bits: (4,1),(4,2),(6,1),(6,2),(6,3),(6,4)"
label variable budget_flag "Any of four search objectives exhausted evaluation budget"
label variable mesh_flag "At least one attempted search start did not converge"
save `"`outdir'/sa_mc_draws.dta"', replace
export delimited using `"`outdir'/sa_mc_draws.csv"', replace
collapse (count) n_attempt=replication n_estimate=error n_interval=covered ///
    (sum) n_point_failure=point_failure n_interval_failure=interval_failure ///
        n_budget_flag=budget_flag n_mesh_flag=mesh_flag ///
    (mean) bias=error mse=squared_error coverage=covered ///
        covered_per_attempt mean_width=interval_width mean_selected=selected_n ///
        mean_seconds=elapsed_seconds mean_validation_gap=max_validation_gap ///
    (sd) sd_error=error sd_squared_error=squared_error, by(design method N drift truth)
generate double rmse=sqrt(mse)
generate double bias_mcse=sd_error/sqrt(n_estimate) if n_estimate>1
generate double rmse_mcse=sd_squared_error/(2*rmse*sqrt(n_estimate)) if n_estimate>1 & rmse>0
generate double coverage_mcse=sqrt(coverage*(1-coverage)/n_interval) if n_interval>0
generate double covered_per_attempt_mcse=sqrt(covered_per_attempt*(1-covered_per_attempt)/n_attempt)
// Exact binomial confidence limits remain informative when coverage is 0 or 1.
generate double coverage_mc_lo=0 if n_interval>0
generate double coverage_mc_hi=1 if n_interval>0
generate long n_covered=round(coverage*n_interval) if n_interval>0
replace coverage_mc_lo=invibeta(n_covered,n_interval-n_covered+1,.025) if n_covered>0 & n_interval>0
replace coverage_mc_hi=invibeta(n_covered+1,n_interval-n_covered,.975) if n_covered<n_interval & n_interval>0
replace n_budget_flag=. if method!="focused_reselect"
replace n_mesh_flag=. if method!="focused_reselect"
order design method N drift truth n_attempt n_estimate n_interval n_point_failure n_interval_failure ///
    bias bias_mcse rmse rmse_mcse coverage coverage_mcse coverage_mc_lo coverage_mc_hi ///
    covered_per_attempt covered_per_attempt_mcse mean_width mean_selected n_budget_flag n_mesh_flag
format bias bias_mcse rmse rmse_mcse coverage coverage_mcse coverage_mc_lo coverage_mc_hi %9.4f
save `"`outdir'/sa_mc_summary.dta"', replace
export delimited using `"`outdir'/sa_mc_summary.csv"', replace
list design method n_estimate n_interval bias rmse coverage coverage_mcse n_interval_failure n_budget_flag, noobs sepby(design) abbreviate(24)
quietly summarize n_point_failure, meanonly
local total_point_failures=r(sum)
quietly summarize n_interval_failure, meanonly
local total_interval_failures=r(sum)
quietly summarize n_budget_flag, meanonly
local total_budget_flags=r(sum)
display as text "All attempts retained. Point failures=`total_point_failures'; interval failures=`total_interval_failures'."
display as text "Reselection intervals with at least one exhausted objective=`total_budget_flags'."
display as text "Coverage is marginal and unconditional over simulated panels, not conditional on the selected set."
display as text "Inspect failure denominators, exact Monte Carlo coverage limits, and inference_diagnostics.csv."
display as text "Successful mesh stopping does not certify global nuisance extrema or theoretical worst-case coverage."
display as result "JWSELECT_SA_MONTE_CARLO_COMPLETED"
display as result "Results: `outdir'"
timer list 91
log close SA_MC_LOG
