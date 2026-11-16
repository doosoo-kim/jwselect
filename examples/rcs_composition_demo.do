// Run after installing jwselect 0.9.1. No external first-stage package needed.
version 17
clear all
set more off
set seed 482026
set obs 6000
generate double x=rnormal()
generate int year=2001+floor(4*runiform())
// Both the treated share and its covariate distribution change with date.
generate double p=invlogit(-.3+.25*x+.2*(year-2000))
generate int first_treat=cond(runiform()<p,2004,0)
generate double y=1+.6*x+.3*(year-2000)+.8*(first_treat==2004) ///
    +2*(first_treat==2004 & year==2004)+rnormal()
jwselect_rcs y x, tvar(year) gvar(first_treat) cohort(2004) ///
    target(2004) baseline(2003) preperiods(2001 2002)
estimates store source
jwselect, select(wbic) candidates(all)
estimates restore source
jwselect, select(focused) candidates(all) inference(none)
// Restore source before requesting another selection/inference procedure.
// Refit jwselect_rcs for another target date: its standardization changes.
