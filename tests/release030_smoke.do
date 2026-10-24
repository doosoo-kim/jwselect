*! Called by run_release_tests.do from the package root.
*! Test-only numerical settings are intentionally small, not publication defaults.
version 17.0
args outdir
if `"`outdir'"' == "" error 198

do tests/release030_assertions.do

do tests/raw_fixture.do
local input bmat(JW_B) vmat(JW_V) cells(JW_C) nunits(500)
local numerical reps(1999) optreps(199) optstarts(2) optmaxeval(80) seed(913)

// First-call autoload, correlated covariance, and baseline/WBIC regression guard.
jwselect, `input' select(manual)
mata: assert(jw_version()=="0.4.0-joint1")
mata: jw030_equal("e(b)","JW_EB_BASE",1e-10)
mata: jw030_equal("e(V)","JW_EV_BASE",1e-10)
jwselect, `input' select(wbic) candidates(all)
mata: jw030_equal("e(b)","JW_EB_WBIC",1e-10)
mata: jw030_equal("e(V)","JW_EV_WBIC",1e-10)
// Test complete joint covariance, not independent cohort approximations.
mata: assert(max(abs(st_matrix("JW_V")[1..6,7..12]))>0)

local old_rng `c(rngstate)'
local caseno 0
foreach family in all nested {
    foreach scope in cohort aggregate both {
        local ++caseno
        local jobdir `"`outdir'/case_`family'_`scope'"'
        local config `input' select(focused) candidates(`family') focus(`scope') `numerical'
        display as text "CASE `caseno': family=`family', focus=`scope'"

        // Serial, uncached: the comparison includes fixed-action SE and inference.
        jwselect, `config' inference(reselect) cachemb(0)
        assert `"`c(rngstate)'"' == `"`old_rng'"'
        matrix JW030_SERIAL = e(table)
        matrix JW030_DIAG = e(inference_diagnostics)
        matrix JW030_EXT = e(inference_extrema)
        matrix JW030_MAP = e(output_map)
        matrix JW030_B = e(b)
        matrix JW030_V = e(V_fixed)
        mata: assert(!hasmissing(st_matrix("JW030_SERIAL")))
        mata: assert(all(st_matrix("JW030_SERIAL")[.,3]:<=st_matrix("JW030_SERIAL")[.,4]))
        mata: assert(all(st_matrix("JW030_SERIAL")[.,5]:>=0) & all(st_matrix("JW030_SERIAL")[.,5]:<=1))
        mata: assert(all(st_matrix("JW030_DIAG")[.,6]:<=80))

        // A positive cache budget may change storage, never the statistical rule.
        jwselect, `config' inference(reselect) cachemb(64)
        mata: jw030_equal("e(table)","JW030_SERIAL",1e-10)
        mata: jw030_equal("e(inference_diagnostics)","JW030_DIAG",1e-10)
        mata: jw030_equal("e(inference_extrema)","JW030_EXT",1e-10)
        assert `"`c(rngstate)'"' == `"`old_rng'"'

        // Preparation freezes inputs/draws but does not run numerical inference.
        jwselect, `config' inference(none) prepare(`"`jobdir'"') cachemb(64)
        matrix JW030_PREP_B = e(b)
        assert `"`c(rngstate)'"' == `"`old_rng'"'
        mata: jw030_equal("e(b)","JW030_B",1e-10)
        mata: jw030_equal("e(V_fixed)","JW030_V",1e-10)
        capture noisily jwselect_collect using `"`jobdir'"'
        assert _rc != 0
        assert `"`c(rngstate)'"' == `"`old_rng'"'

        // Frozen preparations cannot be overwritten accidentally.
        capture noisily jwselect, `config' inference(none) prepare(`"`jobdir'"')
        assert _rc != 0
        mata: jw030_equal("e(b)","JW030_PREP_B",1e-10)

        // A competing session's lock must fail closed rather than duplicate work.
        if `caseno' == 1 {
            mkdir `"`jobdir'/job_1.lock"'
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1)
            assert _rc != 0
            capture noisily jwselect_collect using `"`jobdir'"'
            assert _rc != 0
            rmdir `"`jobdir'/job_1.lock"'
            // A later acquisition failure releases only our earlier lock.
            mkdir `"`jobdir'/job_2.lock"'
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1) jobids(1 2)
            assert _rc==603
            mata: assert(!direxists(st_local("jobdir")+"/job_1.lock"))
            mata: assert(direxists(st_local("jobdir")+"/job_2.lock"))
            rmdir `"`jobdir'/job_2.lock"'
            // Missing READY means incomplete preparation: no worker may proceed.
            copy `"`jobdir'/READY.dta"' `"`jobdir'/READY_test_backup.dta"'
            erase `"`jobdir'/READY.dta"'
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1)
            assert _rc != 0
            copy `"`jobdir'/READY_test_backup.dta"' `"`jobdir'/READY.dta"'
            erase `"`jobdir'/READY_test_backup.dta"'
            // Verify frozen-code integrity before allowing any numerical work.
            copy `"`jobdir'/code/_jwselect_mata.do"' `"`jobdir'/code/test_original.do"'
            tempname corrupt
            file open `corrupt' using `"`jobdir'/code/_jwselect_mata.do"', write append text
            file write `corrupt' _n "// deliberate test-only checksum mismatch" _n
            file close `corrupt'
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1)
            assert _rc != 0
            copy `"`jobdir'/code/test_original.do"' `"`jobdir'/code/_jwselect_mata.do"', replace
            erase `"`jobdir'/code/test_original.do"'
            // A crashed, unmarked partial result must be replaced safely.
            tempname partial
            file open `partial' using `"`jobdir'/job_1.bin"', write binary
            file close `partial'
            // jobids intersects all other filters and the GLOBAL shard index.
            // Of IDs 1/2/3, only 2 is output1/objective1/start2 and in shard2.
            jwselect_worker using `"`jobdir'"', worker(2) workers(2) ///
                jobids(1 2 3) outputs(1) objectives(1) starts(2) cachemb(.002)
            assert r(assigned)==1
            assert r(lock_failures)==0
            confirm file `"`jobdir'/job_2.bin.done"'
            capture confirm file `"`jobdir'/job_1.bin.done"'
            assert _rc != 0
            capture confirm file `"`jobdir'/job_3.bin.done"'
            assert _rc != 0
            // Repeating an explicit completed ID is an idempotent restart.
            jwselect_worker using `"`jobdir'"', worker(1) workers(1) jobids(2) cachemb(0)
            assert r(assigned)==1
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1) jobids(0)
            assert _rc != 0
            assert `"`c(rngstate)'"' == `"`old_rng'"'
        }

        // Exercise objective filters, then out-of-order sharding and restart.
        jwselect_worker using `"`jobdir'"', worker(1) workers(1) objectives(1) cachemb(0)
        capture noisily jwselect_collect using `"`jobdir'"'
        assert _rc != 0
        forvalues worker=4(-1)1 {
            jwselect_worker using `"`jobdir'"', worker(`worker') workers(4) cachemb(64)
            assert r(lock_failures)==0
            assert `"`c(rngstate)'"' == `"`old_rng'"'
        }
        // Repeating completed assignments must safely skip completed jobs.
        jwselect_worker using `"`jobdir'"', worker(1) workers(1) cachemb(0)
        if `caseno' == 1 {
            // A marked result from the wrong job must never be merged silently.
            copy `"`jobdir'/job_2.bin"' `"`jobdir'/job_2_test_original.bin"'
            copy `"`jobdir'/job_1.bin"' `"`jobdir'/job_2.bin"', replace
            // A result-validation error must keep its rc and release our lock.
            capture noisily jwselect_worker using `"`jobdir'"', worker(1) workers(1) jobids(2)
            assert _rc==3498
            mata: assert(!direxists(st_local("jobdir")+"/job_2.lock"))
            assert `"`c(rngstate)'"' == `"`old_rng'"'
            capture noisily jwselect_collect using `"`jobdir'"'
            assert _rc != 0
            copy `"`jobdir'/job_2_test_original.bin"' `"`jobdir'/job_2.bin"', replace
            erase `"`jobdir'/job_2_test_original.bin"'
        }
        jwselect_collect using `"`jobdir'"'
        assert r(N_jobs)>0
        matrix JW030_COLLECT = r(table)
        matrix JW030_CDIAG = r(diagnostics)
        matrix JW030_CEXT = r(extrema)
        matrix JW030_CMAP = r(output_map)
        matrix JW030_JOBS = r(jobs)
        mata: jw030_equal("JW030_COLLECT","JW030_SERIAL",1e-10)
        mata: jw030_equal("JW030_CDIAG","JW030_DIAG",1e-10)
        mata: jw030_equal("JW030_CEXT","JW030_EXT",1e-10)
        mata: jw030_equal("JW030_CMAP","JW030_MAP",1e-10)
        mata: assert(cols(st_matrix("JW030_JOBS"))==4)
        mata: assert(rows(st_matrix("JW030_JOBS"))==8*rows(st_matrix("JW030_MAP")))
        // Worker/collector are rclass; preparation's posted estimates survive.
        mata: jw030_equal("e(b)","JW030_PREP_B",1e-10)
        assert `"`c(rngstate)'"' == `"`old_rng'"'
    }
}

// Request inference for one cohort only. All point estimates/covariances remain.
local config `input' select(focused) candidates(all) focus(cohort) `numerical'
jwselect, `config' inference(reselect) cachemb(64)
matrix JW030_FULL = e(table)
matrix JW030_FULL_B = e(b)
matrix JW030_FULL_V = e(V_fixed)
jwselect, `config' inference(reselect) inferoutputs(1) cachemb(64)
matrix JW030_PART = e(table)
mata: jw030_equal("e(b)","JW030_FULL_B",1e-10)
mata: jw030_equal("e(V_fixed)","JW030_FULL_V",1e-10)
mata: assert(max(abs(st_matrix("JW030_PART")[1,.]-st_matrix("JW030_FULL")[1,.]))<1e-10)
mata: assert(all(vec(st_matrix("JW030_PART")[2..3,3..5]:>=.)))
local jobdir `"`outdir'/one_output"'
jwselect, `config' inference(none) inferoutputs(1) prepare(`"`jobdir'"') cachemb(64)
// Split explicitly by start; output index is an output row, not a cohort year.
jwselect_worker using `"`jobdir'"', worker(1) workers(1) outputs(1) starts(2) cachemb(0)
capture noisily jwselect_collect using `"`jobdir'"'
assert _rc != 0
jwselect_worker using `"`jobdir'"', worker(1) workers(1) outputs(1) starts(1) ///
    jobids(1 3 5 7) cachemb(64)
jwselect_collect using `"`jobdir'"'
assert r(N_jobs)==8
matrix JW030_COLLECT = r(table)
matrix JW030_JOBS = r(jobs)
mata: jw030_equal("JW030_COLLECT","JW030_PART",1e-10)
assert rowsof(JW030_JOBS)==8
assert `"`c(rngstate)'"' == `"`old_rng'"'

// Zero optional moments: no optimizer jobs, normal inference for selected output.
matrix JW_B0 = (0,.5,0,1)
matrix JW_V0 = (0,0,0,0 \ 0,.04,0,.015 \ 0,0,0,0 \ 0,.015,0,.09)
matrix JW_C0 = (3,2 \ 3,3 \ 4,3 \ 4,4)
local zeroinput bmat(JW_B0) vmat(JW_V0) cells(JW_C0) nunits(100)
jwselect, `zeroinput' select(focused) focus(both) inference(reselect) `numerical' inferoutputs(1)
assert e(K)==0
matrix JW030_ZERO = e(table)
local jobdir `"`outdir'/zero_optional"'
jwselect, `zeroinput' select(focused) focus(both) inference(none) `numerical' ///
    inferoutputs(1) prepare(`"`jobdir'"')
jwselect_worker using `"`jobdir'"', worker(1) workers(1)
jwselect_collect using `"`jobdir'"'
assert r(N_jobs)==0
matrix JW030_COLLECT = r(table)
mata: jw030_equal("JW030_COLLECT","JW030_ZERO",1e-10)
assert `"`c(rngstate)'"' == `"`old_rng'"'

// Invalid output selection is rejected and leaves the preceding estimate intact.
matrix JW030_PREV = e(b)
capture noisily jwselect, `input' select(focused) inference(none) inferoutputs(999)
assert _rc != 0
mata: jw030_equal("e(b)","JW030_PREV",1e-10)
assert `"`c(rngstate)'"' == `"`old_rng'"'
display as result "Public API, caching, selection scopes, splitting, and restart: PASS"
