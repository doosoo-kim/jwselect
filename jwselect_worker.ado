*! jwselect 0.5.0 23oct2026
program define jwselect_worker, rclass
    version 17.0
    syntax using/, WORKer(integer) WORKers(integer) ///
        [OUTputs(numlist integer) OBJectives(numlist integer) STARTS(numlist integer) JOBIDS(numlist integer) ///
         CACHEMB(real 64) PROGRESS]
    if `workers'<1 | `worker'<1 | `worker'>`workers' error 198
    if `cachemb'<0 | missing(`cachemb') error 198
    * A partial installation must fail before reserving any job directories.
    quietly findfile _jwselect_lockrelease.ado
    local directory `"`using'"'
    _jwselect_jobguard, directory(`"`directory'"') action(verify)
    _jwselect_load
    tempname jobs
    mata: jw_worker_index_stata()
    local ids
    if `njobs'>0 {
    forvalues j=1/`njobs' {
        local id = el(`jobs',`j',1)
        local o = el(`jobs',`j',2)
        local obj = el(`jobs',`j',3)
        local s = el(`jobs',`j',4)
        local keep = mod(`id'-1,`workers')==`worker'-1
        if "`outputs'"!="" {
            local hit : list o in outputs
            local keep = `keep' & `hit'
        }
        if "`objectives'"!="" {
            local hit : list obj in objectives
            local keep = `keep' & `hit'
        }
        if "`starts'"!="" {
            local hit : list s in starts
            local keep = `keep' & `hit'
        }
        if "`jobids'"!="" {
            local hit : list id in jobids
            local keep = `keep' & `hit'
        }
        if `keep' local ids `ids' `id' 
    }
    }
    * Reserve every assigned search job. Disjoint filtered workers may coexist.
    local locked
    local lock_failures 0
    foreach id of local ids {
        capture mkdir `"`directory'/job_`id'.lock"'
        if _rc {
            foreach old of local locked {
                capture noisily _jwselect_lockrelease, path(`"`directory'/job_`old'.lock"')
                local release_rc = _rc
                if `release_rc' {
                    local ++lock_failures
                    di as error `"Lock cleanup could not finish (r(`release_rc')): `directory'/job_`old'.lock"'
                }
                else if !r(removed) local ++lock_failures
            }
            di as error "Job `id' is already locked. Do not run overlapping assignments."
            di as error "Remove a stale lock only after its original process has stopped."
            exit 603
        }
        local locked `locked' `id'
    }
    local rngstate `c(rngstate)'
    capture noisily mata: jw_jobs_worker_ids(st_local("directory"), ///
        strtoreal(st_local("worker")),strtoreal(st_local("workers")), ///
        strtoreal(tokens(st_local("outputs"))),strtoreal(tokens(st_local("objectives"))), ///
        strtoreal(tokens(st_local("starts"))),strtoreal(tokens(st_local("jobids"))),strtoreal(st_local("cachemb")),st_local("progress")!="")
    local rc = _rc
    quietly set rngstate `rngstate'
    foreach id of local locked {
        capture noisily _jwselect_lockrelease, path(`"`directory'/job_`id'.lock"')
        local release_rc = _rc
        if `release_rc' {
            local ++lock_failures
            di as error `"Lock cleanup could not finish (r(`release_rc')): `directory'/job_`id'.lock"'
        }
        else if !r(removed) local ++lock_failures
    }
    if `rc' exit `rc'
    if `lock_failures' {
        di as error "Search computations completed, but `lock_failures' job lock(s) could not be released."
        di as error "Completed results are retained. Collection requires all job locks to be released."
    }
    return scalar assigned = wordcount("`ids'")
    return scalar worker = `worker'
    return scalar workers = `workers'
    return scalar lock_failures = `lock_failures'
    return local directory `"`directory'"'
end
