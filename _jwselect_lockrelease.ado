*! jwselect 0.3.0 06oct2026
program define _jwselect_lockrelease, rclass
    version 17.0
    * Internal API: call only for a directory acquired by this process.
    * rmdir is nonrecursive. Never erase contents or an adjacent lock.
    syntax, PATH(string) [MAXAttempts(integer 5) PAUSE(integer 500)]
    if `"`path'"'=="" error 198
    if missing(`maxattempts') | `maxattempts'<1 | `maxattempts'>10 error 198
    if missing(`pause') | `pause'<0 | `pause'>1000 error 198
    if (`maxattempts'-1)*`pause'>5000 error 198

    local attempts 0
    local last_rc 0
    local reappeared 0
    * Builtin Mata only: no jwselect loader, data changes, or RNG draws.
    mata: st_local("exists", strofreal(direxists(st_local("path"))))
    while `exists' & `attempts'<`maxattempts' {
        local ++attempts
        capture rmdir `"`path'"'
        local last_rc = _rc
        mata: st_local("exists", strofreal(direxists(st_local("path"))))
        if !`exists' continue, break
        if `last_rc'==0 {
            * Another process may have acquired this name after our removal.
            * Do not attempt to remove a directory whose ownership is unknown.
            local reappeared 1
            continue, break
        }
        if `attempts'<`maxattempts' & `pause'>0 sleep `pause'
    }
    if `exists' {
        if `reappeared' {
            display as error "Lock cleanup warning: directory is present after successful removal; no further removal attempted."
        }
        else {
            display as error "Lock cleanup warning: directory remains after `attempts' attempts (last r(`last_rc'))."
        }
        display as error `"Lock path: `path'"'
        display as error "Confirm all processes using this lock have stopped before inspecting or removing it."
    }
    return scalar removed = !`exists'
    return scalar attempts = `attempts'
    return scalar last_rc = `last_rc'
    return scalar reappeared = `reappeared'
end
