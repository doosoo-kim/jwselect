*! jwselect 0.7.0 23oct2026
program define _jwselect_jobguard, rclass
    version 17.0
    syntax, DIRectory(string) ACTion(string)
    if !inlist("`action'","finalize","verify","unlocked") error 198
    if "`action'"=="unlocked" {
        local locks : dir `"`directory'"' dirs "job_*.lock"
        if `"`locks'"'!="" {
            di as error "A search job is active or has a stale lock. Stop its process before removing a stale lock."
            exit 603
        }
        exit
    }
    if "`action'"=="finalize" {
        confirm file `"`directory'/jw_state.bin"'
        mkdir `"`directory'/code"'
        quietly findfile _jwselect_mata.do
        local loader `"`r(fn)'"'
        copy `"`loader'"' `"`directory'/code/_jwselect_mata.do"'
        mata: st_local("files",invtokens(jw_jobs_files(st_local("directory"))'))
        local files `files' code/_jwselect_mata.do
        tempname mh
        * A type prefix applies to one variable unless the names are grouped.
        * File lengths above 2^24 must not be rounded to default float storage.
        postfile `mh' str244 relative_path double(checksum file_bytes) using `"`directory'/manifest.dta"'
        foreach f of local files {
            quietly checksum `"`directory'/`f'"'
            post `mh' ("`f'") (r(checksum)) (r(filelen))
        }
        postclose `mh'
        quietly checksum `"`directory'/manifest.dta"'
        local manifest_crc = r(checksum)
        local manifest_bytes = r(filelen)
        preserve
        capture noisily {
            clear
            set obs 1
            generate str24 build = "0.7.0-index1"
            generate double manifest_crc = `manifest_crc'
            generate double manifest_bytes = `manifest_bytes'
            save `"`directory'/READY.dta"'
        }
        local rc = _rc
        restore
        if `rc' exit `rc'
        exit
    }
    confirm file `"`directory'/READY.dta"'
    confirm file `"`directory'/manifest.dta"'
    preserve
    capture noisily {
        quietly use `"`directory'/READY.dta"', clear
        assert _N==1 & build=="0.7.0-index1"
        local crc = manifest_crc[1]
        local bytes = manifest_bytes[1]
        quietly checksum `"`directory'/manifest.dta"'
        assert r(checksum)==`crc' & r(filelen)==`bytes'
        quietly use `"`directory'/manifest.dta"', clear
        isid relative_path
        assert _N>=2
        quietly count if relative_path=="jw_state.bin"
        assert r(N)==1
        quietly count if relative_path=="code/_jwselect_mata.do"
        assert r(N)==1
        forvalues i=1/`=_N' {
            local f = relative_path[`i']
            local crc = checksum[`i']
            local bytes = file_bytes[`i']
            quietly checksum `"`directory'/`f'"'
            assert r(checksum)==`crc' & r(filelen)==`bytes'
            if "`f'"=="code/_jwselect_mata.do" {
                local code_crc = `crc'
                local code_bytes = `bytes'
            }
        }
        quietly findfile _jwselect_mata.do
        local loader `"`r(fn)'"'
        quietly checksum `"`loader'"'
        assert r(checksum)==`code_crc' & r(filelen)==`code_bytes'
    }
    local rc = _rc
    restore
    if `rc' {
        di as error "Prepared inputs or installed jwselect code do not match the frozen batch."
        di as error "Use its original package in a fresh session; do not edit prepared files."
        exit `rc'
    }
end
