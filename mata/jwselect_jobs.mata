/* jwselect 0.3.0: portable, resumable objective-by-start jobs.

   Only numeric/string matrices are serialized with fputmatrix, not opaque
   compiled code or raw observations. Prepare and serial inference consume
   random numbers in the same order. Workers consume no RNG. The Stata
   interface protects immutable files with checksums and acquires job locks.
   A .done marker is created only after its result file has been closed.
   An interrupted unmarked result can be replaced; a marked corrupt result
   is an error. No projection cache is written to disk.
*/

struct jw_job_state {
    real rowvector token
    real rowvector settings
    real rowvector numerical
    real matrix index
    real matrix info
    real matrix outmap
    struct jw_reselect_result scalar result
}

string scalar jw_jobs_num(real scalar x)
{
    return(strtrim(strofreal(x,"%18.0f")))
}

string scalar jw_jobs_output_file(string scalar path, real scalar o)
{
    return(path+"/output_"+jw_jobs_num(o)+".bin")
}

string scalar jw_jobs_result_file(string scalar path, real scalar id)
{
    return(path+"/job_"+jw_jobs_num(id)+".bin")
}

void jw_jobs_write_context(real scalar fh, struct jw_rs_context scalar C)
{
    fputmatrix(fh,C.reportB)
    fputmatrix(fh,C.selectB)
    fputmatrix(fh,C.constant)
    fputmatrix(fh,C.component)
    fputmatrix(fh,C.groups)
    fputmatrix(fh,C.EH)
    fputmatrix(fh,C.baseline)
}

struct jw_rs_context scalar jw_jobs_read_context(real scalar fh)
{
    struct jw_rs_context scalar C
    C.reportB = fgetmatrix(fh)
    C.selectB = fgetmatrix(fh)
    C.constant = fgetmatrix(fh)
    C.component = fgetmatrix(fh)
    C.groups = fgetmatrix(fh)
    C.EH = fgetmatrix(fh)
    C.baseline = fgetmatrix(fh)
    C.projection = J(0,0,.)
    if (rows(C.reportB)!=rows(C.selectB) |
        cols(C.reportB)!=cols(C.selectB) |
        rows(C.constant)!=rows(C.reportB) | cols(C.constant)!=1 |
        rows(C.component)!=rows(C.reportB) | cols(C.component)!=1 |
        cols(C.EH)!=cols(C.reportB) |
        rows(C.baseline)!=rows(C.EH) | cols(C.baseline)!=1 |
        hasmissing(C.reportB) | hasmissing(C.selectB) |
        hasmissing(C.constant) | hasmissing(C.EH) |
        hasmissing(C.baseline)) {
        jw_core_error("invalid or incompatible frozen search context")
    }
    return(C)
}

void jw_jobs_write_output(string scalar filename, struct jw_rs_output scalar P)
{
    real scalar fh
    fh = fopen(filename,"w")
    fputmatrix(fh,"jwselect-output-0.3.0")
    fputmatrix(fh,(P.output,P.estimate))
    fputmatrix(fh,P.active)
    fputmatrix(fh,P.hcols)
    fputmatrix(fh,P.fullmu)
    fputmatrix(fh,P.h)
    fputmatrix(fh,P.B)
    fputmatrix(fh,P.starts)
    jw_jobs_write_context(fh,P.C)
    jw_jobs_write_context(fh,P.CV)
    fputmatrix(fh,"complete")
    fclose(fh)
}

struct jw_rs_output scalar jw_jobs_read_output(string scalar path,
                                               real scalar o)
{
    struct jw_rs_output scalar P
    real scalar fh, d
    real rowvector meta
    string scalar magic
    fh = fopen(jw_jobs_output_file(path,o),"r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-output-0.3.0") jw_core_error("incompatible frozen output format")
    meta = fgetmatrix(fh)
    if (cols(meta)!=2 | rows(meta)!=1) jw_core_error("invalid frozen output header")
    P.output = meta[1]
    P.estimate = meta[2]
    P.active = fgetmatrix(fh)
    P.hcols = fgetmatrix(fh)
    P.fullmu = fgetmatrix(fh)
    P.h = fgetmatrix(fh)
    P.B = fgetmatrix(fh)
    P.starts = fgetmatrix(fh)
    P.C = jw_jobs_read_context(fh)
    P.CV = jw_jobs_read_context(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    d = cols(P.active)
    if (magic!="complete" | P.output!=o | d<1 |
        rows(P.h)!=d | cols(P.h)!=1 | rows(P.B)!=d | cols(P.B)!=d |
        rows(P.starts)!=d | cols(P.starts)<1 | missing(P.estimate) |
        hasmissing(P.h) | hasmissing(P.B) | hasmissing(P.starts) |
        cols(P.C.EH)!=d | cols(P.CV.EH)!=d) {
        jw_core_error("invalid or incomplete frozen output")
    }
    return(P)
}

void jw_jobs_write_state(string scalar path, struct jw_job_state scalar JS)
{
    real scalar fh
    fh = fopen(path+"/jw_state.bin","w")
    fputmatrix(fh,"jwselect-jobs-0.3.0")
    fputmatrix(fh,JS.token)
    fputmatrix(fh,JS.settings)
    fputmatrix(fh,JS.numerical)
    fputmatrix(fh,JS.index)
    fputmatrix(fh,JS.info)
    fputmatrix(fh,JS.outmap)
    fputmatrix(fh,JS.result.table)
    fputmatrix(fh,JS.result.diagnostics)
    fputmatrix(fh,JS.result.extrema)
    fputmatrix(fh,"complete")
    fclose(fh)
}

struct jw_job_state scalar jw_jobs_read_state(string scalar path)
{
    struct jw_job_state scalar JS
    real scalar fh
    string scalar magic
    fh = fopen(path+"/jw_state.bin","r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-jobs-0.3.0") jw_core_error("incompatible jobs state format")
    JS.token = fgetmatrix(fh)
    JS.settings = fgetmatrix(fh)
    JS.numerical = fgetmatrix(fh)
    JS.index = fgetmatrix(fh)
    JS.info = fgetmatrix(fh)
    JS.outmap = fgetmatrix(fh)
    JS.result.table = fgetmatrix(fh)
    JS.result.diagnostics = fgetmatrix(fh)
    JS.result.extrema = fgetmatrix(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    if (magic!="complete" | cols(JS.settings)!=9 | rows(JS.settings)!=1 |
        cols(JS.numerical)!=4 | rows(JS.numerical)!=1 |
        cols(JS.index)!=4 | cols(JS.info)!=3 |
        cols(JS.result.table)!=5 | cols(JS.result.diagnostics)!=13 |
        rows(JS.result.diagnostics)!=4*rows(JS.result.table) |
        rows(JS.result.extrema)!=4*rows(JS.result.table) |
        hasmissing(JS.token) | hasmissing(JS.settings) |
        hasmissing(JS.numerical) | hasmissing(JS.index)) {
        jw_core_error("invalid or incomplete jobs state")
    }
    return(JS)
}

real matrix jw_jobs_index(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.index)
}

real matrix jw_jobs_output_map(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.outmap)
}

real matrix jw_jobs_output_info(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.info)
}

real rowvector jw_jobs_settings(string scalar path)
{
    struct jw_job_state scalar JS
    JS = jw_jobs_read_state(path)
    return(JS.settings)
}

string colvector jw_jobs_files(string scalar path)
{
    struct jw_job_state scalar JS
    string colvector files
    real scalar o, n, i
    JS = jw_jobs_read_state(path)
    n = sum(JS.info[,2]:>0)
    files = J(n+1,1,"")
    files[1] = "jw_state.bin"
    i = 1
    for (o=1; o<=rows(JS.info); o++) {
        if (JS.info[o,2]>0) {
            i++
            files[i] = "output_"+jw_jobs_num(JS.info[o,1])+".bin"
        }
    }
    return(files)
}

/* Preparation contains no optimization, and never evaluates an omitted
   output. Empty outputs means all. The caller creates a new directory. */
void jw_jobs_prepare(struct jw_model scalar M, struct jw_fit scalar F,
                    string scalar path, real scalar level,
                    real scalar gamma, real scalar mcerr, real scalar reps,
                    real scalar optstarts, real scalar optmaxeval,
                    real scalar opttol, real scalar optreps,
                    real rowvector outputs, real scalar cachemb,
                    real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    real matrix S, E, EV, index
    real scalar o, j, s, d, initialized, id, nmax, fingerprint
    if (!direxists(path)) jw_core_error("jobs directory does not exist")
    if (fileexists(path+"/jw_state.bin")) {
        jw_core_error("prepared state already exists; use a new directory")
    }
    JS.numerical = jw_rs_validate(M,F,level,gamma,mcerr,reps,optstarts,
                        optmaxeval,opttol,optreps,outputs,cachemb)
    JS.settings = (level,gamma,mcerr,reps,optstarts,optmaxeval,opttol,optreps,cachemb)
    JS.result = jw_rs_result_init(M,F,level)
    JS.outmap = F.outmap
    JS.info = J(rows(F.A),3,0)
    nmax = rows(F.A)*4*optstarts
    index = J(nmax,4,.)
    id = 0
    initialized = 0
    fingerprint = 0
    S = J(0,0,.)
    E = J(0,0,.)
    for (o=1; o<=rows(F.A); o++) {
        JS.info[o,1] = o
        if (cols(outputs)>0) {
            if (!any(outputs:==o)) continue
        }
        P = jw_rs_output_plan(M,F,o,gamma)
        d = cols(P.active)
        if (d==0) {
            jw_rs_zero_output(JS.result,P,jw_fixed_variance(M.V,F.A[o,.]),level)
            continue
        }
        if (JS.numerical[3]<1 | JS.numerical[4]>reps) {
            jw_core_error("reps() too small for finite Monte Carlo quantile bounds; increase reps()")
        }
        if (!initialized) {
            S = jw_infer_sqrt(M.V)
            E = rnormal(optreps,rows(M.z),0,1)*S'
            initialized = 1
            fingerprint = sum(E)
        }
        P.C = jw_rs_prepare_opts(M,F,o,E,P.hcols,0)
        P.starts = jw_rs_starts(P.h,P.B,optstarts)
        EV = rnormal(reps,rows(M.z),0,1)*S'
        P.CV = jw_rs_prepare_opts(M,F,o,EV,P.hcols,0)
        fingerprint = fingerprint+o*sum(EV)+sum(P.starts)
        jw_jobs_write_output(jw_jobs_output_file(path,o),P)
        JS.info[o,.] = (o,d,4*optstarts)
        for (j=1; j<=4; j++) {
            for (s=1; s<=optstarts; s++) {
                id++
                index[id,.] = (id,o,j,s)
            }
        }
        if (progress) {
            printf("Prepared output %g: %g objectives x %g starts, dimension %g\n",
                    o,4,optstarts,d)
            displayflush()
        }
    }
    JS.index = J(0,4,.)
    if (id>0) JS.index = index[|1,1\id,4|]
    /* No random draw is consumed for provenance. Together with the front-
       end's immutable-file checksums this catches copied/mismatched jobs. */
    JS.token = (300,JS.settings[1..8],rows(M.z),M.k,rows(F.A),
                sum(M.z),sum(M.V),sum(F.A),fingerprint)
    jw_jobs_write_state(path,JS)
}

void jw_jobs_save_result(string scalar path, real rowvector token,
                         real rowvector job,
                         struct jw_rs_search_result scalar R)
{
    real scalar fh
    string scalar filename
    filename = jw_jobs_result_file(path,job[1])
    if (fileexists(filename+".done")) {
        jw_core_error("refusing to overwrite a completed search job")
    }
    /* fopen(...,"w") requires a new file. Only an unmarked partial file
       may be removed here, after the front-end acquires this job's lock. */
    unlink(filename)
    fh = fopen(filename,"w")
    fputmatrix(fh,"jwselect-search-0.3.0")
    fputmatrix(fh,token)
    fputmatrix(fh,job)
    fputmatrix(fh,R.x)
    fputmatrix(fh,(R.value,R.evaluations,R.starts,R.converged,R.exhausted))
    fputmatrix(fh,"complete")
    fclose(fh)
    /* The marker is deliberately written last. */
    fh = fopen(filename+".done","w")
    fputmatrix(fh,token)
    fputmatrix(fh,job)
    fclose(fh)
}

struct jw_rs_search_result scalar jw_jobs_load_result(string scalar path,
                    real rowvector token, real rowvector job,
                    real scalar d, real scalar budget)
{
    struct jw_rs_search_result scalar R
    real scalar fh
    real rowvector savedtoken, savedjob, meta
    string scalar filename, magic
    filename = jw_jobs_result_file(path,job[1])
    if (!fileexists(filename+".done")) jw_core_error("required search job is unfinished")
    fh = fopen(filename+".done","r")
    savedtoken = fgetmatrix(fh)
    savedjob = fgetmatrix(fh)
    fclose(fh)
    if (rows(savedtoken)!=1 | cols(savedtoken)!=cols(token) |
        rows(savedjob)!=1 | cols(savedjob)!=4) jw_core_error("invalid job completion marker")
    if (any(savedtoken:!=token) | any(savedjob:!=job)) {
        jw_core_error("job completion marker belongs to a different preparation or job")
    }
    fh = fopen(filename,"r")
    magic = fgetmatrix(fh)
    if (magic!="jwselect-search-0.3.0") jw_core_error("incompatible search result format")
    savedtoken = fgetmatrix(fh)
    savedjob = fgetmatrix(fh)
    R.x = fgetmatrix(fh)
    meta = fgetmatrix(fh)
    magic = fgetmatrix(fh)
    fclose(fh)
    if (magic!="complete" | rows(savedtoken)!=1 | cols(savedtoken)!=cols(token) |
        rows(savedjob)!=1 | cols(savedjob)!=4 | rows(R.x)!=d | cols(R.x)!=1 |
        rows(meta)!=1 | cols(meta)!=5 | hasmissing(R.x) | hasmissing(meta)) {
        jw_core_error("invalid or incomplete completed search result")
    }
    if (any(savedtoken:!=token) | any(savedjob:!=job) |
        sqrt(sum(R.x:^2))>1+1e-10 | meta[2]<1 | meta[2]>budget |
        meta[2]!=floor(meta[2]) | meta[3]!=1 |
        (meta[4]!=0 & meta[4]!=1) | (meta[5]!=0 & meta[5]!=1) |
        meta[4]+meta[5]!=1 | (meta[5]==1 & meta[2]!=budget)) {
        jw_core_error("completed search result is incompatible with its job")
    }
    R.value = meta[1]
    R.evaluations = meta[2]
    R.starts = meta[3]
    R.converged = meta[4]
    R.exhausted = meta[5]
    return(R)
}

void jw_jobs_check_filter(real rowvector filter, real scalar upper,
                           string scalar label)
{
    if (cols(filter)>0) {
        if (hasmissing(filter) | min(filter)<1 | max(filter)>upper |
            any(filter:!=floor(filter))) {
            jw_core_error("invalid "+label+" filter")
        }
    }
}

/* Filters are rowvectors; empty means all. Modulo assignment always uses
   the original global job id, independently of filters. Root owns locks. */
void jw_jobs_worker_ids(string scalar path, real scalar worker,
                    real scalar nworkers, real rowvector outputfilter,
                    real rowvector objectivefilter, real rowvector startfilter,
                    real rowvector jobfilter, real scalar cachemb,
                    real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    struct jw_rs_search_result scalar R
    real scalar i, o, objective, start, lastoutput, budget
    real rowvector job
    JS = jw_jobs_read_state(path)
    if (worker<1 | worker>nworkers | nworkers<1 |
        worker!=floor(worker) | nworkers!=floor(nworkers) |
        missing(worker) | missing(nworkers) | missing(cachemb) | cachemb<0) {
        jw_core_error("invalid worker index, worker count, or cache budget")
    }
    jw_jobs_check_filter(outputfilter,rows(JS.info),"output")
    jw_jobs_check_filter(objectivefilter,4,"objective")
    jw_jobs_check_filter(startfilter,JS.settings[5],"starting-point")
    jw_jobs_check_filter(jobfilter,rows(JS.index),"job")
    lastoutput = 0
    for (i=1; i<=rows(JS.index); i++) {
        job = JS.index[i,.]
        if (mod(job[1]-1,nworkers)!=worker-1) continue
        if (cols(jobfilter)>0) {
            if (!any(jobfilter:==job[1])) continue
        }
        o = job[2]
        objective = job[3]
        start = job[4]
        if (cols(outputfilter)>0) {
            if (!any(outputfilter:==o)) continue
        }
        if (cols(objectivefilter)>0) {
            if (!any(objectivefilter:==objective)) continue
        }
        if (cols(startfilter)>0) {
            if (!any(startfilter:==start)) continue
        }
        budget = jw_rs_start_budget(JS.settings[6],JS.settings[5],start)
        if (fileexists(jw_jobs_result_file(path,job[1])+".done")) {
            R = jw_jobs_load_result(path,JS.token,job,JS.info[o,2],budget)
            if (progress) printf("Search job %g already complete; retained.\n",job[1])
            continue
        }
        if (o!=lastoutput) {
            P = jw_jobs_read_output(path,o)
            jw_rs_cache(P.C,cachemb)
            lastoutput = o
        }
        if (progress) {
            printf("Search job %g/%g: output %g, objective %g, start %g, budget %g\n",
                    job[1],rows(JS.index),o,objective,start,budget)
            displayflush()
        }
        R = jw_rs_search_start(P.C,P.h,P.B,P.starts[,start],JS.numerical[1],
                P.estimate,objective,budget,JS.settings[7],progress)
        jw_jobs_save_result(path,JS.token,job,R)
        if (progress) {
            printf("Search job %g complete: %g evaluations, converged=%g, budget_exhausted=%g\n",
                    job[1],R.evaluations,R.converged,R.exhausted)
            displayflush()
        }
    }
}

/* Compatibility entry point: no explicit job-id filter. */
void jw_jobs_worker(string scalar path, real scalar worker,
                    real scalar nworkers, real rowvector outputfilter,
                    real rowvector objectivefilter, real rowvector startfilter,
                    real scalar cachemb, real scalar progress)
{
    jw_jobs_worker_ids(path,worker,nworkers,outputfilter,objectivefilter,
                       startfilter,J(1,0,.),cachemb,progress)
}

struct jw_reselect_result scalar jw_jobs_collect(string scalar path,
                                                 real scalar progress)
{
    struct jw_job_state scalar JS
    struct jw_rs_output scalar P
    struct jw_rs_search_result scalar R, S
    real matrix located
    real rowvector searchvalues, job
    real scalar o, objective, start, id, budget, d
    JS = jw_jobs_read_state(path)
    id = 0
    /* Inspect every required completion before doing costly validation. */
    for (o=1; o<=rows(JS.index); o++) {
        if (!fileexists(jw_jobs_result_file(path,JS.index[o,1])+".done")) {
            jw_core_error("cannot collect: one or more required search jobs are unfinished")
        }
    }
    for (o=1; o<=rows(JS.info); o++) {
        d = JS.info[o,2]
        if (d==0) continue
        P = jw_jobs_read_output(path,o)
        located = J(d,4,.)
        searchvalues = J(1,4,.)
        for (objective=1; objective<=4; objective++) {
            R.value = .
            R.x = J(d,1,.)
            R.evaluations = 0
            R.starts = 0
            R.converged = 0
            R.exhausted = 0
            for (start=1; start<=JS.settings[5]; start++) {
                id++
                job = JS.index[id,.]
                if (any(job:!=(id,o,objective,start))) jw_core_error("invalid search job order")
                budget = jw_rs_start_budget(JS.settings[6],JS.settings[5],start)
                S = jw_jobs_load_result(path,JS.token,job,d,budget)
                jw_rs_merge(R,S,objective)
            }
            located[,objective] = P.h+P.B*R.x
            searchvalues[objective] = R.value
            jw_rs_record_search(JS.result,P,R,objective,JS.settings[8],JS.settings[4])
        }
        /* CV is deliberately uncached; only four validation points. */
        jw_rs_finish_output(JS.result,P,located,searchvalues,JS.numerical,
                            JS.settings[2],JS.settings[3],progress)
    }
    return(JS.result)
}
