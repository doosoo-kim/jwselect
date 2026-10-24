/* Stata transport only. Numerical modules never depend on st_* state. */
real matrix jw_optional_matrix(string scalar localname)
{
    string scalar name
    name=st_local(localname)
    if (name=="") return(J(0,0,.))
    return(st_matrix(name))
}

void jw_put(string scalar localname, real matrix value)
{
    if (rows(value)>0 & cols(value)>0) st_matrix(st_local(localname),value)
}

void jw_names(string scalar localname, string colvector names, real scalar rownames)
{
    string matrix stripe
    stripe=J(rows(names),1,"") , names
    if (rownames) st_matrixrowstripe(st_local(localname),stripe)
    else st_matrixcolstripe(st_local(localname),stripe)
}

void jw_stata()
{
    struct jw_model scalar M
    struct jw_family scalar C
    struct jw_fit scalar F
    struct jw_reselect_result scalar RI
    real matrix table, Vout, Vcells, joint, pw
    real rowvector outputs
    real colvector estimate
    real scalar i, j, G, m
    string colvector names, cnames, hnames, znames, familycols, diagcols
    string rowvector scorecols
    string scalar selector, focus, inference

    selector=st_local("select")
    focus=st_local("focus")
    inference=st_local("inference")
    pw=J(1,0,.)
    if (st_local("postwindow")!="") pw=strtoreal(tokens(st_local("postwindow")))
    M=jw_design(st_matrix(st_local("rawb")),st_matrix(st_local("rawV")),
        st_matrix(st_local("rawcells")),strtoreal(st_local("nunits")),
        strtoreal(st_local("gap")),strtoreal(st_local("anticipation")),
        st_local("baseline"),pw,st_local("aggregate"),jw_optional_matrix("aggweights"))
    C=jw_candidates(M,selector,st_local("candidates"),jw_optional_matrix("prematrix"),
        jw_optional_matrix("candmatrix"),strtoreal(st_local("maxcand")))
    joint=jw_optional_matrix("jointmatrix")
    F=jw_run_scope(M,C,joint,selector,focus,strtoreal(st_local("maxjoint")),st_local("wbicscope"))
    outputs=J(1,0,.)
    if (st_local("inferoutputs")!="") outputs=strtoreal(tokens(st_local("inferoutputs")))
    if (cols(outputs)>0) {
        if (hasmissing(outputs) | min(outputs)<1 | max(outputs)>rows(F.A) |
            any(outputs:!=floor(outputs))) jw_core_error("invalid inferoutputs() row index")
    }
    if (st_local("prepare")!="") {
        jw_jobs_prepare(M,F,st_local("prepare"),strtoreal(st_local("level")),
            strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),
            strtoreal(st_local("reps")),strtoreal(st_local("optstarts")),
            strtoreal(st_local("optmaxeval")),strtoreal(st_local("opttol")),
            strtoreal(st_local("optreps")),outputs,strtoreal(st_local("cachemb")),
            st_local("progress")!="")
    }
    estimate=F.A*M.z
    Vout=F.A*M.V*F.A'
    Vout=(Vout+Vout')/2
    for (i=1; i<=rows(F.A); i++) Vout[i,i]=jw_fixed_variance(M.V,F.A[i,.])
    Vcells=F.cellA*M.V*F.cellA'
    Vcells=(Vcells+Vcells')/2
    for (i=1; i<=M.p; i++) Vcells[i,i]=jw_fixed_variance(M.V,F.cellA[i,.])
    if (inference=="none") {
        table=J(rows(F.A),5,.)
        table[.,1]=estimate
        for (i=1; i<=rows(F.A); i++) table[i,2]=sqrt(jw_fixed_variance(M.V,F.A[i,.]))
    }
    else if (selector=="focused" & (inference=="auto" | inference=="reselect")) {
        RI=jw_infer_reselect_opts(M,F,strtoreal(st_local("level")),
            strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),
            strtoreal(st_local("reps")),strtoreal(st_local("optstarts")),
            strtoreal(st_local("optmaxeval")),strtoreal(st_local("opttol")),
            strtoreal(st_local("optreps")),outputs,strtoreal(st_local("cachemb")),
            st_local("progress")!="")
        table=RI.table
        jw_put("outdiagnostics",RI.diagnostics)
        jw_put("outextrema",RI.extrema)
    }
    else table=jw_infer(M,F,selector,strtoreal(st_local("level")),
        strtoreal(st_local("gamma")),strtoreal(st_local("mcerr")),strtoreal(st_local("reps")))
    G=rows(M.groups)
    m=rows(M.z)
    names=J(rows(F.A),1,"")
    for (i=1; i<=rows(F.A); i++) {
        if (F.outmap[i,1]==1 | F.outmap[i,1]==5) names[i]="cohort"+strtrim(strofreal(i))
        else if (F.outmap[i,1]==2) names[i]="agg_cohort"
        else if (F.outmap[i,1]==3) names[i]="agg_focused"
        else names[i]="aggregate"
    }
    cnames=J(M.p,1,"")
    for (i=1; i<=M.p; i++) cnames[i]="cell"+strtrim(strofreal(i))
    hnames=J(M.k,1,"")
    for (i=1; i<=M.k; i++) hnames[i]="h"+strtrim(strofreal(i))
    if (selector=="focused" & (inference=="auto" | inference=="reselect")) {
        diagcols=("output"\"objective"\"nactive"\"search_reps"\"validation_reps"\
            "evaluations"\"starts_used"\"converged_starts"\"budget_exhausted"\
            "search_value"\"validated_value"\"validation_gap"\"radius_fraction")
        jw_names("outdiagnostics",diagcols,0)
        if (rows(RI.extrema)>0) jw_names("outextrema",("output"\"objective"\hnames),0)
    }
    znames=cnames\hnames
    jw_put("outb",estimate')
    jw_put("outV",Vout)
    jw_put("outtable",table)
    jw_put("outmap",F.outmap)
    jw_put("outselection",F.selected)
    jw_put("outbases",(M.groups,M.bases))
    jw_put("outpost",M.postmap)
    jw_put("outh",M.hmap)
    jw_put("outbcells",(F.cellA*M.z)')
    jw_put("outVcells",Vcells)
    jw_put("outcellsel",F.cell_selected)
    jw_put("outweights",M.w)
    jw_put("outscores",F.scores)
    jw_put("outfamily",(C.owners,C.ids,C.masks))
    if ((selector=="focused" & focus!="cohort") | (selector=="wbic" & st_local("wbicscope")=="joint")) {
        joint=jw_joint_masks(M,C,joint,strtoreal(st_local("maxjoint")))
        jw_put("outjoint",joint)
        if (M.k>0) jw_names("outjoint",hnames,0)
    }
    jw_put("outz",M.z')
    jw_put("outzV",M.V)
    jw_put("outA",F.A)
    jw_put("outT",F.targets)
    jw_names("outb",names,0)
    jw_names("outV",names,0)
    jw_names("outV",names,1)
    jw_names("outtable",names,1)
    jw_names("outtable",("estimate"\"se_fixed"\"ll"\"ul"\"pvalue"),0)
    jw_names("outmap",names,1)
    jw_names("outmap",("kind"\"cohort"),0)
    jw_names("outbases",("cohort"\"baseline"),0)
    jw_names("outpost",("cohort"\"time"),0)
    jw_names("outpost",cnames,1)
    if (M.k>0) {
        jw_names("outh",("cohort"\"time"),0)
        jw_names("outh",hnames,1)
        jw_names("outselection",names,1)
        jw_names("outselection",hnames,0)
        jw_names("outcellsel",cnames,1)
        jw_names("outcellsel",hnames,0)
    }
    jw_names("outbcells",cnames,0)
    jw_names("outVcells",cnames,0)
    jw_names("outVcells",cnames,1)
    jw_names("outweights",cnames,0)
    scorecols=("scope","cohort","candidate","nrestr","Wald","WBIC","estimate","se_fixed","adjustment","gain","FMSC")
    jw_names("outscores",scorecols',0)
    familycols=("cohort"\"candidate"\hnames)
    jw_names("outfamily",familycols,0)
    jw_names("outz",znames,0)
    jw_names("outzV",znames,0)
    jw_names("outzV",znames,1)
    jw_names("outA",names,1)
    jw_names("outA",znames,0)
    jw_names("outT",names,1)
    jw_names("outT",znames,0)
    st_local("_G",strofreal(G))
    st_local("_K",strofreal(M.k))
    st_local("_P",strofreal(M.p))
}

/* r-class transport for completed divided searches; no e() mutation. */
void jw_worker_index_stata()
{
    real matrix index
    index=jw_jobs_index(st_local("directory"))
    st_local("njobs",strtrim(strofreal(rows(index))))
    if (rows(index)>0) st_matrix(st_local("jobs"),index)
}

void jw_collect_stata()
{
    struct jw_reselect_result scalar R
    string scalar path
    real matrix index
    path=st_local("directory")
    R=jw_jobs_collect(path,st_local("progress")!="")
    st_matrix(st_local("table"),R.table)
    st_matrix(st_local("diagnostics"),R.diagnostics)
    st_matrix(st_local("extrema"),R.extrema)
    st_matrix(st_local("output_map"),jw_jobs_output_map(path))
    index=jw_jobs_index(path)
    st_local("njobs",strtrim(strofreal(rows(index))))
    if (rows(index)==0) index=J(1,4,.)
    st_matrix(st_local("jobs"),index)
    st_matrix(st_local("settings"),jw_jobs_settings(path))
}

/* Joint testing uses the saved canonical contrasts without refitting,
   changing e(), or drawing random numbers. Calendar-date subsets map only
   to eligible optional contrasts; maintained baselines are never tested. */
void jw_pretrend_stata()
{
    real colvector z, ii
    real matrix V, hm, requested, selections
    real rowvector mask, result
    real scalar k, oi, i, h
    string scalar testset

    z = st_matrix("e(z)")'
    V = st_matrix("e(V_z)")
    k = st_numscalar("e(K)")
    mask = J(1,k,1)
    hm = J(0,2,.)
    if (k>0) hm = st_matrix("e(hmap)")
    if (rows(hm)!=k | cols(hm)!=2) {
        jw_core_error("saved optional-preperiod map is incompatible")
    }
    testset = st_local("testset")
    if (testset=="selected") {
        oi = strtoreal(st_local("output"))
        if (k>0) {
            selections = st_matrix("e(selection)")
            if (oi<1 | oi>rows(selections) | cols(selections)!=k) {
                jw_core_error("selected output has no compatible restriction mask")
            }
            mask = selections[oi,.]
        }
    }
    else if (testset=="custom") {
        requested = st_matrix(st_local("prematrix"))
        if (cols(requested)!=2 | hasmissing(requested)) {
            jw_core_error("prematrix() must contain finite cohort and calendar-time pairs")
        }
        mask = J(1,k,0)
        for (i=1; i<=rows(requested); i++) {
            h = jw_findcell(hm,requested[i,1],requested[i,2],st_numscalar("e(gap)"))
            if (h==0) {
                jw_core_error("prematrix() contains an ineligible date or maintained baseline")
            }
            if (mask[h]==1) jw_core_error("prematrix() contains duplicate cohort-time pairs")
            mask[h] = 1
        }
    }
    result = jw_pretrend_test(z,V,mask,k)
    st_matrix(st_local("testtable"),result)
    if (k>0) st_matrix(st_local("testmask"),mask)
    if (sum(mask)>0) {
        ii = vec(selectindex(mask'))
        st_matrix(st_local("testperiods"),hm[ii,.])
    }
}

/* Completion marker must be LAST: a partially compiled bridge is not ready. */
string scalar jw_version()
{
    return("0.4.0-joint1")
}
