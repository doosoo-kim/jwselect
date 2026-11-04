*! jwselect 0.9.0 04nov2026: population-share Sun-Abraham scalar inputs
version 16.1
mata:

void jw_sa_iw_error(string scalar s)
{
    errprintf("jwselect_sa_iw: %s\n",s)
    _error(3498)
}

/* Rows are complete independent panels, ordered by id and calendar time.
   The first joint coordinate is the population IW target; remaining
   coordinates are individual cohort/preperiod contrasts. Sampling V has
   no regression degrees-of-freedom multiplier. */
void jw_sa_iw_inputs(string scalar yvar, string scalar gvar,
    string scalar touse, real scalar nt, real scalar first,
    real rowvector window, real colvector requested, real scalar delta,
    string scalar baseline, string scalar dname, string scalar hname,
    string scalar vname, string scalar hmname, string scalar wname,
    string scalar twname, string scalar basename, string scalar targetmapname,
    string scalar kname, string scalar sharevname, string scalar fixedvname)
{
    real matrix raw,Y,H,IF,weights,targetweights
    real colvector G,groups,eligible,ci,gi,dd,psi,share,base,ns,dk
    real rowvector h
    real scalar n,j,g,t,b,lo,hi,k,pos,pK,d,w,muG,muC

    raw=st_data(.,(yvar,gvar),touse)
    n=rows(raw)/nt
    if (n!=floor(n) | n<=1) jw_sa_iw_error("invalid balanced panel")
    Y=colshape(raw[,1],nt)
    G=raw[(0::(n-1)):*nt:+1,2]
    groups=uniqrows(sort(select(G,G:>0),1))
    ci=selectindex(G:==0)
    if (length(ci)<2) jw_sa_iw_error("at least two never-treated units are required")
    eligible=select(groups,(groups:+window[1]:>=first):&
        (groups:+window[2]:<=first+nt-1))
    if (length(requested)) {
        if (rows(uniqrows(sort(requested,1)))!=rows(requested))
            jw_sa_iw_error("cohorts() must not repeat a cohort")
        for (j=1;j<=length(requested);j++) {
            if (sum(eligible:==requested[j])!=1)
                jw_sa_iw_error("a requested cohort lacks the complete reporting window")
        }
        groups=requested
    }
    else groups=eligible
    if (!length(groups)) jw_sa_iw_error("no cohort supports the complete reporting window")
    base=groups:-delta:-1
    if (baseline=="common") base=J(rows(groups),1,min(groups)-delta-1)
    if (min(base)<first) jw_sa_iw_error("a maintained baseline is outside observed support")
    H=J(0,2,.)
    weights=J(rows(groups),5,.)
    targetweights=J(0,3,.)
    ns=J(rows(groups),1,.)
    dk=J(rows(groups),1,.)
    for (j=1;j<=rows(groups);j++) {
        g=groups[j]
        ns[j]=sum(G:==g)
        if (ns[j]<2) jw_sa_iw_error("each included cohort needs at least two units")
        for (t=first;t<g-delta;t++) if (t!=base[j]) H=H\(g,t)
    }
    k=rows(H)
    pK=sum(ns)/n
    IF=J(n,k+1,0)
    h=J(1,k,0)
    psi=J(n,1,0)
    share=J(n,1,0)
    d=0
    pos=0
    for (j=1;j<=rows(groups);j++) {
        g=groups[j]
        gi=selectindex(G:==g)
        b=base[j]-first+1
        lo=g+window[1]-first+1
        hi=g+window[2]-first+1
        dd=rowsum(Y[,lo..hi])/(hi-lo+1)-Y[,b]
        muG=mean(dd[gi]); muC=mean(dd[ci])
        dk[j]=muG-muC
        w=ns[j]/sum(ns)
        d=d+w*dk[j]
        psi[gi]=psi[gi]+w*n/ns[j]:*(dd[gi]:-muG)
        psi[ci]=psi[ci]-w*n/length(ci):*(dd[ci]:-muC)
        weights[j,.]=(g,ns[j],ns[j]/n,w,dk[j])
        for (t=lo;t<=hi;t++) targetweights=targetweights\(g,t+first-1,w/(hi-lo+1))
        for (t=first;t<g-delta;t++) {
            if (t==base[j]) continue
            pos++
            dd=Y[,t-first+1]-Y[,b]
            muG=mean(dd[gi]); muC=mean(dd[ci])
            h[pos]=muG-muC
            IF[gi,pos+1]=n/ns[j]:*(dd[gi]:-muG)
            IF[ci,pos+1]=-n/length(ci):*(dd[ci]:-muC)
        }
    }
    /* d(weights)/dP contributes this term. Keeping the whole IF retains
       all cross terms; no orthogonality is assumed by the calculation. */
    for (j=1;j<=rows(groups);j++) {
        gi=selectindex(G:==groups[j])
        share[gi]=J(length(gi),1,(dk[j]-d)/pK)
    }
    IF[,1]=psi+share
    st_matrix(dname,d)
    st_matrix(hname,h)
    st_matrix(vname,quadcross(IF,IF)/n^2)
    st_matrix(hmname,H)
    st_matrix(wname,weights)
    st_matrix(twname,targetweights)
    st_matrix(basename,(groups,base))
    st_matrix(targetmapname,(groups[1],groups[1]+window[1]))
    st_numscalar(kname,k)
    st_numscalar(sharevname,quadcross(share,share)/n^2)
    st_numscalar(fixedvname,quadcross(psi,psi)/n^2)
}

string scalar jw_sa_iw_version() return("0.9.0")
end
