*! jwselect 0.9.0 04nov2026: transported cell AIPW scores and full sandwich
version 16.1
mata:
mata set matastrict on

real colvector jwrcs_logistic(real colvector a)
{
    real colvector p, ix
    p=J(rows(a),1,.)
    ix=selectindex(a:>=0)
    if (rows(ix)) p[ix]=1:/(1:+exp(-a[ix]))
    ix=selectindex(a:<0)
    if (rows(ix)) p[ix]=exp(a[ix]):/(1:+exp(a[ix]))
    return(p)
}

real scalar jwrcs_ll(real matrix X, real colvector a, real colvector b)
{
    real colvector z
    z=X*b
    return(mean(a:*z-rowmax((z,J(rows(z),1,0)))-ln(1:+exp(-abs(z)))))
}

real matrix jwrcs_inverse(real matrix Q, string scalar what)
{
    real rowvector ev
    ev=symeigenvalues((Q+Q')/2)
    if (hasmissing(ev) | min(ev)<=1e-12*max((1,max(ev)))) {
        errprintf("RCS %s is singular or numerically unidentified; no observations are trimmed\n",what)
        _error(498)
    }
    return(invsym(Q))
}

real colvector jwrcs_logitfit(real matrix X, real colvector a, real scalar maxiter)
{
    real colvector b, p, score, step, bn
    real matrix Q
    real scalar iter, ll, lln, scale, tries, converged
    b=J(cols(X),1,0)
    b[cols(X)]=ln(mean(a)/(1-mean(a)))
    converged=0
    for (iter=1;iter<=maxiter;iter++) {
        p=jwrcs_logistic(X*b)
        score=quadcross(X,a-p)/rows(X)
        Q=quadcross(X,p:*(1:-p),X)/rows(X)
        step=jwrcs_inverse(Q,"propensity information")*score
        if (max(abs(score))<1e-10) {
            converged=1
            break
        }
        ll=jwrcs_ll(X,a,b)
        scale=1
        for (tries=1;tries<=40;tries++) {
            bn=b+scale*step
            lln=jwrcs_ll(X,a,bn)
            if (!missing(lln) & lln>=ll-1e-13) break
            scale=scale/2
        }
        if (tries>40) break
        b=bn
    }
    if (!converged) {
        errprintf("RCS propensity fit did not reach its score tolerance; no result posted\n")
        _error(430)
    }
    return(b)
}

void jwrcs_fit(string scalar yname, string scalar onames, string scalar pnames,
    string scalar tname, string scalar gname, string scalar sample,
    real rowvector dates, real scalar gg, real scalar ss, real scalar bb,
    real scalar maxiter, string scalar sampling, string scalar bname, string scalar vname,
    string scalar cname, string scalar mname, string scalar mvname,
    string scalar diagname)
{
    real colvector y,t,g,A,cell,ix,ip,beta,gamma,m,p,w,e,rg,psi
    real matrix Xo,Xp,Q,IFbeta,IFgamma,F,L,C,means,diagnostics,V,MF
    real rowvector theta,dO,dP,pA,cn
    real scalar n,k,j,d,nj,aidx,bidx,mu,rbar,kappa,q,isg,scoremax
    y=st_data(.,yname,sample)
    t=st_data(.,tname,sample)
    g=st_data(.,gname,sample)
    n=rows(y);k=cols(dates)
    Xo=J(n,0,.)
    Xp=J(n,0,.)
    if (strlen(strtrim(onames))) Xo=st_data(.,tokens(onames),sample)
    if (strlen(strtrim(pnames))) Xp=st_data(.,tokens(pnames),sample)
    Xo=(Xo,J(n,1,1));Xp=(Xp,J(n,1,1))
    A=(g:==gg):&(t:==ss)
    q=mean(A)
    if (q==0) {
        errprintf("No target cohort/date observations\n")
        _error(2000)
    }
    theta=J(1,2*k,.)
    F=J(n,2*k,0)
    means=J(2*k,3,.)
    diagnostics=J(2*k,7,.)
    ix=selectindex(dates':==ss);aidx=ix[1]
    ix=selectindex(dates':==bb);bidx=ix[1]
    for (j=1;j<=2*k;j++) {
        isg=(j<=k)
        d=mod(j-1,k)+1
        cell=(g:==(isg*gg)) :& (t:==dates[d])
        nj=sum(cell)
        means[j,1..2]=(isg*gg,dates[d])
        diagnostics[j,1..3]=(isg*gg,dates[d],nj)
        if (j==aidx) {
            theta[j]=sum(A:*y)/sum(A)
            F[,j]=A:*(y:-theta[j])/q
            continue
        }
        if (nj<=cols(Xo)) {
            errprintf("Insufficient observations in cohort %g, date %g for outcome regression\n",isg*gg,dates[d])
            _error(2001)
        }
        ix=selectindex(cell)
        Q=quadcross(Xo,cell,Xo)/n
        beta=jwrcs_inverse(Q,"outcome-regression Gram matrix")*quadcross(Xo,cell:*y)/n
        m=Xo*beta
        e=y-m
        IFbeta=(Xo:*(cell:*e))*jwrcs_inverse(Q,"outcome-regression Gram matrix")
        ip=selectindex(A:+cell)
        gamma=jwrcs_logitfit(Xp[ip,],A[ip],maxiter)
        p=jwrcs_logistic(Xp*gamma)
        if (min(p[ip])<=1e-12 | max(p[ip])>=1-1e-12 | hasmissing(p)) {
            errprintf("RCS fitted pairwise propensity is numerically at a support boundary; no clipping or trimming is permitted\n")
            _error(498)
        }
        Q=quadcross(Xp,(A:+cell):*p:*(1:-p),Xp)/n
        IFgamma=(Xp:*((A:+cell):*(A-p)))*jwrcs_inverse(Q,"propensity information")
        w=J(n,1,0)
        w[ix]=exp(Xp[ix,]*gamma)
        if (hasmissing(w) | min(w[ix])<=0) {
            errprintf("RCS propensity density ratios are nonfinite; no clipping or trimming is permitted\n")
            _error(498)
        }
        kappa=mean(w)
        mu=sum(A:*m)/sum(A)
        rbar=sum(w:*e)/sum(w)
        theta[j]=mu+rbar
        dO=quadcross(A,Xo)/sum(A)-quadcross(w,Xo)/sum(w)
        dP=quadcross(w:*(e:-rbar),Xp)/sum(w)
        psi=A:*(m:-mu)/q+w:*(e:-rbar)/kappa+IFbeta*dO'+IFgamma*dP'
        F[,j]=psi
        scoremax=max(abs(quadcross(Xp,(A:+cell):*(A-p))/n))
        diagnostics[j,4..7]=(min(p[ip]),max(p[ip]),max(w),scoremax)
    }
    means[,3]=theta'
    L=J(k,2*k,0)
    C=(J(k,1,gg),dates')
    for (d=1;d<=k;d++) {
        if (d==bidx) continue
        L[d,d]=1
        L[d,bidx]=-1
        L[d,k+d]=-1
        L[d,k+bidx]=1
    }
    // All nuisance fits share observations. F contains FULL estimated-score
    // influences, including target-distribution and both ratio normalizers.
    if (sampling=="stratified") {
        // Condition on fixed date counts. Center the COMPLETE influence,
        // after nuisance and normalizer corrections, within each date.
        for (d=1;d<=k;d++) {
            ix=selectindex(t:==dates[d])
            F[ix,]=F[ix,]:-mean(F[ix,])
        }
    }
    else F=F:-mean(F)
    MF=F*L'
    V=quadcross(MF,MF)/(n*n)
    st_matrix(bname,theta*L')
    st_matrix(vname,(V+V')/2)
    st_matrix(cname,C)
    st_matrix(mname,means)
    st_matrix(mvname,quadcross(F,F)/(n*n))
    st_matrix(diagname,diagnostics)
}
// Completion marker is last so a partially compiled module is never reused.
string scalar jwrcs_version() return("0.9.0")
end
