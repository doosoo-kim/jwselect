// Independent finite-difference stacked-equation sandwich for RCS validation.
// No analytic nuisance derivatives from the production IF are used here.
mata:
real matrix rcs_test_m(real rowvector th, real matrix X, real colvector y,
    real colvector cell, real scalar acell, real scalar nc)
{
    real matrix M
    real colvector A,I,p,r,m,b,g
    real scalar j,o,k,mu,R
    A=(cell:==acell);k=cols(X)
    M=J(rows(X),cols(th),0)
    M[,1]=A:*(y:-th[1])
    o=2
    for(j=1;j<=nc;j++) {
        if(j==acell) continue
        I=(cell:==j)
        b=th[o..o+k-1]';g=th[o+k..o+2*k-1]'
        mu=th[o+2*k];R=th[o+2*k+1]
        p=1:/(1:+exp(-X*g));r=exp(X*g);m=X*b
        M[,o..o+k-1]=X:*(I:*(y-m))
        M[,o+k..o+2*k-1]=X:*((A+I):*(A-p))
        M[,o+2*k]=A:*(m:-mu)
        M[,o+2*k+1]=I:*r:*(y-m:-R)
        o=o+2*k+2
    }
    return(M)
}
void rcs_test_reference()
{
    real matrix X,Jac,M,IF,L,C,V
    real rowvector th,tp,tm,btrue
    real colvector y,t,g,A,I,cell,b,gam,m,w,ip
    real scalar n,k,j,o,h,q,mu,R,delta,acell
    y=st_data(.,"y");t=st_data(.,"time");g=st_data(.,"cohort")
    X=(st_data(.,"x"),J(rows(y),1,1));n=rows(y);k=cols(X)
    cell=t:+4*(g:==0);acell=4;A=(cell:==acell)
    th=J(1,1+7*(2*k+2),0)
    th[1]=sum(A:*y)/sum(A)
    L=J(4,cols(th),0);L[4,1]=1
    o=2
    for(j=1;j<=8;j++) {
        if(j==acell) continue
        I=(cell:==j)
        b=qrsolve(select(X,I),select(y,I))
        ip=selectindex(A+I)
        gam=jwrcs_logitfit(X[ip,],A[ip],100)
        m=X*b;w=I:*exp(X*gam)
        mu=sum(A:*m)/sum(A);R=sum(w:*(y-m))/sum(w)
        th[o..o+2*k+1]=(b',gam',mu,R)
        C=J(4,1,0)
        if(j==1) C[1]=1
        if(j==2) C[2]=1
        if(j==3) C=( -1\-1\0\-1 )
        if(j==5) C[1]=-1
        if(j==6) C[2]=-1
        if(j==7) C=( 1\1\0\1 )
        if(j==8) C[4]=-1
        L[,o+2*k]=C;L[,o+2*k+1]=C
        o=o+2*k+2
    }
    M=rcs_test_m(th,X,y,cell,acell,8)
    assert(max(abs(mean(M)))<1e-8)
    Jac=J(cols(th),cols(th),.)
    for(h=1;h<=cols(th);h++) {
        delta=1e-5*(1+abs(th[h]));tp=th;tm=th
        tp[h]=tp[h]+delta;tm[h]=tm[h]-delta
        Jac[,h]=(mean(rcs_test_m(tp,X,y,cell,acell,8))-
            mean(rcs_test_m(tm,X,y,cell,acell,8)))'/(2*delta)
    }
    IF=-M*luinv(Jac)'*L'
    V=quadcross(IF,IF)/(n*n)
    st_matrix("rcs_reference_b",th*L')
    st_matrix("rcs_reference_V",V)
    for(j=1;j<=4;j++) {
        ip=selectindex(t:==j)
        IF[ip,]=IF[ip,]:-mean(IF[ip,])
    }
    st_matrix("rcs_reference_V_stratified",quadcross(IF,IF)/(n*n))
    st_numscalar("rcs_reference_error",max(abs(V-st_matrix("e(V)"))))
    st_numscalar("rcs_reference_b_error",max(abs(th*L'-st_matrix("e(b)"))))
}
end
