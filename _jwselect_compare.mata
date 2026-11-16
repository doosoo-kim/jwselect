*! jwselect 0.9.1 16nov2026
version 16.1
mata:
/* A reporting-only decomposition. No selection is rerun and no e() result
   is replaced. Each target is split into its cohort components for own MD.
   Their loadings are added BEFORE applying the full covariance, retaining
   all cross-cohort covariance in aggregate standard errors. */
void jw_compare_loadings(real matrix V, real matrix T, real matrix U,
                         real matrix postmap, real matrix hmap,
                         real matrix Aown, real matrix Ajoint)
{
    real scalar p, k, i, g
    real colvector groups, ii
    real rowvector component, ownmask
    p=rows(postmap)
    k=rows(hmap)
    if (cols(T)!=p+k | rows(V)!=p+k | cols(V)!=p+k |
        rows(U)!=rows(T) | cols(U)!=k) {
        jw_core_error("incompatible dimensions in estimator comparison")
    }
    if (k>0) {
        if (max(abs(T[.,(p+1)..(p+k)]))>0) {
            jw_core_error("estimator comparison requires targets defined on post-treatment cells")
        }
    }
    groups=uniqrows(sort(postmap[.,1],1))
    Aown=J(rows(T),cols(T),0)
    Ajoint=J(rows(T),cols(T),0)
    for (i=1; i<=rows(T); i++) {
        Ajoint[i,.]=jw_loading(V,T[i,.],U[i,.],k)
        for (g=1; g<=rows(groups); g++) {
            ii=selectindex(postmap[.,1]:==groups[g])
            component=J(1,cols(T),0)
            component[ii]=T[i,ii]
            if (k>0) ownmask=U[i,.]:*(hmap[.,1]':==groups[g])
            else ownmask=J(1,0,0)
            Aown[i,.]=Aown[i,.]+jw_loading(V,component,ownmask,k)
        }
    }
}

void jw_compare_put(string scalar localname, real matrix X,
                    string matrix rs, string matrix cs)
{
    string scalar name
    name=st_local(localname)
    st_matrix(name,X)
    st_matrixrowstripe(name,rs)
    st_matrixcolstripe(name,cs)
}

void jw_compare_stata()
{
    real scalar p, k, i, output, cells
    real matrix V, T, U, H, post, map, cs, S, Aown, Ajoint,
                Vbase, Vown, Vjoint, tab
    real rowvector primary
    real colvector z, bbase, bown, bjoint, dbase, down, djoint
    string matrix rs, zs, hs, maps, ts
    cells=st_local("cells")!=""
    output=strtoreal(st_local("output"))
    z=st_matrix("e(z)")'
    V=st_matrix("e(V_z)")
    post=st_matrix("e(postmap)")
    p=rows(post)
    k=st_numscalar("e(K)")
    H=J(0,2,.)
    primary=J(1,k,0)
    if (k>0) {
        H=st_matrix("e(hmap)")
        cs=st_matrix("e(cell_selection)")
        primary=colsum(cs):>0
    }
    if (cells) {
        T=(I(p),J(p,k,0))
        map=post
        rs=st_matrixrowstripe("e(postmap)")
        maps=st_matrixcolstripe("e(postmap)")
        U=J(p,1,1)*primary
    }
    else {
        T=st_matrix("e(targets)")
        map=st_matrix("e(output_map)")
        rs=st_matrixrowstripe("e(targets)")
        maps=st_matrixcolstripe("e(output_map)")
        U=J(rows(T),1,1)*primary
        /* Kind 3 is the independently aggregate-focused target. Other
           primary rows share the union, including cohort-focused rows. */
        if (k>0) {
            S=st_matrix("e(selection)")
            for (i=1; i<=rows(T); i++) {
                if (map[i,1]==3) U[i,.]=S[i,.]
            }
        }
        if (output>0) {
            T=T[output,.]
            U=U[output,.]
            map=map[output,.]
            rs=rs[output,.]
        }
    }
    jw_compare_loadings(V,T,U,post,H,Aown,Ajoint)
    Vbase=T*V*T'
    Vown=Aown*V*Aown'
    Vjoint=Ajoint*V*Ajoint'
    bbase=T*z
    bown=Aown*z
    bjoint=Ajoint*z
    dbase=diagonal(Vbase)
    down=diagonal(Vown)
    djoint=diagonal(Vjoint)
    if (min((dbase\down\djoint)) < -1e-10*max((1\abs(diagonal(V))))) {
        jw_core_error("negative variance in estimator comparison")
    }
    /* Remove only floating-point negative zeros. */
    dbase=rowmax((dbase,J(rows(T),1,0)))
    down=rowmax((down,J(rows(T),1,0)))
    djoint=rowmax((djoint,J(rows(T),1,0)))
    tab=(bbase,sqrt(dbase),bown,sqrt(down),bjoint,sqrt(djoint),
         bown-bbase,bjoint-bbase,bjoint-bown)
    ts=(J(9,1,""),("baseline"\"se_baseline"\"own"\"se_own"\
        "joint"\"se_joint"\"own_minus_base"\"joint_minus_base"\"joint_minus_own"))
    zs=st_matrixcolstripe("e(z)")
    jw_compare_put("table",tab,rs,ts)
    jw_compare_put("Vbaseline",Vbase,rs,rs)
    jw_compare_put("Vown",Vown,rs,rs)
    jw_compare_put("Vjoint",Vjoint,rs,rs)
    jw_compare_put("Abaseline",T,rs,zs)
    jw_compare_put("Aown",Aown,rs,zs)
    jw_compare_put("Ajoint",Ajoint,rs,zs)
    jw_compare_put("map",map,rs,maps)
    if (k>0) {
        hs=st_matrixrowstripe("e(hmap)")
        jw_compare_put("mask",U,rs,hs)
    }
}
end
