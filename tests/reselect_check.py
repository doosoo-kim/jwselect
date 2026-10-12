#!/usr/bin/env python3
"""Independent focused-reselection checks. This does NOT execute Stata/Mata.

Uses restricted GLS solves and explicitly selects each candidate in each draw,
independently of the production candidate-loading and optimized-error routines.
Run from any directory with Python, NumPy, and SciPy installed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.stats import chi2, binom, beta

HERE = Path(__file__).resolve().parent
Q = np.array([[.04,.012,.003],[.012,.09,.015],[.003,.015,.06]])
C = np.array([[.03,.01,.004],[.015,-.025,-.002],[.01,.06,.02],[-.012,.018,.032]])
R = np.diag([.03, .08, .05, .07]) + .005
V = np.block([[R+C@np.linalg.solve(Q,C.T), C], [C.T, Q]])
Z = np.array([1.,1.5,2.,2.5,-.27,.45,.22])
P = 4
T = np.array([[.5,.5,0,0,0,0,0], [0,0,.5,.5,0,0,0]])
W = np.array([.5,0,.1,.4,0,0,0])
LOCAL = (np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]]), np.array([[0,0,0],[0,0,1]]))
JOINT = np.array([a+b for a in LOCAL[0] for b in LOCAL[1]])
MU = np.array([.07,-.11,.09])
E = np.vstack([np.zeros(7), np.random.default_rng(654221).multivariate_normal(np.zeros(7),V,15)])


def gls(target, mask, covariance=V, p=P):
    """Constrained GLS normal equations, not the production Schur formula."""
    free=np.setdiff1d(np.arange(len(covariance)),p+np.flatnonzero(mask))
    design=np.eye(len(covariance))[:,free]
    ivx=np.linalg.solve(covariance,design)
    beta=np.linalg.solve(design.T@ivx,ivx.T)
    out=target[free]@beta
    # Exact structural zeros in the post block avoid artificial tie breaking.
    out[:p]=target[:p]
    out[np.abs(out)<1e-14]=0
    return out


def fmsc(a,t,y,v=V):
    correction=t-a
    return (correction@y)**2 + t@v@t - 2*t@v@correction


def choose(target,masks,y,v=V,p=P):
    actions=np.array([gls(target,mask,v,p) for mask in masks])
    scores=np.array([fmsc(a,target,y,v) for a in actions])
    return int(np.argmin(scores)), actions


def reference_errors(eps=E,mu=MU,scope='both',v=V,z=Z,t=T,w=W,local=LOCAL,joint=JOINT,p=P):
    """Directly refit/reselect simulated samples, then subtract true target."""
    true=np.r_[z[:p],mu]
    results=[]
    chosen=[]
    for draw in eps:
        sample=true+draw
        if scope in ('cohort','both'):
            selected=[choose(t[g],local[g],sample,v,p)[0] for g in range(len(t))]
            aa=np.array([gls(t[g],local[g][selected[g]],v,p) for g in range(len(t))])
            cellrows=[]
            for g in range(len(t)):
                for j in np.flatnonzero(t[g,:p]):
                    cellrows.append((j,gls(np.eye(len(v))[j],local[g][selected[g]],v,p)))
            cells=np.array([a for j,a in sorted(cellrows)])
            actions=np.vstack([aa,w[:p]@cells])
        else:
            selected=[]
            actions=np.empty((0,len(v)))
        if scope in ('aggregate','both'):
            j,_=choose(w,joint,sample,v,p)
            if scope=='aggregate':
                actions=np.vstack([*(gls(tt,joint[j],v,p) for tt in t),gls(w,joint[j],v,p)])
            else:
                actions=np.vstack([actions,gls(w,joint[j],v,p)])
            selected.append(j)
        targets=np.vstack([t,w,w]) if scope=='both' else np.vstack([t,w])
        results.append(actions@sample-targets@true)
        chosen.append(selected)
    return np.array(results),np.array(chosen)


def check_error_identity_and_reselection():
    direct,choices=reference_errors()
    true=np.r_[Z[:P],MU]
    for r,e in enumerate(E):
        idx,_=choose(W,JOINT,true+e)
        a=gls(W,JOINT[idx])
        assert_allclose(direct[r,-1],a@e+(a-W)[P:]@MU,atol=1e-13)
    assert len(np.unique(choices[:,-1]))>1
    fixed,_=choose(W,JOINT,Z)
    a=gls(W,JOINT[fixed])
    wrongly_fixed=E@a+(a-W)[P:]@MU
    assert np.max(np.abs(wrongly_fixed-direct[:,-1]))>.01
    # Omitting nuisance bias from the selected error changes the distribution.
    wrong=[]
    for e in E:
        j,aa=choose(W,JOINT,true+e)
        wrong.append(aa[j]@e)
    assert np.max(np.abs(np.array(wrong)-direct[:,-1]))>.001


def check_custom_weights_and_joint_selection():
    direct,choices=reference_errors()
    true=np.r_[Z[:P],MU]
    wrong_report_rule=[]
    wrong_aggregate_cohort_rule=[]
    for e in E:
        sample=true+e
        total=0
        for g in range(2):
            component=np.where(T[g]!=0,W,0)
            j,aa=choose(component,LOCAL[g],sample)
            total+=aa[j]@sample-component@true
        wrong_report_rule.append(total)
        j,aa=choose(T[0],JOINT,sample)
        wrong_aggregate_cohort_rule.append(aa[j]@sample-T[0]@true)
    assert np.max(np.abs(np.array(wrong_report_rule)-direct[:,2]))>.01
    aggregate,aggregate_choices=reference_errors(scope='aggregate')
    assert_allclose(direct[:,-1],aggregate[:,-1],atol=1e-13)
    assert np.max(np.abs(np.array(wrong_aggregate_cohort_rule)-aggregate[:,0]))>.01
    # Both cohort outputs in aggregate mode inherit ONE aggregate winner.
    for r,e in enumerate(E):
        j=aggregate_choices[r,0]
        for g in range(2):
            a=gls(T[g],JOINT[j])
            assert_allclose(aggregate[r,g],a@e+(a-T[g])[P:]@MU,atol=1e-13)


def check_zero_nuisance_and_post_translation():
    direct,choices=reference_errors(mu=np.zeros(3))
    true=np.r_[Z[:P],np.zeros(3)]
    for r,e in enumerate(E):
        j,aa=choose(W,JOINT,true+e)
        assert_allclose(direct[r,-1],aa[j]@e,atol=1e-13)
    zz=Z.copy(); zz[:P]+=np.array([33,-52,71,13])
    translated,choices2=reference_errors(z=zz)
    original,choices1=reference_errors()
    assert_allclose(translated,original,atol=5e-14)
    assert np.array_equal(choices1,choices2)


def check_first_min_tie_with_different_reporting_target():
    # The average has zero h covariance, but the first-cell report does not.
    v=np.array([[3.,0.,1.],[0.,3.,-1.],[1.,-1.,1.]])
    t=np.array([.5,.5,0.]); report=np.array([1.,0.,0.])
    masks=np.array([[0],[1]])
    assert_allclose(gls(t,masks[0],v,2),gls(t,masks[1],v,2),atol=1e-15)
    for h in (-100.,-.1,0.,.1,100.):
        j,aa=choose(t,masks,np.array([1,2,h]),v,2)
        assert j==0
        assert_allclose(gls(report,masks[j],v,2),report,atol=1e-15)
    assert np.max(np.abs(gls(report,masks[1],v,2)-report))==1


def check_ellipsoid_and_candidate_envelope_bound():
    rng=np.random.default_rng(1024)
    radius=np.sqrt(chi2.ppf(.99,3)); root=np.linalg.cholesky(Q)
    aa=np.array([gls(W,mask) for mask in JOINT])
    b=(aa-W)[:,P:]
    support=radius*np.sqrt(np.einsum('ij,jk,ik->i',b,Q,b))
    center=E@aa.T+(b@Z[P:])[None,:]
    lower=np.min(center-support,axis=1); upper=np.max(center+support,axis=1)
    for _ in range(75):
        direction=rng.normal(size=3); direction*=rng.random()**(1/3)/np.linalg.norm(direction)
        mu=Z[P:]+radius*root@direction
        assert (mu-Z[P:])@np.linalg.solve(Q,mu-Z[P:])<=radius**2+1e-11
        err,_=reference_errors(mu=mu,scope='aggregate')
        assert np.all(lower<=err[:,-1]+1e-13)
        assert np.all(err[:,-1]<=upper+1e-13)
    # Standardized ball boundary/center exactly matches the stated ellipsoid.
    for x in (np.zeros(3),np.array([1.,0.,0.]),np.array([0.,-1.,0.])):
        mu=Z[P:]+radius*root@x
        assert_allclose((mu-Z[P:])@np.linalg.solve(Q,mu-Z[P:]),radius**2*(x@x),atol=1e-13)


def check_cross_cohort_covariance_and_scale():
    a=gls(T[0],LOCAL[0][-1]); b=gls(T[1],LOCAL[1][-1])
    assert abs(a@V@b)>.001
    assert_allclose((a+b)@V@(a+b),a@V@a+b@V@b+2*a@V@b)
    value,choice=reference_errors()
    factor=7.5
    scaled,choice2=reference_errors(E*factor,MU*factor,v=V*factor**2,z=Z*factor)
    assert_allclose(scaled,value*factor,atol=1e-12)
    assert np.array_equal(choice,choice2)


def mata_literal(value):
    a=np.atleast_2d(np.asarray(value,dtype=float))
    return '('+' \\ '.join(', '.join(format(v,'.17g') for v in row) for row in a)+')'


def write_stata_fixture():
    """Small cross-language fixture: no simulated RNG, no existing data changes."""
    both,_=reference_errors()
    agg,_=reference_errors(scope='aggregate')
    zero,_=reference_errors(mu=np.zeros(3))
    lines=['*! Generated by tests/reselect_check.py --write-fixture.',
           '*! Requires current jwselect Mata code loaded; leaves data/RNG unchanged.',
           'version 16.1','capture mata: mata drop jw_test_reselect_fixture()',
           'capture mata: mata drop jw_test_reselect_tie()',
           'mata:', 'void jw_test_reselect_fixture()', '{',
           'struct jw_model scalar M','struct jw_family scalar C',
           'struct jw_fit scalar F','real matrix E, expected','real colvector mu',
           'real scalar o']
    for field,value in [('z',Z[:,None]),('V',V),('postmap',[[3,3],[3,4],[4,4],[4,5]]),
                        ('hmap',[[3,0],[3,1],[4,1]]),('groups',[[3],[4]]),('bases',[[2],[3]]),
                        ('T',T),('w',W[:P])]:
        lines.append(f'M.{field} = {mata_literal(value)}')
    lines+=['M.p = 4','M.k = 3','M.nunit = 500','M.gap = 1']
    for field,value in [('masks',np.vstack(LOCAL)),('owners',[[3],[3],[3],[3],[4],[4]]),
                        ('ids',[[1],[2],[3],[4],[1],[2]])]:
        lines.append(f'C.{field} = {mata_literal(value)}')
    lines += [f'E = {mata_literal(E)}',f'mu = {mata_literal(MU[:,None])}',
              'F = jw_run(M,C,J(0,M.k,.),"focused","both",4096)',
              f'expected = {mata_literal(both)}',
              'for (o=1; o<=rows(F.A); o++) {',
              '    assert(max(abs(jw_rs_errors(M,F,o,E,mu)-expected[,o])) < 1e-10)',
              '    assert(abs(jw_rs_errors(M,F,o,E[1,.],mu)-expected[1,o]) < 1e-10)',
              '}',
              f'expected = {mata_literal(zero)}',
              'for (o=1; o<=rows(F.A); o++) {',
              '    assert(max(abs(jw_rs_errors(M,F,o,E,J(M.k,1,0))-expected[,o])) < 1e-10)',
              '}',
              'F = jw_run(M,C,J(0,M.k,.),"focused","aggregate",4096)',
              f'expected = {mata_literal(agg)}',
              'for (o=1; o<=rows(F.A); o++) {',
              '    assert(max(abs(jw_rs_errors(M,F,o,E,mu)-expected[,o])) < 1e-10)',
              '    assert(abs(jw_rs_errors(M,F,o,E[1,.],mu)-expected[1,o]) < 1e-10)',
              '}',
              '}',
              'void jw_test_reselect_tie()', '{',
              'struct jw_model scalar M','struct jw_family scalar C','struct jw_fit scalar F',
              'real matrix E','real scalar j',
              'M.z = (1 \\ 2 \\ .1)', 'M.V = (3,0,1 \\ 0,3,-1 \\ 1,-1,1)',
              'M.postmap = (3,3 \\ 3,4)', 'M.hmap = (3,1)','M.groups = 3','M.bases = 2',
              'M.T = (.5,.5,0)','M.w = (1,0)','M.p = 2','M.k = 1',
              'M.nunit = 100','M.gap = 1','C.masks = (0 \\ 1)','C.owners = (3 \\ 3)',
              'C.ids = (1 \\ 2)','E = (1,2,-100 \\ -1,1,-.1 \\ 0,0,0 \\ 1,-1,.1 \\ -2,-1,100)',
              'F = jw_run(M,C,J(0,1,.),"focused","cohort",4096)',
              '// Cohort selection scores exactly tie even though reporting actions differ.',
              'assert(max(abs(F.selA-F.selT)) == 0)',
              'for (j=1; j<=5; j++) {',
              '    assert(max(abs(jw_rs_errors(M,F,2,E,(j-3))-E[,1])) < 1e-12)',
              '}',
              '}',
              'jw_test_reselect_fixture()','jw_test_reselect_tie()',
              'mata drop jw_test_reselect_fixture()','mata drop jw_test_reselect_tie()',
              'end','display as result "jwselect deterministic reselection fixture: PASS"']
    (HERE/'reselect_smoke.do').write_text('\n'.join(lines)+'\n')


def write_search_fixture():
    n=1000; p=.0195; delta=.001/16
    ranks=np.arange(1,n+1)
    lower=ranks[binom.cdf(ranks-1,n,p)<=delta][-1]
    upper=ranks[binom.sf(ranks-1,n,p)<=delta][0]
    cp=beta.ppf(1-delta,4,997)
    lines=['capture mata: mata drop jw_test_reselect_search()', 'mata:',
           'void jw_test_reselect_search()', '{',
           'struct jw_rs_context scalar C','struct jw_rs_search_result scalar L, U',
           'real matrix starts','real colvector h','real rowvector ranks',
           'C.reportB = (1,2)','C.selectB = (0,0)','C.constant = 0',
           'C.component = 1','C.groups = 1','C.EH = J(5,2,0)',
           'C.baseline = (-2 \\ -1 \\ 0 \\ 1 \\ 2)','C.projection = J(5,1,0)',
           'h = (.2 \\ -.1)','starts = (0,1,-1 \\ 0,0,0)',
           '// Linear quantile shift: exact extrema are +/- (2+sqrt(5)).',
           'L = jw_rs_search(C,h,I(2),starts,.1,0,1,4000,1e-5)',
           'U = jw_rs_search(C,h,I(2),starts,.1,0,2,4000,1e-5)',
           'assert(abs(L.value-(-2-sqrt(5))) < 1e-4)',
           'assert(abs(U.value-(2+sqrt(5))) < 1e-4)',
           'assert(sum(L.x:^2)<=1+1e-12 & sum(U.x:^2)<=1+1e-12)',
           'assert(L.evaluations<=4000 & U.evaluations<=4000)',
           'assert(L.starts==3 & U.starts==3)',
           f'ranks = jw_rs_quantile_ranks({n},{p:.17g},{delta:.17g})',
           f'assert(ranks[1]=={lower} & ranks[2]=={upper})',
           f'assert(abs(jw_rs_probability_upper(3,1000,{delta:.17g})-{cp:.17g}) < 1e-10)',
           'ranks = jw_rs_quantile_ranks(100,.0195,.001/16)',
           'assert(ranks[1]==0)',
           '}', 'jw_test_reselect_search()', 'mata drop jw_test_reselect_search()',
           'end','display as result "jwselect numerical-search and MC-bound fixtures: PASS"']
    with (HERE/'reselect_smoke.do').open('a') as f:
        f.write('\n'+'\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-fixture',action='store_true')
    args=parser.parse_args()
    checks=[check_error_identity_and_reselection,check_custom_weights_and_joint_selection,
            check_zero_nuisance_and_post_translation,check_first_min_tie_with_different_reporting_target,
            check_ellipsoid_and_candidate_envelope_bound,check_cross_cohort_covariance_and_scale]
    for check in checks:
        check(); print('PASS',check.__name__)
    if args.write_fixture:
        write_stata_fixture(); write_search_fixture(); print('Wrote tests/reselect_smoke.do')
    print(f'PASS: {len(checks)} independent reselection checks; Stata/Mata NOT executed')

if __name__=='__main__':
    main()
