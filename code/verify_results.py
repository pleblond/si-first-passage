#!/usr/bin/env python3
"""Fast integrity/numerical guards on the supplied result files."""
from pathlib import Path
import json
import numpy as np
from run_benchmarks import mds, ROOT
R=ROOT/'results'
def load(name):return json.loads((R/(name+'.json')).read_text())
checks=[]
assert len(load('si_baseline'))==6
assert {r['seed'] for r in load('delaunay_si')}==set(range(5))
assert len(load('delaunay'))==10
assert len(load('rg'))==24
assert len(load('disk_curvature'))==6
assert len(load('dimension_controls'))==6
checks.append('All configured result groups and seeds are present.')
A=np.array([[0.,1.,4.],[1.,0.,1.],[4.,1.,0.]])
assert np.array_equal(A,A.T) and A[0,2]>A[0,1]+A[1,2]
checks.append('Independent-row counterexample is symmetric and violates a triangle inequality.')
M=np.load(R/'si_matrices.npz')
for seed in range(5):
 D=M[f'mean_seed{seed}']
 assert np.array_equal(D,D.T) and np.all(np.diag(D)==0)
 assert max(np.max(D-D[:,k,None]-D[None,k,:]) for k in range(len(D)))==0
checks.append('All five square-grid SI mean matrices satisfy every triangle inequality.')
W=np.load(R/'weak_shortcuts.npz');base=W['baseline_mean'];last=np.zeros_like(base)
for a in [1,4,20,100]:
 D=W[f'mean_factor{a}'];assert np.all(D>=last) and np.all(D<=base);last=D
checks.append('Weak-shortcut stored means are pointwise monotone and bounded by the coupled intact mean.')
# General finite-matrix perturbation inequality, using a convex mixture of Euclidean metrics.
rng=np.random.default_rng(77);X=rng.normal(size=(15,2));Y=rng.normal(size=(15,3))
D0=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=2);D1=np.linalg.norm(Y[:,None,:]-Y[None,:,:],axis=2);D=.98*D0+.02*D1
J=np.eye(len(D))-np.ones_like(D)/len(D);diff=-.5*J@(D**2-D0**2)@J
delta=np.max(abs(D-D0));bound=len(D)/2*delta*(2*D0.max()+delta)
assert np.linalg.norm(diff,ord=2)<=bound+1e-12
checks.append('The explicit Gram perturbation bound holds on a nontrivial metric mixture.')
# All generated LaTeX table files exist and are nonempty.
for name in ['static','shells','si','weak','rg','delaunay','dimensions','disk']:
 text=(ROOT/'tables'/f'{name}.tex').read_text();assert '\\begin{tabular}' in text and '\\bottomrule' in text
checks.append('All eight data tables are present.')
(R/'final_checks.json').write_text(json.dumps({'checks_passed':checks,'gram_observed_error':float(np.linalg.norm(diff,ord=2)),'gram_bound':float(bound)},indent=2)+'\n')
for c in checks:print('PASS:',c)
