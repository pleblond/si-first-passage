#!/usr/bin/env python3
"""Secondary diagnostic batch; chronology disclosed in extension_design.json."""
from __future__ import annotations
import json, hashlib, os, time
import numpy as np
import networkx as nx
from scipy.spatial import Delaunay, distance
from run_benchmarks import ROOT,OUT,apsp,mds,geometric_fields,save_json
cfg=json.loads((ROOT/'code/extension_design.json').read_text())
(OUT/'extension_design.sha256').write_text(hashlib.sha256((ROOT/'code/extension_design.json').read_bytes()).hexdigest()+'\n')
dc=cfg['dimension_controls'];out=[]
for dim_str,side in (dc['side_by_dimension'].items() if not os.environ.get('ONLY_DELAUNAY') else []):
    dim=int(dim_str);original=nx.grid_graph(dim=[side]*dim)
    X=np.array(sorted(original),dtype=float).reshape(-1,dim)
    g=nx.convert_node_labels_to_integers(original,ordering='sorted');hop=apsp(g)
    for method,D in [('hop',hop),('SI',geometric_fields(g,np.full(g.number_of_edges(),dc['probability']),dc['trials'],np.random.default_rng(dc['field_seed'])))]:
        rec,ev,V,Y=mds(D);tol=1e-10*max(abs(ev).max(),1);pos=ev[ev>tol]
        Y=V[:,:dim]*np.sqrt(np.maximum(ev[:dim],0))[None,:];Xc=X-X.mean(0);Yc=Y-Y.mean(0)
        U,s,VT=np.linalg.svd(Yc.T@Xc,full_matrices=False);aligned=s.sum()/np.sum(Yc**2)*Yc@(U@VT)
        rec.update(dimension=dim,side=side,n=len(g),method=method,gof_d=float(pos[:dim].sum()/pos.sum()),
                   dimension_gap=float(ev[dim-1]/ev[dim]) if ev[dim]>tol else None,
                   coordinate_nrmse_d=float(np.linalg.norm(aligned-Xc)/np.linalg.norm(Xc)))
        out.append(rec)
    print(f'Dimension control d={dim} complete.',flush=True)
if out: save_json('dimension_controls',out)
cfg=cfg['delaunay_SI'];out=[];mats={}
if os.environ.get('GRAPH_SEED') and (OUT/'delaunay_si.json').exists(): out=json.loads((OUT/'delaunay_si.json').read_text())
for seed in cfg['graph_seeds']:
    if os.environ.get('GRAPH_SEED') is not None and seed!=int(os.environ['GRAPH_SEED']): continue
    start=time.monotonic()
    X=np.random.default_rng(seed).random((cfg['n'],2))*cfg['box_side'];tri=Delaunay(X)
    edges=sorted(set(tuple(sorted((int(t[i]),int(t[(i+1)%3])))) for t in tri.simplices for i in range(3)))
    g=nx.Graph();g.add_nodes_from(range(len(X)));g.add_edges_from(edges)
    D=geometric_fields(g,np.full(len(edges),cfg['probability']),cfg['trials'],np.random.default_rng(seed+cfg['field_seed_offset']),edges)
    E=distance.squareform(distance.pdist(X));r2=np.sum((X-cfg['box_side']/2)**2,axis=1)
    rec,ev,V,Y=mds(D,X);tort=(D/np.where(E>0,E,np.inf)).sum(1)/(len(X)-1)
    scale=float(np.sum(D*E)/np.sum(D*D));residual=scale*D-E;pairerr=residual.mean(1);bands=np.quantile(r2,[0,.25,.5,.75,1]);bm=[]
    for k in range(4):
        mask=(r2>=bands[k]) & ((r2<bands[k+1]) if k<3 else (r2<=bands[k+1]));bm.append(float(pairerr[mask].mean()/E.mean()))
    rec.update(seed=seed,method='shared_field_SI',trials=cfg['trials'],radial_correlation=float(abs(np.corrcoef(V[:,2],r2)[0,1])),
               tortuosity_correlation=float(abs(np.corrcoef(V[:,2],tort)[0,1])),
               pair_distance_nrmse=float(np.linalg.norm(residual)/np.linalg.norm(E)),radial_quartile_mean_residuals=bm)
    out=[r for r in out if r['seed']!=seed]+[rec]
    save_json('delaunay_si', sorted(out,key=lambda r:r['seed']))
    if seed==0:mats={'coordinates':X,'SI_D':D,'SI_eigen':ev,'SI_v3':V[:,2],'rho2':r2}
    if seed==0: np.savez_compressed(OUT/'delaunay_si_matrices.npz',**mats)
    print(f'{time.monotonic()-start:.2f}s Delaunay SI seed {seed}: GoF2={rec["gof2"]:.4f}, radial corr={rec["radial_correlation"]:.4f}',flush=True)
save_json('delaunay_si',sorted(out,key=lambda r:r['seed']))
