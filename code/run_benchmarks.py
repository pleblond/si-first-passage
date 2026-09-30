#!/usr/bin/env python3
"""Self-contained computations for Paper 2 (2026-09-30).

All computed tables are NEW runs, not outputs copied from the historical branch.
The swap and detour definitions independently implement the algorithms documented
by Philippe Leblond in alltoall-black-holes, branch cursor/d1-dynamics-u-2031.
SI uses shared undirected geometric edge fields, NOT independently simulated rows.
"""
from __future__ import annotations
import argparse, hashlib, json, os, platform, random, sys, time
from pathlib import Path
from typing import Any
import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import Delaunay, distance
from scipy.optimize import linprog
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
CFG = json.loads((ROOT/'code/design.json').read_text())
OUT.mkdir(exist_ok=True)

def save_json(name: str, obj: Any) -> None:
    (OUT/f'{name}.json').write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')

def grid(side: int) -> tuple[nx.Graph, np.ndarray]:
    g = nx.convert_node_labels_to_integers(nx.grid_2d_graph(side,side), ordering='sorted')
    return g, np.array([(i,j) for i in range(side) for j in range(side)],dtype=float)

def apsp(g: nx.Graph, weights: np.ndarray|None=None, edges: list|None=None) -> np.ndarray:
    if list(g.nodes()) != list(range(len(g))):
        raise ValueError('APSP requires consecutive integer node labels in order.')
    ed = list(g.edges()) if edges is None else edges
    if weights is None: weights = np.ones(len(ed))
    weights = np.asarray(weights, dtype=float)
    if len(weights)!=len(ed) or np.any(weights<=0) or not np.all(np.isfinite(weights)):
        raise ValueError('Passage costs must be finite, positive and aligned to edges.')
    u,v = np.array(ed,dtype=int).T
    adj = coo_matrix((np.r_[weights,weights],(np.r_[u,v],np.r_[v,u])),shape=(len(g),len(g))).tocsr()
    out = dijkstra(adj, directed=False)
    if not np.all(np.isfinite(out)): raise ValueError('Disconnected graph: no finite all-pairs metric.')
    return out

def mds(D: np.ndarray, X: np.ndarray|None=None) -> tuple[dict,np.ndarray,np.ndarray,np.ndarray]:
    if not np.allclose(D,D.T) or np.any(D<0) or not np.allclose(np.diag(D),0):
        raise ValueError('Expected a symmetric nonnegative dissimilarity matrix with zero diagonal.')
    Q = D*D
    B = -.5*(Q-Q.mean(0)[None,:]-Q.mean(1)[:,None]+Q.mean())
    ev,V = np.linalg.eigh(B); ev=ev[::-1]; V=V[:,::-1]
    tol = CFG['spectral_relative_tolerance'] * max(abs(ev).max(),1.)
    pos=ev[ev>tol]; neg=ev[ev < -tol]
    Y=V[:,:2]*np.sqrt(np.maximum(ev[:2],0.))[None,:]
    rec = {'gof2':float(pos[:2].sum()/pos.sum()),
           'gap23':float(ev[1]/ev[2]) if ev[2]>tol else None,
           'gap34':float(ev[2]/ev[3]) if ev[3]>tol else None,
           'negative_mass':float(-neg.sum()/pos.sum()),
           'participation_ratio':float(pos.sum()**2/np.dot(pos,pos)),
           'positive_rank':int(len(pos)), 'negative_rank':int(len(neg)),
           'lambda3_over_lambda1':float(ev[2]/ev[0])}
    embedded=distance.squareform(distance.pdist(Y))
    # Raw metric stress: no external coordinates, no post-hoc scale fitting.
    rec['mds_stress2']=float(np.linalg.norm(embedded-D)/np.linalg.norm(D))
    if X is not None:
        Xc=X-X.mean(0); Yc=Y-Y.mean(0)
        U,s,VT=np.linalg.svd(Yc.T@Xc,full_matrices=False)
        scale=s.sum()/np.sum(Yc*Yc)
        aligned=scale*Yc@(U@VT)
        rec['coordinate_nrmse']=float(np.linalg.norm(aligned-Xc)/np.linalg.norm(Xc))
    return rec,ev,V,Y

def inject_swaps(g: nx.Graph, count: int, seed: int=0) -> nx.Graph:
    rng=random.Random(seed); h=g.copy(); ed=list(h.edges()); m=len(ed)
    for _ in range(count):
        i,j=rng.sample(range(m),2); (a,b),(c,d)=ed[i],ed[j]
        if len({a,b,c,d})<4 or h.has_edge(a,d) or h.has_edge(c,b): continue
        h.remove_edge(a,b); h.remove_edge(c,d); h.add_edge(a,d); h.add_edge(c,b)
        ed[i],ed[j]=(a,d),(c,b)
    if not nx.is_connected(h): raise ValueError('Swap benchmark became disconnected.')
    return h

def detour(g: nx.Graph,u: int,v: int, radius: int=3) -> int:
    seen={u}; frontier=[u]
    for depth in range(1,radius+1):
        nxt=[]
        for x in frontier:
            for y in g[x]:
                if (x==u and y==v) or (x==v and y==u): continue
                if y==v:return depth
                if y not in seen:seen.add(y);nxt.append(y)
        frontier=nxt
    return radius+1

def event_ratio(d0: np.ndarray, d1: np.ndarray, centers: list[int]|None=None) -> dict:
    """Exact event-radius volume audit; unit vertex masses include the center."""
    centers=list(range(len(d0))) if centers is None else centers
    max_ratio=0.; witness=None
    for c in centers:
        a=np.sort(d0[c]); b=np.sort(d1[c]); radii=np.unique(np.r_[a,b])
        vals=np.searchsorted(b,radii,side='right')/np.searchsorted(a,radii,side='right')
        j=int(np.argmax(vals))
        if vals[j]>max_ratio:max_ratio=float(vals[j]);witness=[int(c),float(radii[j])]
    return {'maximum_ball_ratio':max_ratio,'witness_center_radius':witness,'centers_audited':len(centers),
            'contracted_unordered_pairs':int(np.count_nonzero(np.triu(d1<d0-1e-10,1)))}

def do_static() -> None:
    cfg=CFG['damage']; base,X=grid(cfg['side']); D0=apsp(base); out=[]
    mats={'reference_hop':D0,'coordinates':X}
    for swaps in cfg['attempted_swaps']:
        h=inject_swaps(base,swaps,cfg['seed']); ed=list(h.edges())
        mask=np.array([detour(h,u,v,cfg['long_detour_radius'])>cfg['long_detour_radius'] for u,v in ed])
        cost_list=cfg['long_costs'] if swaps==10 else [1]
        for cost in cost_list:
            w=np.where(mask,cost,1.); D=apsp(h,w,ed); rec,ev,V,Y=mds(D,X)
            rec.update(swaps=swaps,long_count=int(mask.sum()),long_cost=cost,
                       underpriced_edges=int(sum(w[j]<D0[u,v]-1e-10 for j,(u,v) in enumerate(ed))))
            rec.update(event_ratio(D0,D));out.append(rec)
            mats[f'D_swaps{swaps}_cost{cost}']=D
            mats[f'eigen_swaps{swaps}_cost{cost}']=ev
        if swaps==10:
            np.savetxt(OUT/'damage_edges.csv',np.c_[ed,mask.astype(int)],delimiter=',',header='u,v,is_long',comments='',fmt='%d')
    shells=[]
    for r in CFG['shell_radii']:
        P=np.array([(x,y) for x in range(-r,r+1) for y in range(-r,r+1) if abs(x)+abs(y)==r])
        D=distance.cdist(P,P,metric='cityblock'); rec,*_=mds(D,P)
        rec.update(radius=r,shell_vertices=len(P),ball_vertices=1+2*r*(r+1),cut_edges=8*r+4)
        shells.append(rec)
    save_json('static',out);save_json('shells',shells)
    np.savez_compressed(OUT/'static_matrices.npz',**mats)
    print('Static and shell diagnostics complete.',flush=True)

def geometric_fields(g: nx.Graph, probs: np.ndarray, K: int, rng: np.random.Generator,
                     edges: list|None=None) -> np.ndarray:
    edges=list(g.edges()) if edges is None else edges
    if np.any(probs<=0) or np.any(probs>1):raise ValueError('Probabilities must lie in (0,1].')
    acc=np.zeros((len(g),len(g)))
    for _ in range(K):
        acc += apsp(g,rng.geometric(probs),edges)
    return acc/K

def do_si() -> None:
    cfg=CFG['si']; g,X=grid(cfg['side']); D0=apsp(g); out=[]; mats={'hop':D0,'coordinates':X}
    rec,*_=mds(D0,X);rec.update(method='hop',seed=None,trials=0);out.append(rec)
    for seed in cfg['seeds']:
        D=geometric_fields(g,np.full(g.number_of_edges(),cfg['probability']),cfg['trials_per_replication'],np.random.default_rng(seed))
        rec,ev,V,Y=mds(D,X);rec.update(method='shared_field_SI',seed=seed,trials=cfg['trials_per_replication'])
        maxviol=max(float(np.max(D-D[:,k,None]-D[None,k,:])) for k in range(len(g)))
        assert maxviol < 1e-10
        rec['maximum_triangle_violation']=maxviol
        out.append(rec);mats[f'mean_seed{seed}']=D;mats[f'eigen_seed{seed}']=ev
        print(f'SI baseline seed {seed}: GoF2={rec["gof2"]:.4f}',flush=True)
    save_json('si_baseline',out);np.savez_compressed(OUT/'si_matrices.npz',**mats)
    cfg=CFG['weak_shortcuts'];g,X=grid(cfg['side']);Dref=apsp(g);base_edges=list(g.edges());rng=np.random.default_rng(cfg['graph_seed'])
    longs=[]
    while len(longs)<cfg['number']:
        u,v=sorted(rng.choice(len(g),2,replace=False).tolist())
        if not g.has_edge(u,v) and (u,v) not in longs and Dref[u,v]>=cfg['minimum_reference_span']:longs.append((u,v))
    h=g.copy();h.add_edges_from(longs);all_edges=base_edges+longs
    costs=cfg['cost_factors'];sums=[np.zeros((len(g),len(g))) for _ in costs];baseline=np.zeros((len(g),len(g)))
    rng=np.random.default_rng(cfg['field_seed'])
    for trial in range(cfg['trials']):
        u=rng.random(len(all_edges)); u=np.maximum(u,np.finfo(float).tiny)
        base_tau=np.ceil(np.log(u[:len(base_edges)])/np.log1p(-cfg['base_probability']))
        D0=apsp(g,base_tau,base_edges);baseline+=D0;last=np.zeros_like(D0)
        for k,cost in enumerate(costs):
            tau_long=np.ceil(np.log(u[len(base_edges):])/np.log1p(-cfg['base_probability']/cost))
            D=apsp(h,np.r_[base_tau,tau_long],all_edges)
            assert np.all(D>=last-1e-10) and np.all(D<=D0+1e-10)
            sums[k]+=D;last=D
    baseline/=cfg['trials'];out=[];mats={'baseline_mean':baseline,'coordinates':X,'shortcut_edges':np.array(longs)}
    brec,*_=mds(baseline,X);brec.update(cost_factor=None,long_probability=0,relative_metric_error=0.,reference=True);out.append(brec)
    for cost,S in zip(costs,sums):
        D=S/cfg['trials'];rec,*_=mds(D,X)
        rec.update(cost_factor=cost,long_probability=cfg['base_probability']/cost,reference=False,
                   relative_metric_error=float(np.linalg.norm(D-baseline)/np.linalg.norm(baseline)))
        out.append(rec);mats[f'mean_factor{cost}']=D
    save_json('weak_shortcuts',out);np.savez_compressed(OUT/'weak_shortcuts.npz',**mats)
    print('Coupled weak-shortcut benchmark complete.',flush=True)

def do_delaunay() -> None:
    cfg=CFG['delaunay'];out=[];mats={}
    for seed in cfg['seeds']:
        X=np.random.default_rng(seed).random((cfg['n'],2))*cfg['box_side'];tri=Delaunay(X)
        edges=sorted(set(tuple(sorted((int(t[i]),int(t[(i+1)%3])))) for t in tri.simplices for i in range(3)))
        g=nx.Graph();g.add_nodes_from(range(len(X)));g.add_edges_from(edges)
        E=distance.squareform(distance.pdist(X));rho2=np.sum((X-cfg['box_side']/2)**2,axis=1)
        for method,w in [('hop',np.ones(len(edges))),('Euclidean_edge_oracle',np.array([E[u,v] for u,v in edges]))]:
            D=apsp(g,w,edges);rec,ev,V,Y=mds(D,X)
            tort=(D/np.where(E>0,E,np.inf)).sum(1)/(len(X)-1)
            def corr(a,b):return float(abs(np.corrcoef(a,b)[0,1]))
            # Euclidean scale is fitted ONLY for diagnostic residual evaluation.
            scale=float(np.sum(D*E)/np.sum(D*D)); residual=scale*D-E
            pairerr=np.mean(residual,axis=1)
            bands=np.quantile(rho2,[0,.25,.5,.75,1]);bandmean=[]
            for k in range(4):
                mask=(rho2>=bands[k]) & ((rho2<bands[k+1]) if k<3 else (rho2<=bands[k+1]))
                bandmean.append(float(pairerr[mask].mean()/E.mean()))
            rec.update(seed=seed,method=method,edges=len(edges),radial_correlation=corr(V[:,2],rho2),
                       tortuosity_correlation=corr(V[:,2],tort),
                       pair_distance_nrmse=float(np.linalg.norm(residual)/np.linalg.norm(E)),
                       radial_quartile_mean_residuals=bandmean)
            out.append(rec)
            if seed==0:
                mats[method+'_D']=D;mats[method+'_eigen']=ev;mats[method+'_v3']=V[:,2];mats[method+'_tortuosity']=tort
        if seed==0:mats['coordinates']=X;mats['rho2']=rho2;mats['edges']=np.array(edges)
        print(f'Delaunay seed {seed} complete.',flush=True)
    save_json('delaunay',out);np.savez_compressed(OUT/'delaunay_matrices.npz',**mats)

def coarsen(g: nx.Graph,rule: str) -> nx.Graph:
    h=nx.Graph();h.add_nodes_from(sorted(set((u[0]//2,u[1]//2) for u in g)))
    bins={}
    for u,v,a in g.edges(data=True):
        x=(u[0]//2,u[1]//2);y=(v[0]//2,v[1]//2)
        if x==y:continue
        key=tuple(sorted((x,y)));bins.setdefault(key,[]).append(a['cost'])
    for (x,y),vals in sorted(bins.items()):h.add_edge(x,y,cost=min(vals) if rule=='min' else float(np.mean(vals)))
    return h

def do_rg() -> None:
    cfg=CFG['rg'];base=nx.grid_2d_graph(cfg['side'],cfg['side']);dam=base.copy()
    nx.connected_double_edge_swap(dam,cfg['swaps'],seed=cfg['seed'])
    # Persist the realized graph: connected_double_edge_swap is version-sensitive.
    np.savetxt(OUT/'rg_initial_edges.csv',np.array([[*u,*v] for u,v in dam.edges()]),delimiter=',',header='ux,uy,vx,vy',comments='',fmt='%d')
    out=[]
    for price in cfg['long_costs']:
        for rule in cfg['rules']:
            h=dam.copy()
            for u,v in h.edges():h[u][v]['cost']=price if sum(abs(u[k]-v[k]) for k in (0,1))>1 else 1
            for lev in range(cfg['levels']):
                spans=np.array([sum(abs(u[k]-v[k]) for k in (0,1)) for u,v in h.edges()]);w=np.array([a['cost'] for _,_,a in h.edges(data=True)])
                mask=spans>1;viol=w<spans-1e-10;fabric=~mask
                out.append(dict(price=price,rule=rule,level=lev,nodes=len(h),edges=h.number_of_edges(),
                    long_count=int(mask.sum()),underpriced_long_count=int((viol&mask).sum()),
                    long_fraction=float(mask.mean()),underpriced_fraction=float(viol.mean()),
                    phi=float(viol[mask].mean()) if mask.any() else 0.,
                    mean_long_cost_span=float(np.mean(w[mask]/spans[mask])) if mask.any() else None,
                    fabric_unit_fraction=float(np.mean(np.isclose(w[fabric],1))) if fabric.any() else None,
                    connected=nx.is_connected(h)))
                h=coarsen(h,rule)
    save_json('rg',out);print('Weighted coarsening benchmark complete.',flush=True)

def ollivier_curvature(g: nx.Graph,D: np.ndarray,u: int,v: int) -> float:
    nu=list(g[u]);nv=list(g[v]);m,n=len(nu),len(nv);C=D[np.ix_(nu,nv)]
    A=np.zeros((m+n,m*n))
    for i in range(m):A[i,i*n:(i+1)*n]=1
    for j in range(n):A[m+j,j::n]=1
    result=linprog(C.ravel(),A_eq=A,b_eq=np.r_[np.full(m,1/m),np.full(n,1/n)],bounds=(0,None),method='highs')
    if not result.success:raise RuntimeError(result.message)
    return float(1-result.fun/D[u,v])

def do_disk() -> None:
    cfg=CFG['disk'];g,X=grid(cfg['side']);center=np.array([cfg['side']//2]*2)
    plug=set(map(int,np.flatnonzero(np.sum((X-center)**2,axis=1)<=cfg['radius']**2)));out=[]
    for name in ['clique_disk','king_disk']:
        h=g.copy();P=sorted(plug)
        for i,u in enumerate(P):
            for v in P[i+1:]:
                if name=='clique_disk' or np.max(np.abs(X[u]-X[v]))<=1:h.add_edge(u,v)
        D=apsp(h);classes={'internal':[],'boundary':[],'fabric':[]}
        for u,v in h.edges():
            key='internal' if u in plug and v in plug else 'boundary' if (u in plug)!=(v in plug) else 'fabric'
            classes[key].append((u,v))
        rng=random.Random(cfg['seed'])
        for key,ed in classes.items():
            sample=rng.sample(ed,min(cfg['edges_per_class'],len(ed)))
            values=[ollivier_curvature(h,D,u,v) for u,v in sample]
            out.append(dict(graph=name,edge_class=key,plug_vertices=len(plug),sample_count=len(sample),
                            class_size=len(ed),mean=float(np.mean(values)),minimum=min(values),maximum=max(values),
                            edges=sample,curvatures=values))
        print(f'Curvature {name} complete.',flush=True)
    save_json('disk_curvature',out)

def do_checks() -> None:
    checks=[]
    # Complete graphs yield a high-rank equilateral simplex, not a 2D state.
    for n in [8,16,64]:
        D=np.ones((n,n))-np.eye(n);rec,ev,*_=mds(D)
        assert rec['positive_rank']==n-1 and np.allclose(ev[:n-1],.5)
        assert np.isclose(rec['gof2'],2/(n-1))
        checks.append({'check':'complete_graph','n':n,'gof2':rec['gof2'],'rank':rec['positive_rank']})
    # Rescaled binary paths have exact 1D Euclidean realizations.
    n=51;P=np.arange(n)/(n-1);D=np.abs(P[:,None]-P[None,:]);rec,ev,*_=mds(D)
    assert rec['positive_rank']==1 and rec['negative_rank']==0
    checks.append({'check':'rescaled_path','n':n,'positive_rank':rec['positive_rank']})
    # Geometric triangle at p=1/2: E min(X,Y+Z)=16/9, while min sum E tau=2.
    rng=np.random.default_rng(12345);tau=rng.geometric(.5,size=(200000,3));T=np.minimum(tau[:,0],tau[:,1]+tau[:,2])
    checks.append({'check':'geometric_triangle','analytic_mean':16/9,'mc_mean':float(T.mean()),
                   'mc_standard_error':float(T.std(ddof=1)/np.sqrt(len(T))),'trials':len(T)})
    assert abs(T.mean()-16/9)<6*T.std(ddof=1)/np.sqrt(len(T))
    # Shell counts and sum-of-balls multiplicities.
    for R in range(1,31):
        actual=sum(1+2*r*(r+1) for r in range(R+1));expected=(R+1)+2*R*(R+1)*(R+2)//3
        assert actual==expected
    # Samplewise edge price domination and intact-backbone equality.
    g,X=grid(8);D0=apsp(g);h=g.copy();h.add_edges_from([(0,63),(7,56),(3,59)])
    weights=np.array([D0[u,v] for u,v in h.edges()]);Dw=apsp(h,weights)
    assert np.array_equal(D0,Dw)
    checks.append({'check':'intact_backbone_exact_recovery','maximum_error':float(np.max(abs(D0-Dw)))})
    # A positive spanning-tree SI metric need not be Euclidean in the plane.
    s=nx.star_graph(8);D=2*apsp(s);rec,*_=mds(D)
    checks.append({'check':'star_mean_SI_p_half','gof2':rec['gof2'],'positive_rank':rec['positive_rank'],'negative_mass':rec['negative_mass']})
    save_json('verification',checks)
    env={'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'networkx':nx.__version__,
         'platform':platform.platform(),'design_sha256':hashlib.sha256((ROOT/'code/design.json').read_bytes()).hexdigest()}
    save_json('environment',env);print('Exact checks complete.',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--part',choices=['all','static','si','delaunay','rg','disk','checks'],default='all');args=parser.parse_args()
    jobs={'checks':do_checks,'static':do_static,'si':do_si,'delaunay':do_delaunay,'rg':do_rg,'disk':do_disk}
    for key,fn in jobs.items():
        if args.part in ['all',key]:
            start=time.monotonic();fn();print(f'{key}: {time.monotonic()-start:.2f} seconds',flush=True)
