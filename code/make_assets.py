#!/usr/bin/env python3
"""Render figures and tables directly from saved results (no recomputation)."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=ROOT/'figures';T=ROOT/'tables'
F.mkdir(exist_ok=True);T.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':10,'legend.fontsize':9,'pdf.fonttype':42})
def load(n):return json.loads((R/(n+'.json')).read_text())
def finish(name,fig):
    fig.tight_layout();fig.savefig(F/(name+'.pdf'),bbox_inches='tight');fig.savefig(F/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig)
def tab(name,header,rows,spec):
    s='\\begin{tabular}{'+spec+'}\n\\toprule\n'+' & '.join(header)+r' \\'+'\n\\midrule\n'
    s+='\n'.join(' & '.join(map(str,row))+r' \\' for row in rows)
    s+='\n\\bottomrule\n\\end{tabular}\n';(T/(name+'.tex')).write_text(s)
def f(v,k=3):return '--' if v is None else f'{v:.{k}f}'
static=load('static');shells=load('shells');si=load('si_baseline');weak=load('weak_shortcuts');rg=load('rg')
delrows=load('delaunay');delsi=load('delaunay_si');dims=load('dimension_controls');disk=load('disk_curvature')
tab('static',['Swaps',r'$\ell_{\rm long}$','Longs',r'$G_2$',r'$g_{23}$',r'$r_+$',r'$\max R$','Bad edges'],
    [[r['swaps'],r['long_cost'],r['long_count'],f(r['gof2']),f(r['gap23'],2),f(r['participation_ratio'],2),f(r['maximum_ball_ratio'],2),r['underpriced_edges']] for r in static], 'rrrrrrrr')
tab('shells',[r'$r$',r'$|S_r|$',r'$|B_r|$',r'$G_2$',r'$n_+$',r'$r_+$'],
    [[r['radius'],r['shell_vertices'],r['ball_vertices'],f(r['gof2']),r['positive_rank'],f(r['participation_ratio'])] for r in shells],'rrrrrr')
fields=['gof2','gap23','negative_mass','mds_stress2','coordinate_nrmse'];samples=[r for r in si if r['method']=='shared_field_SI']
means={k:float(np.mean([r[k] for r in samples])) for k in fields};sds={k:float(np.std([r[k] for r in samples],ddof=1)) for k in fields}
tab('si',['Metric',r'$G_2$',r'$g_{23}$',r'$\eta_-$',r'$S_2$',r'$E_X$'],
    [['Hop']+[f(si[0][k],4 if k!='gap23' else 2) for k in fields],['SI mean']+[f(means[k],4 if k!='gap23' else 2) for k in fields],['SI sample SD']+[f(sds[k],4 if k!='gap23' else 2) for k in fields]],'lrrrrr')
tab('weak',[r'$a$',r'$p_{\rm long}$',r'$E_D$',r'$G_2$',r'$g_{23}$',r'$E_X$'],
    [[r['cost_factor'] if r['cost_factor'] else r'$\infty$',f(r['long_probability'],3),f(r['relative_metric_error'],4),f(r['gof2'],4),f(r['gap23'],2),f(r['coordinate_nrmse'],4)] for r in weak[1:]+weak[:1]],'rrrrrr')
tab('rg',[r'Level','$N$','$|E|$',r'$\lambda_{\rm long}$',r'$\Phi$',r'$\langle\ell/s\rangle_{\rm long}$'],
    [[r['level'],r['nodes'],r['edges'],f(r['long_fraction']),f(r['phi']),f(r['mean_long_cost_span'])] for r in rg if r['price']==10 and r['rule']=='min'],'rrrrrr')
rows=[];summaries={}
for meth,label,data in [('hop','Hop',delrows),('shared_field_SI','SI',delsi),('Euclidean_edge_oracle','Length oracle',delrows)]:
    rr=[r for r in data if r['method']==meth];keys=['gof2','gap34','radial_correlation','pair_distance_nrmse','coordinate_nrmse'];summaries[meth]={k:float(np.mean([r[k] for r in rr])) for k in keys}
    rows.append([label]+[f(summaries[meth][k]) for k in keys])
tab('delaunay',['Metric',r'$G_2$',r'$g_{34}$',r'$|\rho(v_3,r^2)|$',r'$E_{\rm pair}$',r'$E_X$'],rows,'lrrrrr')
tab('dimensions',['$d$','$N$','Metric',r'$G_d$',r'$g_{d,d+1}$',r'$E_X^{(d)}$'],
    [[r['dimension'],r['n'],r['method'],f(r['gof_d'],4),f(r['dimension_gap'],2) if r['dimension_gap'] is not None else r'$\infty$',f(r['coordinate_nrmse_d'],4)] for r in dims],'rrlrrr')
tab('disk',['Plug','Edge class','$n$',r'Mean $\kappa$',r'Minimum',r'Maximum'],
    [[{'clique_disk':'Clique','king_disk':'King'}[r['graph']],r['edge_class'],r['sample_count'],f(r['mean']),f(r['minimum']),f(r['maximum'])] for r in disk],'llrrrr')
# Figure 1: damage destroys concentration rather than creating rank two.
a=np.load(R/'static_matrices.npz');fig,ax=plt.subplots(figsize=(6.3,3.35))
for key,label in [('eigen_swaps0_cost1','Healthy grid'),('eigen_swaps10_cost1','10 swaps, unit cost'),('eigen_swaps10_cost20','10 swaps, long cost 20')]:
    e=a[key];pos=e[e>1e-10*max(abs(e))];ax.plot(np.arange(1,21),pos[:20]/pos.sum(),marker='o',markersize=3,label=label)
ax.set(xlabel='Positive eigenvalue index',ylabel='Fraction of positive spectral mass',yscale='log',xticks=[1,2,3,5,10,15,20]);ax.legend(frameon=False);finish('damage_spectrum',fig)
# Figure 2: SI square-grid calibration.
a=np.load(R/'si_matrices.npz');Dhop=a['hop'];Q=Dhop**2;B=-.5*(Q-Q.mean(0)[None,:]-Q.mean(1)[:,None]+Q.mean());hop=np.linalg.eigvalsh(B)[::-1]
fig,ax=plt.subplots(figsize=(6.3,3.35))
for ev,label in [(hop,'Hop metric'),(a['eigen_seed0'],'Mean SI, 256 shared fields')]:
    pos=ev[ev>1e-10*max(abs(ev))];ax.plot(np.arange(1,16),pos[:15]/pos.sum(),marker='o',markersize=3,label=label)
ax.set(xlabel='Positive eigenvalue index',ylabel='Fraction of positive spectral mass',yscale='log',xticks=[1,2,3,5,10,15]);ax.legend(frameon=False);finish('si_spectrum',fig)
# Figure 3: explicit weak-shortcut limit.
fig,ax=plt.subplots(figsize=(6.3,3.2));rr=weak[1:];ax.plot([r['cost_factor'] for r in rr],[r['relative_metric_error'] for r in rr],marker='o')
ax.set(xlabel=r'Shortcut waiting-time factor $a$ ($p_{\rm long}=0.5/a$)',ylabel='Relative error to intact SI metric',xscale='log',xticks=[1,4,20,100]);ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter());finish('weak_shortcut_recovery',fig)
# Figure 4: coarsening changes prices relative to spans, not a proven fixed point.
fig,ax=plt.subplots(figsize=(6.3,3.2))
for price in [1,3,10]:
    rr=[r for r in rg if r['price']==price and r['rule']=='min'];ax.plot([r['level'] for r in rr],[r['phi'] for r in rr],marker='o',label=f'Initial long cost {price}')
ax.set(xlabel='Coordinate-block coarsening level',ylabel=r'Underpriced fraction among longs $\Phi$',xticks=[0,1,2,3],ylim=(-.03,1.08));ax.legend(frameon=False);finish('rg_pricing',fig)
# Figure 5: boundary/radial residuals; every seed reported.
fig,ax=plt.subplots(figsize=(6.3,3.3))
for method,label,data in [('hop','Hop',delrows),('shared_field_SI','SI',delsi),('Euclidean_edge_oracle','Euclidean-edge oracle',delrows)]:
    A=np.array([r['radial_quartile_mean_residuals'] for r in data if r['method']==method]);ax.errorbar(np.arange(1,5),A.mean(0),yerr=A.std(0,ddof=1),marker='o',capsize=3,label=label)
ax.axhline(0,linestyle=':',linewidth=.8);ax.set(xlabel='Quartile of squared distance from box center',ylabel='Mean normalized pair-distance residual',xticks=[1,2,3,4]);ax.legend(frameon=False);finish('delaunay_radial_residual',fig)
# Figure 6: radial mode of the unweighted graph.
a=np.load(R/'delaunay_matrices.npz');v=a['hop_v3'];r=a['rho2'];v=v*np.sign(np.corrcoef(v,r)[0,1])
fig,ax=plt.subplots(figsize=(6.3,3.35));ax.scatter(r,v,s=8,alpha=.6);ax.set(xlabel=r'Squared Euclidean radius $r^2$ (evaluation only)',ylabel='Third normalized MDS eigenvector');finish('delaunay_bowl',fig)
# New values for the manuscript, and a machine-readable aggregation.
nums={'SImeanGof':f(means['gof2'],4),'SIsdGof':f(sds['gof2'],4),'SImeanGap':f(means['gap23'],2),'SImeanNeg':f(means['negative_mass'],4),'SImeanCoord':f(means['coordinate_nrmse'],4)}
(T/'numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in nums.items())+'\n')
(R/'aggregates.json').write_text(json.dumps({'si_mean':means,'si_sample_sd':sds,'delaunay_means':summaries},indent=2)+'\n')
print('Saved 6 figures, 8 tables, and machine-readable aggregates.')
