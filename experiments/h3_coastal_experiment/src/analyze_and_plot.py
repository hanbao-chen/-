"""Read saved counts, independently recompute metrics, and make scientific figures."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
deps=ROOT.parents[1]/'work'/'deps'
if deps.exists():sys.path.insert(0,str(deps))
import csv, hashlib, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from shapely.geometry import shape
from scipy.stats import poisson

RESULTS=ROOT/'results'/'stage1'
FIGURES=ROOT/'figures'
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.unicode_minus':False,'savefig.facecolor':'white','figure.facecolor':'white'})
summary=json.loads((RESULTS/'summary.json').read_text('utf-8'))
cfg=json.loads((RESULTS/'config.json').read_text('utf-8'))
geo=json.loads((RESULTS/'geometry.json').read_text('utf-8'))
small,full=cfg['small_hotspot_index'],cfg['full_hotspot_index']
checks=[]
for rel,expected in json.loads((RESULTS/'manifest.json').read_text('utf-8'))['sha256'].items():
    actual=hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()
    assert actual==expected, rel
checks.append('all original run manifest SHA-256 entries matched')
names={'full_area':'除以完整格面积','sea_area':'除以实际海域面积','shrink':'海域面积＋保守修正','exclude_small':'不发布海域占比<10%的格子'}
detail={}
for scenario in ['uniform','hotspots']:
    data=np.load(RESULTS/f'{scenario}_counts.npz')
    N,A,H,m,truth=[data[k] for k in ['counts','sea_area','full_area','prior_mean','truth']]
    est={'full_area':N/H,'sea_area':N/A,'shrink':(N+m[:,None])/(A+1),'exclude_small':np.where(A/H>=.1,N/A,np.nan)}
    bg=truth<=2
    detail[scenario]={}
    for method,D in est.items():
        P=D>2
        fp=P[:,bg].sum(1)
        expected=summary['scenarios'][scenario][method]['empirical']
        assert np.isclose(fp.mean(),expected['false_positive_cells'])
        if scenario=='hotspots':
            assert np.isclose(P[:,small].mean(),expected['small_hotspot_detected'])
            assert np.isclose(P[:,full].mean(),expected['full_hotspot_detected'])
        detail[scenario][method]={'small_cell_positive_fraction':float(P[:,small].mean()),'full_cell_positive_fraction':float(P[:,full].mean()),'snapshots_with_any_false_positive_fraction':float((fp>0).mean()),'false_positive_cells_per_snapshot':float(fp.mean()),'total_false_positive_decisions':int(fp.sum()),'small_cell_positive_count':int(P[:,small].sum()),'full_cell_positive_count':int(P[:,full].sum())}
    detail[scenario]['input']={'observed_current_vessel_total':int(N.sum()),'historic_vessel_total':int(data['history_counts'].sum()),'n_repeats':len(N),'small_cell_expected_current_count':float(truth[small]*A[small]),'full_cell_expected_current_count':float(truth[full]*A[full]),'small_cell_count_histogram':np.bincount(N[:,small]).tolist(),'prior_mean_min':float(m.min()),'prior_mean_max':float(m.max())}
checks.append('all false-positive means and all hotspot recalls independently recomputed from saved matrices')

# 1. Outcome bars. Two main methods, plus full-area context in a separate note.
fig,axes=plt.subplots(1,3,figsize=(12.8,4.8))
methods=['sea_area','shrink'];colors=['#286c8f','#ca7833']
panels=[('uniform','false_positive_cells','均匀海域：误报','每张快照平均误报格数',1),('hotspots','small_hotspot_detected','小海域真实热点：检出','检出率（%）',100),('hotspots','full_hotspot_detected','完整海域真实热点：检出','检出率（%）',100)]
for ax,(sc,key,title,ylabel,scale) in zip(axes,panels):
    vals=[summary['scenarios'][sc][m]['empirical'][key]*scale for m in methods]
    bars=ax.bar([0,1],vals,color=colors,width=.57)
    ax.set_xticks([0,1],['按海域面积','保守修正']);ax.set_title(title,pad=15);ax.set_ylabel(ylabel)
    ax.set_ylim(0,1.15 if scale==1 else 112)
    ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
    for bar,v in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,v+(0.025 if scale==1 else 2),f'{v:.3f}' if scale==1 else f'{v:.2f}%',ha='center',va='bottom',fontweight='bold')
fig.suptitle('保守修正减少误报，也可能压低真实小热点',fontsize=16,y=.99)
fig.text(.5,.015,'真实 H3 r7 网格＋人工岸线；每场景 2,000 个独立模拟快照；热点阈值 >2 艘/km²；修正强度 τ=1 km²',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.065,1,.92]);fig.savefig(FIGURES/'01_main_results.png',dpi=190);plt.close(fig)

# 2. Sensitivity: a trade-off view, not optimized on these test results.
with (RESULTS/'sensitivity.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
u=[r for r in rows if r['scenario']=='uniform'];h=[r for r in rows if r['scenario']=='hotspots']
fig,ax=plt.subplots(figsize=(8.5,5.3))
x=[float(r['false_positive_cells']) for r in u]
y=[float(r['small_hotspot_recall'])*100 for r in h]
ax.plot(x,y,'o-',color='#286c8f',lw=2,ms=8)
for a,b,r in zip(x,y,u):ax.annotate(f"τ={float(r['tau_km2']):g}",(a,b),xytext=(8,5),textcoords='offset points')
ax.set_xlim(-.025,1.1);ax.set_ylim(-3,85);ax.set_xlabel('均匀海域：每张快照平均误报格数（越少越好）');ax.set_ylabel('小海域真实热点检出率（%，越高越好）');ax.set_title('修正强度不同，误报与漏报的取舍也不同',pad=14);ax.grid(alpha=.18)
fig.text(.5,.015,'τ=0 为直接按海域面积计算。各点来自预先列出的强度，不据此宣称普遍最优参数。',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(FIGURES/'02_tradeoff.png',dpi=190);plt.close(fig)

# 3. Mechanism with empirical count frequencies and exact model PMF.
d=np.load(RESULTS/'hotspots_counts.npz');N=d['counts'];A=d['sea_area'];m=d['prior_mean']
fig,ax=plt.subplots(figsize=(8.8,5.1))
k=np.arange(8);observed=np.array([(N[:,small]==i).mean() for i in k])*100
ax.bar(k,observed,color=['#aeb8bf' if i==0 else ('#75a5be' if i<3 else '#ca7833') for i in k],width=.65,label='2,000 次模拟的观测频率')
ax.plot(k,poisson.pmf(k,4*A[small])*100,'o--',color='#303f49',ms=5,label='泊松模型理论概率')
for i,v in zip(k,observed):
    if v>.1:ax.text(i,v+.6,f'{v:.1f}%',ha='center',fontsize=9)
ax.axvline(.5,color='#286c8f',ls=':',lw=2,label='直接面积法：至少 1 艘即判热点')
ax.axvline(2.5,color='#ca7833',ls=':',lw=2,label='本次保守修正：至少 3 艘才判热点')
ax.set_xticks(k);ax.set_xlabel('这个小海域格内，单张快照看到的船数');ax.set_ylabel('快照比例（%）');ax.set_ylim(0,53);ax.set_title(f'真实热点也常常只看到 0–2 艘船\n海域面积 {A[small]:.3f} km²，真密度 4 艘/km²，平均仅 {4*A[small]:.2f} 艘/快照',pad=12);ax.legend(frameon=False,fontsize=9,loc='upper right');ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(FIGURES/'03_why_hotspot_is_missed.png',dpi=190);plt.close(fig)

# 4. True field and first predesignated snapshot. All panels share one density scale.
polys=[shape(f['geometry']) for f in geo['features']]
verts=[];index=[]
for i,p in enumerate(polys):
    parts=[p] if p.geom_type=='Polygon' else list(p.geoms)
    for part in parts:
        if part.geom_type=='Polygon':verts.append(np.asarray(part.exterior.coords));index.append(i)
idx=np.array(index)
raw=N[0]/A;shrunk=(N[0]+m[0])/(A+1);truth=d['truth']
fig,axes=plt.subplots(1,3,figsize=(12.8,6.2),sharex=True,sharey=True)
for ax,values,title in zip(axes,[truth,raw,shrunk],['模型中的真实密度','首个种子的直接面积估计','同一快照的保守修正']):
    pc=PolyCollection(verts,array=values[idx],cmap='YlOrRd',norm=Normalize(0,6),edgecolors='#858585',linewidths=.25)
    ax.add_collection(pc)
    for j,marker in [(small,'*'),(full,'o')]:
        p=polys[j].centroid
        ax.scatter([p.x],[p.y],marker=marker,s=135 if j==small else 80,facecolor='none',edgecolor='#206694',linewidth=1.5,zorder=4)
    ax.set_xlim(-7,10);ax.set_ylim(-11,6);ax.set_aspect('equal');ax.set_xlabel('相对东向距离（km）');ax.set_title(title,pad=10);ax.set_facecolor('#e9edef')
axes[0].set_ylabel('相对北向距离（km）')
fig.subplots_adjust(left=.065,right=.89,bottom=.29,top=.84,wspace=.12)
cax=fig.add_axes([.91,.31,.016,.48]);fig.colorbar(pc,cax=cax,label='艘/km²（上限 6）',extend='max')
fig.suptitle('人工岸线下的两个同强度真实热点（局部放大）',fontsize=16,y=.97)
fig.legend(handles=[Line2D([0],[0],marker='*',color='none',markeredgecolor='#206694',markerfacecolor='none',markersize=12,label='小海域真实热点'),Line2D([0],[0],marker='o',color='none',markeredgecolor='#206694',markerfacecolor='none',markersize=9,label='完整海域真实热点')],loc='lower center',bbox_to_anchor=(.5,.105),ncol=2,frameon=False)
fig.text(.5,.025,'灰色为陆地或图外；海岸线为人工构造。例图固定用种子 2026111900，结论依据全部重复。',ha='center',fontsize=9)
fig.savefig(FIGURES/'04_synthetic_map.png',dpi=190);plt.close(fig)

detail['verification']=checks
(RESULTS/'independent_verification.json').write_text(json.dumps(detail,ensure_ascii=False,indent=2),encoding='utf-8')
manifest={'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [RESULTS/'summary.json',RESULTS/'sensitivity.csv',RESULTS/'uniform_counts.npz',RESULTS/'hotspots_counts.npz',Path(__file__)]},'figures':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in FIGURES.glob('*.png')}}
(FIGURES/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(detail,ensure_ascii=False,indent=2))
