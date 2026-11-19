"""Independently replay selected saved sequences and display the matched results."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
from scipy.stats import poisson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent;R=ROOT/'results'
config=json.loads((R/'config.json').read_text('utf-8'))
rows=json.loads((R/'summary.json').read_text('utf-8'))
checks=[]
manifest=json.loads((ROOT/'manifest.json').read_text('utf-8'))
for rel,digest in manifest['files_sha256'].items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==digest,rel
checks.append('all original experiment file SHA-256 fingerprints matched')
selected=[0,1,2,3,4,5,6,7,8,9,1999]
for mode in ['independent','repeated_5']:
    for scenario in ['uniform','hotspots','boundary']:
        tag=f'{mode}_{scenario}'
        data=np.load(R/(tag+'_data.npz'))
        B=int(data['block_size']);A=data['sea_area'];N=data['block_counts'][selected]
        decisions=np.load(R/(tag+'_decisions.npz'))
        methods=[k[:-6] for k in decisions.files if k.endswith('_state')]
        for method in methods:
            S=np.zeros((len(selected),len(A)),np.int8)
            first=np.full(S.shape,121,int)
            for t in config['looks']:
                # Explicit reconstruction of frames, different from the optimized
                # weighted-block sum in the experimental implementation.
                frames=np.repeat(N,B,axis=1)[:,:t]
                K=frames.sum(1,dtype=np.int64);D=K/(A*t)
                if method=='raw_cumulative':S=np.where(D>2,1,-1)
                elif method=='shrink_cumulative':S=np.where((K+.5)/(A*t+1)>2,1,-1)
                elif 'fixed' not in method or t==120:
                    block=method.endswith('_block')
                    if block:
                        independent=frames[:,::B]
                        KK=independent.sum(1,dtype=np.int64);E=A*independent.shape[1]
                    else:KK=K;E=A*t
                    alpha=.05/(2*len(A)*(1 if 'fixed' in method else len(config['looks'])))
                    P=poisson.sf(KK-1,2*E);Q=poisson.cdf(KK,2*E)
                    old=S.copy()
                    S[(old==0)&(P<=alpha)]=1
                    S[(old==0)&(Q<=alpha)]=-1
                first[(first==121)&(S==1)]=t
            assert np.array_equal(S,decisions[method+'_state'][selected]),(tag,method)
            assert np.array_equal(first,decisions[method+'_first'][selected]),(tag,method,'time')
        checks.append(tag+': all methods replayed with direct p-values on 11 predetermined sequences')
for row in rows:
    for k in ['small','full']:
        assert np.isclose(sum(row[k+'_'+s] for s in ['confirmed','pending','negative']),1)
checks.append('all three-state proportions sum to one for both target cells at every look')
(ROOT/'verification.json').write_text(json.dumps({'checks':checks,'passed':True},ensure_ascii=False,indent=2),encoding='utf-8')

plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
fig,ax=plt.subplots(figsize=(9.2,5.8))
groups=[('independent','poisson_sequential','独立快照：校正后的逐步检验','#286c8f','-'),('repeated_5','poisson_sequential_block','每5张重复：按独立船组校正','#cb7833','-'),('repeated_5','poisson_sequential','每5张重复：错误当作独立证据','#7b7b7b','--')]
for mode,method,label,color,ls in groups:
    v=sorted([r for r in rows if r['mode']==mode and r['scenario']=='hotspots' and r['method']==method],key=lambda r:r['look'])
    ax.plot([r['look'] for r in v],[r['small_confirmed']*100 for r in v],marker='o',lw=2,color=color,ls=ls,label=label)
ax.set_xlim(0,125);ax.set_ylim(0,108);ax.set_xticks([1,20,40,80,120]);ax.set_xlabel('实际观察轮次');ax.set_ylabel('小海域真实热点确认率（%）');ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=10,loc='upper left')
ax.set_title('重复照片会让证据显得比实际更充分',fontsize=15,pad=13)
fig.text(.5,.015,'每组2,000条模拟序列；观察时间相同。灰线的错误模型在边界场景中失去误报控制，不能据其高确认率选优。',ha='center',fontsize=8.5)
fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(ROOT/'continuous_confirmation.png',dpi=180);plt.close(fig)

print('VERIFIED',len(checks))
for mode,method in [('independent','poisson_sequential'),('repeated_5','poisson_sequential_block')]:
    print(mode,method)
    for r in rows:
        if r['mode']==mode and r['scenario']=='hotspots' and r['method']==method:
            print({k:r[k] for k in ['look','small_confirmed','small_pending','small_negative','pending_cells_mean','high_density_pending_cells_mean','background_candidates_ever_mean','candidate_cell_rounds_mean']})
print('PAIRED_CIS')
for p in sorted(R.glob('*comparisons.json')):
    print(p.stem,json.loads(p.read_text('utf-8')))
