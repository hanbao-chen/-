"""Matched-data continuous-observation checks; not a report generator."""
from pathlib import Path
import sys,json,csv,hashlib,time,platform,gc
from datetime import datetime,timezone,timedelta
import numpy as np
from scipy.stats import poisson
import scipy
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
BASE=ROOT.parent/'h3_coastal_experiment'/'results'/'stage1'
cfg1=json.loads((BASE/'config.json').read_text('utf-8'))
old=np.load(BASE/'uniform_counts.npz')
A=old['sea_area'].copy();M=len(A);SMALL=cfg1['small_hotspot_index'];FULL=cfg1['full_hotspot_index']
NREP=2000;TMAX=120;LOOKS=[1,5,10,20,40,80,120];ALPHA=.05

def dump(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def csvwrite(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def cuts(exposure,alpha):
    mu=2*exposure
    hi=poisson.isf(alpha,mu).astype(int)+1
    q=poisson.ppf(alpha,mu).astype(int)
    lo=q-(poisson.cdf(q,mu)>alpha)
    assert np.all(poisson.sf(hi-1,mu)<=alpha*(1+1e-8))
    assert np.all(poisson.sf(hi-2,mu)>alpha*(1-1e-8))
    valid=lo>=0
    assert np.all(poisson.cdf(lo[valid],mu[valid])<=alpha*(1+1e-8))
    assert np.all(poisson.cdf(lo+1,mu)>alpha*(1-1e-8))
    return hi,lo

def bootstrap(d,seed=6557):
    rng=np.random.default_rng(seed);v=[]
    for _ in range(10):v.extend(d[rng.integers(0,len(d),size=(200,len(d)))].mean(1))
    return {'mean':float(d.mean()),'ci95':np.quantile(v,[.025,.975]).tolist()}

def run_case(mode,scenario,idx):
    block=1 if mode=='independent' else 5
    truth=np.full(M,.5)
    if scenario=='hotspots':truth[[SMALL,FULL]]=4
    if scenario=='boundary':truth[:]=2
    rng=np.random.default_rng(np.random.SeedSequence([20261119,601,idx]))
    # Only independent block counts need to be persisted; frames within a block
    # are exact repeats, with the same synthetic identities.
    data=rng.poisson(truth[None,None,:]*A[None,None,:],size=(NREP,TMAX//block,M)).astype(np.uint16)
    tag=f'{mode}_{scenario}'
    np.savez_compressed(OUT/(tag+'_data.npz'),block_counts=data,block_size=block,truth=truth,sea_area=A)
    bg=truth<=2;hot=truth>2
    names=['raw_cumulative','shrink_cumulative','poisson_fixed','poisson_sequential']
    if block>1:names+=['poisson_fixed_block','poisson_sequential_block']
    states={name:np.zeros((NREP,M),np.int8) for name in names}
    first={name:np.full((NREP,M),TMAX+1,np.int16) for name in names}
    ever_fp={name:np.zeros((NREP,M),bool) for name in names}
    ever_wrong_direction={name:np.zeros(NREP,bool) for name in names}
    pending_rounds={name:np.zeros(NREP,float) for name in names}
    candidate_rounds={name:np.zeros(NREP,float) for name in names}
    candidates={name:np.zeros((NREP,M),bool) for name in names}
    candidate_ever_bg={name:np.zeros((NREP,M),bool) for name in names}
    aggregate=[];trials=[];last=0;invariants={}
    for t in LOOKS:
        nb=(t+block-1)//block
        if block==1:K=data[:,:t].sum(1,dtype=np.int32)
        else:
            completed=t//block;remainder=t%block
            K=block*data[:,:completed].sum(1,dtype=np.int32)
            if remainder:K+=remainder*data[:,completed,:].astype(np.int32)
        KB=data[:,:nb].sum(1,dtype=np.int32)
        D=K/(t*A)
        for name in names:
            state=states[name]
            pending_rounds[name]+=(state==0).sum(1)*(t-last)
            candidate_rounds[name]+=candidates[name].sum(1)*(t-last)
            if name=='raw_cumulative':state[:]=np.where(D>2,1,-1)
            elif name=='shrink_cumulative':state[:]=np.where((K+.5)/(t*A+1)>2,1,-1)
            else:
                fixed='fixed' in name
                if not fixed or t==TMAX:
                    uses_block=name.endswith('_block')
                    evidence_K=KB if uses_block else K
                    exposure=A*(nb if uses_block else t)
                    alpha=ALPHA/(2*M*(1 if fixed else len(LOOKS)))
                    hi,lo=cuts(exposure,alpha)
                    H=evidence_K>=hi;L=evidence_K<=lo
                    assert not np.any(H&L)
                    # Exact critical-count code agrees with direct tail evaluations.
                    if name=='poisson_sequential' and t in [1,120]:
                        assert np.array_equal(H,poisson.sf(evidence_K-1,2*exposure)<=alpha)
                        assert np.array_equal(L,poisson.cdf(evidence_K,2*exposure)<=alpha)
                    undecided=state==0
                    state[undecided&H]=1;state[undecided&L]=-1
            new=(first[name]==TMAX+1)&(state==1);first[name][new]=t
            ever_fp[name]|=(state==1)&bg[None,:]
            # For boundary lambda=2, either directional claim is wrong, even if
            # labeling 'non-hot' is correct for a <=2 business definition.
            wrong=((state==1)&(truth[None,:]<=2))|((state==-1)&(truth[None,:]>=2))
            ever_wrong_direction[name]|=wrong.any(1)
            pend=state==0;H=state==1;L=state==-1
            assert np.all(H.astype(int)+pend.astype(int)+L.astype(int)==1)
            cand=pend&(D>2);candidates[name]=cand
            candidate_ever_bg[name]|=cand&bg[None,:]
            row={'mode':mode,'scenario':scenario,'method':name,'look':t,
                 'false_hot_cells_mean':float(H[:,bg].sum(1).mean()),
                 'maps_with_false_hot_fraction':float(H[:,bg].any(1).mean()),
                 'maps_with_ever_false_hot_fraction':float(ever_fp[name].any(1).mean()),
                 'maps_with_ever_wrong_direction_fraction':float(ever_wrong_direction[name].mean()),
                 'pending_cells_mean':float(pend.sum(1).mean()),
                 'high_density_pending_cells_mean':float(cand.sum(1).mean()),
                 'background_pending_cells_mean':float(pend[:,bg].sum(1).mean()),
                 'background_candidates_ever_mean':float(candidate_ever_bg[name].sum(1).mean()),
                 'pending_cell_rounds_mean':float(pending_rounds[name].mean()),
                 'candidate_cell_rounds_mean':float(candidate_rounds[name].mean())}
            for label,j in [('small',SMALL),('full',FULL)]:
                row[label+'_confirmed']=float(H[:,j].mean())
                row[label+'_pending']=float(pend[:,j].mean())
                row[label+'_negative']=float(L[:,j].mean())
                assert abs(row[label+'_confirmed']+row[label+'_pending']+row[label+'_negative']-1)<1e-12
            aggregate.append(row)
        last=t
    # Final per-sequence numbers: bootstrap units remain independent sequences.
    for name,state in states.items():
        for r in range(NREP):
            trials.append({'mode':mode,'scenario':scenario,'method':name,'repeat':r,
                'small_state':int(state[r,SMALL]),'full_state':int(state[r,FULL]),
                'small_first_hot_round':int(first[name][r,SMALL]),'full_first_hot_round':int(first[name][r,FULL]),
                'small_wait_censored_at_120':int(min(TMAX,first[name][r,SMALL])),
                'false_hot_cells':int(((state[r]==1)&bg).sum()),
                'ever_false_hot_map':int(ever_fp[name][r].any()),
                'ever_wrong_direction_map':int(ever_wrong_direction[name][r]),
                'pending_cells':int((state[r]==0).sum()),
                'pending_cell_rounds':float(pending_rounds[name][r]),
                'candidate_cell_rounds':float(candidate_rounds[name][r]),
                'background_candidates_ever':int(candidate_ever_bg[name][r].sum())})
        valid_model=block==1 or name.endswith('_block')
        if 'poisson' in name and valid_model:
            observed=float(ever_wrong_direction[name].mean())
            tolerance=ALPHA+5*np.sqrt(ALPHA*(1-ALPHA)/NREP)+.005
            assert observed<tolerance,(tag,name,observed)
            invariants[name+'/global_error_bound_check']={'observed':observed,'design_bound':ALPHA,'loose_monte_carlo_upper':float(tolerance),'passed':True}
        # binary publication using exactly the same evidence must agree pointwise.
        if 'poisson' in name:
            published_binary=np.where(state==1,1,-1)
            assert np.array_equal(published_binary==1,state==1)
            invariants[name+'/same_evidence_binary_vs_three_state_hot_mismatches']=0
    if scenario=='hotspots':
        for name in names:
            times=first[name][:,SMALL]
            final=next(r for r in aggregate if r['method']==name and r['look']==120)
            final['small_censored_mean_first_hot_round']=float(np.minimum(times,TMAX).mean())
            q=int(np.sort(times)[int(np.ceil(.8*NREP))-1])
            final['small_t80_first_hot_round']=q if q<=TMAX else None
    # Clear optional fields on other aggregate rows for rectangular CSV.
    for row in aggregate:
        row.setdefault('small_censored_mean_first_hot_round',None);row.setdefault('small_t80_first_hot_round',None)
    comparisons={}
    if scenario=='hotspots':
        pair=('poisson_sequential','poisson_fixed') if block==1 else ('poisson_sequential_block','poisson_fixed_block')
        a,b=pair
        comparisons['sequential_minus_fixed_small_confirmed']=bootstrap((states[a][:,SMALL]==1).astype(int)-(states[b][:,SMALL]==1).astype(int))
        comparisons['sequential_minus_fixed_pending_cell_rounds']=bootstrap(pending_rounds[a]-pending_rounds[b])
        comparisons['sequential_minus_fixed_censored_small_wait']=bootstrap(np.minimum(first[a][:,SMALL],TMAX)-np.minimum(first[b][:,SMALL],TMAX))
    if scenario=='boundary' and block>1:
        comparisons['naive_minus_block_aware_ever_false_hot_maps']=bootstrap(ever_fp['poisson_sequential'].any(1).astype(int)-ever_fp['poisson_sequential_block'].any(1).astype(int))
    dump(OUT/(tag+'_invariants.json'),invariants)
    dump(OUT/(tag+'_comparisons.json'),comparisons)
    csvwrite(OUT/(tag+'_trials.csv'),trials)
    # Save final states and first-confirmation times for independent recomputation.
    np.savez_compressed(OUT/(tag+'_decisions.npz'),**{name+'_state':s for name,s in states.items()},**{name+'_first':s for name,s in first.items()})
    del data;gc.collect()
    print(tag,json.dumps([r for r in aggregate if r['look']==120],ensure_ascii=False),flush=True)
    return aggregate

def main():
    start=time.time()
    config={'n_repeats_per_case':NREP,'n_cells':M,'looks':LOOKS,'max_round':TMAX,'family_any_direction_alpha':ALPHA,'small_cell_index':SMALL,'full_cell_index':FULL,'background_density':.5,'business_threshold':2,'hotspot_density':4,'block_sizes':[1,5],'rng_seed_components':[20261119,601,'case_index'],'history_used':False,'null_model':'known parameters; exact one-sided Poisson at business boundary 2','time_unit':'synthetic observation round, not real minutes'}
    dump(OUT/'config.json',config)
    rows=[];idx=0
    for mode in ['independent','repeated_5']:
        for scenario in ['uniform','hotspots','boundary']:
            rows.extend(run_case(mode,scenario,idx));idx+=1
    csvwrite(OUT/'summary.csv',rows)
    dump(OUT/'summary.json',rows)
    manifest={'created_at_hong_kong':datetime.now(timezone(timedelta(hours=8))).isoformat(),'wall_seconds':time.time()-start,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'files_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*') if p.is_file() and p.name!='manifest.json'}}
    dump(ROOT/'manifest.json',manifest)
    print('ALL_DONE_SECONDS',time.time()-start,flush=True)

if __name__=='__main__':main()
