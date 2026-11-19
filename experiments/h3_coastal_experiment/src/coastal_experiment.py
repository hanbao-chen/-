"""Reproducible stage-one H3 coastal density experiment; synthetic coastline."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[1]
LOCAL_DEPS = ROOT.parents[1] / 'work' / 'deps'
if LOCAL_DEPS.exists():
    sys.path.insert(0, str(LOCAL_DEPS))
import argparse, csv, gzip, hashlib, json, math, platform, time
from datetime import datetime, timezone, timedelta
import numpy as np
import h3, shapely
from shapely.geometry import Polygon, mapping
from scipy.stats import poisson

R = 6371.007180918475
LAT0, LNG0 = 30.5, 122.0
METHODS = ['full_area', 'sea_area', 'shrink', 'exclude_small']
TAUS = [0, 0.25, 0.5, 1, 2]

def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')

def project(latlng):
    pts = np.radians(np.asarray(latlng, float))
    lat, lon = pts[:,0], pts[:,1]
    a, b = np.radians([LAT0, LNG0])
    dl = lon-b
    k = np.sqrt(2/(1+np.sin(a)*np.sin(lat)+np.cos(a)*np.cos(lat)*np.cos(dl)))
    x = R*k*np.cos(lat)*np.sin(dl)
    y = R*k*(np.cos(a)*np.sin(lat)-np.sin(a)*np.cos(lat)*np.cos(dl))
    return np.column_stack([x,y])

def dense_boundary(cell):
    boundary=np.radians(np.asarray(h3.cell_to_boundary(cell),float))
    xyz=np.column_stack([np.cos(boundary[:,0])*np.cos(boundary[:,1]),np.cos(boundary[:,0])*np.sin(boundary[:,1]),np.sin(boundary[:,0])])
    out=[]
    for v,w in zip(xyz,np.roll(xyz,-1,axis=0)):
        omega=np.arccos(np.clip(v@w,-1,1))
        steps=max(1,int(np.ceil(R*omega/0.08)))
        t=np.arange(steps)/steps
        p=(np.sin((1-t)*omega)[:,None]*v+np.sin(t*omega)[:,None]*w)/np.sin(omega)
        ll=np.degrees(np.column_stack([np.arcsin(p[:,2]),np.arctan2(p[:,1],p[:,0])]))
        out.extend(project(ll))
    return np.asarray(out)

def geometry():
    cells=sorted(h3.grid_disk(h3.latlng_to_cell(LAT0,LNG0,7),10))
    polys=[Polygon(dense_boundary(c)) for c in cells]
    xy=np.concatenate([np.asarray(p.exterior.coords) for p in polys])
    lo,hi=xy.min(0)-10,xy.max(0)+10
    ys=np.linspace(lo[1],hi[1],1600)
    xs=0.5+2.6*np.sin(ys/3)+0.12*ys
    sea=Polygon(list(zip(xs,ys))+[(hi[0],hi[1]),(hi[0],lo[1])])
    records=[]; features=[]
    for cell,p in zip(cells,polys):
        q=p.intersection(sea)
        if q.area <= 1e-9:
            continue
        h=float(p.area); a=float(q.area)
        center=q.centroid
        records.append({'cell':cell,'full_area_km2':h,'sea_area_km2':a,'sea_fraction':a/h,'x_km':center.x,'y_km':center.y,'h3_area_km2':h3.cell_area(cell,'km^2')})
        features.append({'type':'Feature','properties':{'cell':cell},'geometry':mapping(q),'full_geometry':mapping(p)})
    A=np.array([r['sea_area_km2'] for r in records]); H=np.array([r['full_area_km2'] for r in records])
    F=A/H
    candidates=np.where((F>=.02)&(F<.10))[0]
    if not len(candidates):
        raise RuntimeError('No predefined small-area geometry candidate. Do not change after inspecting results.')
    small=int(min(candidates,key=lambda i: records[i]['x_km']**2+records[i]['y_km']**2))
    full=int(min(np.where(F>1-1e-8)[0],key=lambda i: (records[i]['x_km']-records[small]['x_km'])**2+(records[i]['y_km']-records[small]['y_km'])**2))
    for i,r in enumerate(records):
        r['design_role']='small_hotspot' if i==small else ('full_hotspot' if i==full else 'background')
    rel=np.abs(H-np.array([r['h3_area_km2'] for r in records]))/H
    assert rel.max()<1e-4, rel.max()
    assert np.all((A>0)&(A<=H*(1+1e-12)))
    return records,features,sea,small,full,float(rel.max())

def estimates(counts,A,H,m,tau=1):
    raw=counts/A
    return {'full_area':counts/H,'sea_area':raw,'shrink':(counts+tau*m[:,None])/(A+tau),'exclude_small':np.where(A/H>=.1,raw,np.nan)}

def bootstrap_paired(d,seed=4501):
    rng=np.random.default_rng(seed)
    means=[]
    for _ in range(20):
        ix=rng.integers(0,len(d),size=(250,len(d)))
        means.extend(d[ix].mean(axis=1))
    return {'mean':float(np.mean(d)),'ci95':np.quantile(means,[.025,.975]).tolist(),'unit':'independent simulated snapshot repeat'}

def write_csv(path,rows):
    if not rows: return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def run(out,n):
    start=time.time();out.mkdir(parents=True,exist_ok=True)
    records,features,sea,small,full,areaerr=geometry()
    write_csv(out/'cells.csv',records)
    dump(out/'geometry.json',{'coordinate_system':'spherical Lambert azimuthal equal-area, km; synthetic coast','features':features,'sea_polygon':mapping(sea)})
    A=np.array([r['sea_area_km2'] for r in records]);H=np.array([r['full_area_km2'] for r in records]);F=A/H
    config={'run_type':'stage_one_fixed_design','n_repeats_per_scenario':n,'first_seed':2026111900,'h3_resolution':7,'disk_k':10,'projection_center_lat_lng':[LAT0,LNG0],'baseline_density':.5,'hotspot_density':4.,'hotspot_threshold_strictly_greater_than':2.,'historical_snapshots':20,'main_tau_km2':1.,'sensitivity_tau_km2':TAUS,'small_hotspot_index':small,'full_hotspot_index':full,'small_hotspot_cell':records[small],'full_hotspot_cell':records[full]}
    dump(out/'config.json',config)
    summary={'geometry':{'sea_cells':len(A),'full_disk_cells':331,'small_fraction_cells':int(np.sum(F<.1)),'total_sea_area_km2':float(A.sum()),'max_projection_area_relative_error':areaerr},'scenarios':{}}
    trialrows=[];sensitivity=[];comparisons={};checks={}
    for scidx,scenario in enumerate(['uniform','hotspots']):
        truth=np.full(len(A),.5)
        if scenario=='hotspots':truth[[small,full]]=4.
        bg=truth<=2;hot=~bg
        counts=np.empty((n,len(A)),dtype=np.int32)
        hist=np.empty_like(counts)
        for j in range(n):
            rs=np.random.SeedSequence([2026111900+j,scidx])
            rtest,rhist=[np.random.default_rng(s) for s in rs.spawn(2)]
            counts[j]=rtest.poisson(truth*A)
            hist[j]=rhist.poisson(20*truth*A)
        m=hist.sum(axis=1)/(20*A.sum())
        np.savez_compressed(out/f'{scenario}_counts.npz',counts=counts,history_counts=hist,prior_mean=m,truth=truth,sea_area=A,full_area=H,seeds=np.arange(2026111900,2026111900+n))
        methods=estimates(counts,A,H,m)
        scenario_summary={};per={}
        for method,est in methods.items():
            pred=est>2
            fp=pred[:,bg].sum(axis=1)
            metrics={'false_positive_cells':fp,'background_false_positive_rate':fp/bg.sum(),'small_hotspot_detected':pred[:,small].astype(int) if scenario=='hotspots' else np.full(n,np.nan),'full_hotspot_detected':pred[:,full].astype(int) if scenario=='hotspots' else np.full(n,np.nan)}
            valid=np.isfinite(est[0])
            sq=(est[:,valid]-truth[valid])**2
            metrics['mse_cell_equal_weight_published']=sq.mean(axis=1)
            metrics['mse_area_weight_published']=(sq*A[valid]).sum(axis=1)/A[valid].sum()
            cov=float(valid.mean())
            means={k:float(np.mean(v)) if np.isfinite(v).all() else None for k,v in metrics.items()}
            means['published_cell_fraction']=cov
            means['published_sea_area_fraction']=float(A[valid].sum()/A.sum())
            # Independent Poisson-tail evaluation, conditional on each history-derived prior.
            if method=='shrink':cut=2*(A[None,:]+1)-m[:,None]
            elif method=='full_area':cut=np.broadcast_to(2*H,counts.shape)
            else:cut=np.broadcast_to(2*A,counts.shape)
            prob=poisson.sf(np.floor(cut),truth[None,:]*A[None,:])
            if method=='exclude_small':prob[:,F<.1]=0
            exact={'expected_false_positive_cells':float(prob[:,bg].sum(axis=1).mean()),'small_hotspot_detection_probability':float(prob[:,small].mean()) if scenario=='hotspots' else None,'full_hotspot_detection_probability':float(prob[:,full].mean()) if scenario=='hotspots' else None}
            scenario_summary[method]={'empirical':means,'poisson_tail_check':exact}
            per[method]=metrics
            # Allow five estimated standard errors plus a small finite-sample margin.
            for key,val in [('false_positive_cells',exact['expected_false_positive_cells']),('small_hotspot_detected',exact['small_hotspot_detection_probability']),('full_hotspot_detected',exact['full_hotspot_detection_probability'])]:
                if val is None:continue
                obs=metrics[key]
                delta=float(abs(obs.mean()-val));bound=float(5*np.std(obs,ddof=1)/np.sqrt(n)+.005)
                assert delta<bound,(scenario,method,key,delta,bound)
                checks[f'{scenario}/{method}/{key}']={'absolute_difference':delta,'allowed_monte_carlo_difference':bound,'passed':True}
            for j in range(n):
                row={'scenario':scenario,'seed':2026111900+j,'method':method,'prior_mean':float(m[j])}
                row.update({k:float(v[j]) if np.isfinite(v[j]) else '' for k,v in metrics.items()})
                trialrows.append(row)
        summary['scenarios'][scenario]=scenario_summary
        comparisons[f'{scenario}_sea_minus_shrink_false_positive_cells']=bootstrap_paired(per['sea_area']['false_positive_cells']-per['shrink']['false_positive_cells'])
        if scenario=='hotspots':
            for key in ['small_hotspot_detected','full_hotspot_detected']:
                comparisons[f'{scenario}_sea_minus_shrink_{key}']=bootstrap_paired(per['sea_area'][key]-per['shrink'][key])
        for tau in TAUS:
            est=(counts+tau*m[:,None])/(A+tau);pred=est>2
            sensitivity.append({'scenario':scenario,'tau_km2':tau,'false_positive_cells':float(pred[:,bg].sum(axis=1).mean()),'small_hotspot_recall':float(pred[:,small].mean()) if scenario=='hotspots' else '', 'full_hotspot_recall':float(pred[:,full].mean()) if scenario=='hotspots' else ''})
        # Analytical unbiasedness and variance; empirical means are a stochastic check.
        direct=counts/A
        bias=direct.mean(axis=0)-truth
        checks[scenario+'/raw_density_expectation']={'max_abs_standardized_mean_error':float(np.max(np.abs(bias)/np.sqrt(truth/A/n))),'theory':'E[N/A]=lambda; Var[N/A]=lambda/A','passed':bool(np.max(np.abs(bias)/np.sqrt(truth/A/n))<6)}
        assert checks[scenario+'/raw_density_expectation']['passed']
        assert np.array_equal(estimates(counts,A,H,m,0)['shrink'],direct)
        print(scenario,json.dumps(scenario_summary,ensure_ascii=False),flush=True)
    summary['paired_bootstrap']=comparisons
    write_csv(out/'trials.csv',trialrows);write_csv(out/'sensitivity.csv',sensitivity)
    dump(out/'summary.json',summary);dump(out/'checks.json',checks)
    files=[p for p in out.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += list((ROOT/'src').glob('*.py'))+[ROOT/'EXPERIMENT_PLAN.md',ROOT/'requirements.txt']
    manifest={'created_at_hong_kong':datetime.now(timezone(timedelta(hours=8))).isoformat(),'wall_seconds':time.time()-start,'environment':{'python':platform.python_version(),'platform':platform.platform(),'numpy':np.__version__,'h3':h3.versions(),'shapely':shapely.__version__},'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    dump(out/'manifest.json',manifest)
    print('GEOMETRY',json.dumps(summary['geometry']));print('HOTSPOT_CELLS',json.dumps({'small':records[small],'full':records[full]}));print('PAIRED',json.dumps(comparisons));print('WALL_SECONDS',time.time()-start)

if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'results'/'stage1');parser.add_argument('--repeats',type=int,default=2000)
    args=parser.parse_args();run(args.output,args.repeats)
