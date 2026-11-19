"""Post-hoc diagnostic requested by the student; no change to stage-one results."""
from pathlib import Path
import sys,json
import numpy as np
from scipy.stats import poisson
sys.stdout.reconfigure(encoding='utf-8')
root=Path('outputs/h3_coastal_experiment/results/stage1')
save=Path('work/threshold_check');save.mkdir(exist_ok=True)
u=np.load(root/'uniform_counts.npz')
A=u['sea_area'];L=.5;T=2;tau=1
cfg=json.loads((root/'config.json').read_text('utf-8'))
s,f=cfg['small_hotspot_index'],cfg['full_hotspot_index']

# Calibrate a fixed threshold using ONLY the known null model. No hotspot
# observations or false-positive observations from the saved test arrays used.
# Integrate over the distribution of the independent historical total.
hist_mean=20*L*A.sum()
totals=np.arange(max(0,int(poisson.ppf(1e-14,hist_mean))),int(poisson.ppf(1-1e-14,hist_mean))+1)
weights=poisson.pmf(totals,hist_mean)
assert abs(weights.sum()-1)<1e-10
prior=totals/(20*A.sum())
cutoff=T*(A[None,:]+tau)-tau*prior[:,None]
conditional_fp=poisson.sf(np.floor(cutoff),L*A[None,:]).sum(1)
target=float(weights@conditional_fp)
def fixed_fp(c):return float(poisson.sf(np.floor(c*A),L*A).sum())
low,high=2.,4.
while fixed_fp(high)>target:high*=2
for _ in range(100):
    mid=(low+high)/2
    if fixed_fp(mid)>target:low=mid
    else:high=mid
matched=high
result={'comparison_type':'post-hoc diagnostic; all saved stage-one observations retained',
        'fixed_threshold_calibration':'known Poisson uniform null model only; oracle diagnostic, not a fitted deployment method',
        'shrink_null_expected_fp':target,'fixed_threshold':matched,'fixed_null_expected_fp':fixed_fp(matched),
        'fixed_threshold_low_bracket':low,'fixed_null_expected_fp_low_bracket':fixed_fp(low),
        'scenarios':{}}
for scenario in ['uniform','hotspots']:
    z=np.load(root/(scenario+'_counts.npz'));N=z['counts'];m=z['prior_mean'];D=N/A
    shrink=(N+tau*m[:,None])/(A+tau)>T
    effective=T+tau*(T-m[:,None])/A
    adaptive=D>effective
    mismatch=int((shrink!=adaptive).sum())
    assert mismatch==0
    assert np.all(m<T)
    assert np.all(~shrink | (D>T))
    estimates=[('raw_threshold_2',D>2),('raw_threshold_4',D>4),('raw_threshold_8',D>8),('raw_fixed_null_matched',D>matched),('shrink_tau_1_threshold_2',shrink),('raw_exact_adaptive_threshold',adaptive)]
    truth=z['truth'];bg=truth<=T
    rows={}
    for name,pred in estimates:
        rows[name]={'false_positive_cells_per_snapshot':float(pred[:,bg].sum(1).mean()),
            'small_cell_positive_fraction':float(pred[:,s].mean()),
            'full_cell_positive_fraction':float(pred[:,f].mean()),
            'different_decisions_from_shrink':int((pred!=shrink).sum())}
    result['scenarios'][scenario]={'identical_adaptive_threshold_mismatches':mismatch,'decision_count':int(shrink.size),'methods':rows}

expected_m=(z['truth']*A).sum()/A.sum()
def illustrative(area):
    return {'area':float(area),'data_weight':float(area/(area+1)),
       'expected_shrunk_density_at_true_4':float((4*area+expected_m)/(area+1)),
       'equivalent_raw_threshold':float(2+(2-expected_m)/area)}
result['mechanism']={'expected_historical_regional_mean':float(expected_m),'small':illustrative(A[s]),'full':illustrative(A[f]),
    'small_count_examples':[{'count':n,'raw_density':float(n/A[s]),'shrunk_density':float((n+expected_m)/(A[s]+1))} for n in [1,2,3]]}
(save/'diagnostic.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
