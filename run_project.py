"""Verify immutable outputs or reproduce into a fresh directory."""
from pathlib import Path
import argparse,hashlib,json,sys,shutil,subprocess,csv

ROOT=Path(__file__).resolve().parent

def normalized(root,relative):return root.joinpath(*relative.replace('\\','/').split('/'))

def verify():
    count=0
    manifests=[(ROOT/'experiments/h3_coastal_experiment',ROOT/'experiments/h3_coastal_experiment/results/stage1/manifest.json','sha256'),(ROOT/'experiments/h3_continuous_checks',ROOT/'experiments/h3_continuous_checks/manifest.json','files_sha256')]
    for base,path,key in manifests:
        for relative,expected in json.loads(path.read_text('utf-8'))[key].items():
            p=normalized(base,relative)
            if not p.is_file():raise RuntimeError(f'Missing original file: {p}')
            if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise RuntimeError(f'Hash mismatch: {p}')
            count+=1
    rows=json.loads((ROOT/'experiments/h3_continuous_checks/results/summary.json').read_text('utf-8'))
    for r in rows:
        for target in ['small','full']:
            if abs(sum(r[target+'_'+s] for s in ['confirmed','pending','negative'])-1)>1e-10:raise RuntimeError('State accounting failed')
    result={'passed':True,'original_fingerprints_verified':count,'continuous_summary_rows':len(rows),'checks':['original experiment output/source fingerprints','three-state accounting at all looks']}
    print(json.dumps(result,ensure_ascii=False,indent=2))

def reproduce(destination):
    destination=destination.resolve()
    if destination.exists():raise RuntimeError('Choose a new output directory; existing results are preserved.')
    destination.mkdir(parents=True)
    scratch=destination/'experiments'
    src1=ROOT/'experiments/h3_coastal_experiment'
    one=scratch/'h3_coastal_experiment'
    (one/'src').mkdir(parents=True)
    for name in ['coastal_experiment.py','analyze_and_plot.py']:shutil.copy2(src1/'src'/name,one/'src'/name)
    for name in ['EXPERIMENT_PLAN.md','requirements.txt']:shutil.copy2(src1/name,one/name)
    subprocess.run([sys.executable,str(one/'src/coastal_experiment.py')],check=True)
    subprocess.run([sys.executable,str(one/'src/analyze_and_plot.py')],check=True)
    src2=ROOT/'experiments/h3_continuous_checks'
    two=scratch/'h3_continuous_checks';two.mkdir()
    for name in ['PLAN.md','run_checks.py']:shutil.copy2(src2/name,two/name)
    subprocess.run([sys.executable,str(two/'run_checks.py')],check=True)
    shutil.copy2(src2/'verify_and_plot.py',two/'verify_and_plot.py')
    subprocess.run([sys.executable,str(two/'verify_and_plot.py')],check=True)
    # The archived post-hoc diagnostic originally used workspace-relative paths.
    # Recreate that small layout in the new run rather than editing its source.
    facade=destination/'threshold_workspace'
    shutil.copytree(one,facade/'outputs/h3_coastal_experiment')
    (facade/'work').mkdir()
    diagnostic=ROOT/'experiments/threshold_diagnostic/original_threshold_diagnostic.py'
    subprocess.run([sys.executable,str(diagnostic)],cwd=facade,check=True)
    # Compare deterministic numeric outputs rather than timestamps or archive bytes.
    old1=json.loads((src1/'results/stage1/summary.json').read_text('utf-8'))
    new1=json.loads((one/'results/stage1/summary.json').read_text('utf-8'))
    old2=json.loads((src2/'results/summary.json').read_text('utf-8'))
    new2=json.loads((two/'results/summary.json').read_text('utf-8'))
    old3=json.loads((ROOT/'experiments/threshold_diagnostic/diagnostic.json').read_text('utf-8'))
    new3=json.loads((facade/'work/threshold_check/diagnostic.json').read_text('utf-8'))
    result={'stage1_summary_identical':old1==new1,'stage2_summary_identical':old2==new2,'threshold_diagnostic_identical':old3==new3,'output_directory':str(destination)}
    (destination/'reproduction_check.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    if not all([result['stage1_summary_identical'],result['stage2_summary_identical'],result['threshold_diagnostic_identical']]):raise RuntimeError('Reproduction summaries differ; inspect the output and environment.')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('verify');r=sub.add_parser('reproduce');r.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.command=='verify':verify()
    else:reproduce(a.output)
