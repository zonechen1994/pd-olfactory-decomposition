"""Reproduce the current analyses from licensed PPMI inputs.

Run ``python code/run_analysis.py`` from the repository root after setting
PPMI_DATA_DIR as described in README.md.
"""
from pathlib import Path
import os,subprocess,sys
R=Path(__file__).resolve().parent.parent;os.chdir(R)
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
for name in ['01_fit_model.py','02_survival_cohort.py','03_discrimination.py','04_residual.py','05_baseline_table.py','06_ph_test.py','07_channel_followup.py']:
    subprocess.run([sys.executable,str(R/'code'/name)],check=True)
A=R/'code/validated'
for local,target in [('results/intermediate',R/'results/intermediate'),('model_output',R/'model_output')]:
    p=A/local;p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():p.symlink_to(target,target_is_directory=True)
subprocess.run([sys.executable,str(A/'audit_analysis.py')],check=True)
subprocess.run([sys.executable,str(A/'check_secondary.py')],check=True)
(A/'figures').mkdir(exist_ok=True)
os.environ['PANEL_CASE']='upper';os.environ['SUPP_FIG_DIR']=str(A/'figures')
subprocess.run([sys.executable,str(A/'mk_main_figures.py'),'en',str(A/'figures')],check=True)
subprocess.run([sys.executable,str(A/'mk_supplementary.py')],check=True)
subprocess.run([sys.executable,str(A/'mk_source_data.py')],check=True)
print('Validated results are in code/validated/results/repro; figures in code/validated/figures.')
