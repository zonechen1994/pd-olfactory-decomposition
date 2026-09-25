"""Run from the extracted package root after placing the licensed PPMI inputs there."""
from pathlib import Path
import os,subprocess,sys
R=Path(__file__).resolve().parent;os.chdir(R)
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
for name in ['01_fit_model.py','02_survival_cohort.py','03_discrimination.py','04_residual.py','05_baseline_table.py','06_ph_test.py','07_channel_followup.py']:
    subprocess.run([sys.executable,str(R/'code'/name)],check=True)
A=R/'audit/JTM_20260924'
subprocess.run([sys.executable,str(A/'audit_analysis.py')],check=True)
subprocess.run([sys.executable,str(A/'check_secondary.py')],check=True)
for local,target in [('results/intermediate',R/'results/intermediate'),('model_output',R/'model_output')]:
    p=A/local;p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():p.symlink_to(target,target_is_directory=True)
(A/'figures').mkdir(exist_ok=True)
os.environ['PANEL_CASE']='upper';os.environ['SUPP_FIG_DIR']=str(A/'figures')
subprocess.run([sys.executable,str(A/'821_results/mkfig234.py'),'en',str(A/'figures')],check=True)
subprocess.run([sys.executable,str(A/'mk_supplementary.py')],check=True)
subprocess.run([sys.executable,str(A/'mk_source_data.py')],check=True)
print('Corrected JTM results are in audit/JTM_20260924/results/repro; figures in its figures/ directory.')
