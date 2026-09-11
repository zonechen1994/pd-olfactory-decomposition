"""Shared setup for the analysis scripts.

Every script starts with `from common import *`. This module fixes the working
directory to the repository root, resolves the PPMI source directory, imports the
libraries the analyses use, and provides `save_state` / `load_state` so that the
participant-level frames built by one script are available to the next without
being recomputed. State files are written under results/intermediate/ and are
never distributed (they contain participant-level PPMI data).
"""
import warnings, json, pickle, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import mean_absolute_error, roc_auc_score
from xgboost import XGBRegressor
from lightgbm import LGBMRegressor
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
np.random.seed(42)

# ---------------------------------------------------------------- paths
ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)                                   # all relative paths below are from the repo root
BASE = Path(os.environ.get('PPMI_DATA_DIR', ROOT))   # directory holding the PPMI source files
OUT = ROOT / 'results' / 'repro'
FIG = OUT / 'figures'
FIG_DIR = 'results/repro/figures'
MODELS = ROOT / 'model_output' / 'repro_sbr_models'
STATE = ROOT / 'results' / 'intermediate'
for d in (OUT, FIG, MODELS, STATE):
    d.mkdir(parents=True, exist_ok=True)

SBR_FILE = BASE / 'Xing_Core_Lab_-_Quant_SBR_18Mar2026.csv'
CLIN_FILE = BASE / 'PPMI_Curated_Data_Cut_Public_20260223.xlsx'
CLIN_SHEET = '20260105'
BIO_FILE = BASE / 'molecular' / 'Current_Biospecimen_Analysis_Results_24Mar2026.csv'
NULISA_CNS_FILE = BASE / 'molecular' / 'PPMI_Project_282_NULISAseq_CNSDiseasePanel_NPQCounts_20260120.xlsx'

plt.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.format': 'pdf', 'savefig.bbox': 'tight',
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'font.size': 10,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': False})
COL_3W = {'upsit': '#6b6b6b', 'OIS': '#2e3a5c', 'OMI': '#d4624a'}

# ---------------------------------------------------------------- state passing
def save_state(name, **objs):
    with open(STATE / f'{name}.pkl', 'wb') as f:
        pickle.dump(objs, f)
    print(f'-> state saved: results/intermediate/{name}.pkl ({", ".join(objs)})')


def load_state(name):
    p = STATE / f'{name}.pkl'
    if not p.exists():
        sys.exit(f'missing {p}. Run the earlier scripts first (see README).')
    with open(p, 'rb') as f:
        return pickle.load(f)


def make_pct_expected_putamen(norm_beta, age_col='age'):
    """Percent of age/sex-expected lowest putamen SBR (PPMI canonical metric, Chahine 2021).

    norm_beta = (intercept, beta_age, beta_sex) from the healthy-control regression fitted in
    02_survival_cohort.py. Identical to the function defined there; rebuilt here so that later
    scripts can use it without re-fitting.
    """
    def pct_expected_putamen(frame, age_override=None):
        a = (frame[age_col].astype(float) if age_override is None
             else pd.Series(np.asarray(age_override, dtype=float), index=frame.index))
        exp = norm_beta[0] + norm_beta[1] * a + norm_beta[2] * frame['sex_num'].astype(float)
        return frame[['PUTAMEN_L_REF_CWM', 'PUTAMEN_R_REF_CWM']].min(axis=1) / exp
    return pct_expected_putamen


def merge_ledger(sec):
    """Merge a dict of results into results/repro/section16_stats.json (create if absent)."""
    p = OUT / 'section16_stats.json'
    cur = json.load(open(p)) if p.exists() else {}
    cur.update(sec)
    with open(p, 'w') as f:
        json.dump(cur, f, indent=2, default=float)
    print(f'-> ledger results/repro/section16_stats.json: {len(cur)} keys')
