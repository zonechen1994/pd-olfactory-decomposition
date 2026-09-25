# Decomposition of the olfactory score by dopamine transporter imaging

Analysis code for the manuscript *Decomposition of the olfactory score by dopamine transporter imaging improves Parkinson's disease risk stratification in hyposmic individuals* (PPMI cohort).

**Current JTM submission:** The corrected comparison and landmark analyses, their aggregate audit outputs, and instructions for reproducing the submitted results are documented in [README_JTM.md](README_JTM.md). Run [RUN_JTM.py](RUN_JTM.py) to generate the submitted analysis from licensed PPMI inputs. The original `code/` scripts and root `results/repro/` remain as a base snapshot; use the corrected outputs under `audit/JTM_20260924/results/repro/` for the submission. The committed `Source_Data.xlsx` matches the submitted Supplementary Material 2.

The UPSIT olfactory score is regressed on the 33 striatal dopamine transporter binding ratios released by the PPMI imaging core laboratory, plus age, sex and years of education, by ridge regression fitted in 1,398 participants with Parkinson's disease or healthy control status. The fitted value is the olfactory imaging score (OIS) and the residual is the olfactory mismatch index (OMI), so that UPSIT = OIS + OMI. OIS is then evaluated for conversion to Parkinson's disease in 1,759 prodromal PPMI participants (152 conversions), and the residual is characterised against CSF tau and GBA1 against LRRK2 genotype.

## What is here

| Path | Content |
|---|---|
| `code/01_fit_model.py` | Loads the two PPMI files, compares six regression families under five-fold cross-validation, fits the ridge model, writes OIS and OMI for every participant, leave-one-centre-out validation over 44 centres |
| `code/02_survival_cohort.py` | Prodromal survival frame (1,759 / 152), temporal and centre splits, age/sex-expected putamen metric validated against the official PPMI staging field, test-retest ICC, landmark analyses, comparison with Cox models trained on conversion labels, calibration, variance decomposition. Starts the numerical ledger |
| `code/03_discrimination.py` | Pre-specified paired comparisons, interaction test, concordance with bootstrap intervals for every stratum and reading (Tables 2 to 4), continuous against dichotomised readings, calibration and age within strata, group flow, median split (Fig. 3d, 3e), permutation-corrected three-group split, IPCW time-dependent AUC (Fig. 3b) |
| `code/04_residual.py` | Genotype contrast (Fig. 4b, 4c), staging axes, CSF tau and neurofilament (Fig. 4a), residual by biological stage, training-cohort alternatives, prodromal seed amplification, genotype premise, cognitive endpoints, age in the residual, partial correlations and tau tertile means (Fig. 1d) |
| `code/05_baseline_table.py` | Table 1 |
| `code/06_ph_test.py` | Schoenfeld proportional-hazards test for the reported Cox models (Supplementary Table 5) |
| `code/07_channel_followup.py` | Follow-up by recruitment cohort (Supplementary Table 2) |
| `code/08_main_figures.py` | Figures 2 to 4 |
| `code/09_fig1_panels.py` | The three data panels of Fig. 1d (the rest of Fig. 1 is artwork) |
| `code/10_supplementary.py` | `Supplementary_Information_EN.md` and Supplementary Figs 1 to 5 |
| `code/11_source_data.py` | `Source_Data.xlsx`, one sheet per figure panel with the values each panel draws |
| `code/common.py` | Paths, shared imports, state passing between scripts |
| `results/repro/section16_stats.json` | The original numerical ledger, retained as the base snapshot |
| `results/repro/*.csv`, `*.json` | Original aggregate tables, retained as the base snapshot |
| `audit/JTM_20260924/` | Corrections and checks for the submitted analysis; see `README_JTM.md` |
| `audited_aggregate_outputs/` | Aggregate results and checks from the corrected server run |
| `model_output/repro_sbr_models/` | The fitted ridge model, its standardisation scaler and feature metadata (all that is needed to compute OIS in new individuals) |
| `figures/main/`, `figures/Supplementary/` | Figures as generated (PDF with editable Type 1 text, plus PNG) |
| `Supplementary_Information_EN.md` | Supplementary Information as generated |
| `Source_Data.xlsx` | Source data for Figs 1d, 2, 3, 4 and Supplementary Figs 1 to 5 (aggregate values only) |
| `TRIPOD_AI_checklist.md` | Completed TRIPOD+AI checklist |
| `data/README.md` | The PPMI source files required, with checksums |

## Data

No PPMI data are included. The PPMI Data Use Agreement does not permit redistribution, including of participant-level derived values. Apply for access at www.ppmi-info.org/access-data-specimens/download-data, download the files listed in `data/README.md`, and set `PPMI_DATA_DIR` to the directory holding them (default: the repository root).

The manuscript states this as follows.

> **Data availability.** All data used in this study are from the Parkinson's Progression Markers Initiative (PPMI) database (www.ppmi-info.org/access-data-specimens/download-data, RRID:SCR_006431) and are available on application. Clinical data are the curated data cut of 23 February 2026 and imaging quantification is the core laboratory binding-ratio file of 18 March 2026. Every derived value reported here is released with the code as an aggregate. The PPMI Data Use Agreement does not permit redistribution of participant-level data, including participant-level derived values such as the two components, and those are regenerated by the released scripts from the source files.

> **Code availability.** The base analysis code, scripts and aggregate audit outputs for the corrected analyses, the fitted ridge model and standardisation parameters are available in this repository. Aggregate figure source data are available here and in the submitted Supplementary Material 2. The corrected analyses are described in `README_JTM.md`. An archived version with a DOI will be deposited on acceptance.

## Running

The commands below reproduce the original base analysis. For the submitted corrected analysis, follow `README_JTM.md` and use `python RUN_JTM.py` instead.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export PPMI_DATA_DIR=/path/to/ppmi/files

python code/01_fit_model.py          # about 1 min
python code/02_survival_cohort.py    # about 3 min
python code/03_discrimination.py     # about 15 min (bootstrap and permutation loops)
python code/04_residual.py           # about 2 min (reads the 177 MB biospecimen file)
python code/05_baseline_table.py
python code/06_ph_test.py
python code/07_channel_followup.py
python code/08_main_figures.py
python code/09_fig1_panels.py
python code/10_supplementary.py
python code/11_source_data.py
```

The scripts are run in this order. Scripts 01 and 02 write participant-level frames to `results/intermediate/` that the later scripts read, so they only need to be run once. All paths are resolved from the repository root whatever the working directory. Random seeds are fixed (five-fold splitting 41, bootstrap 42, centre split 14). On the PPMI release listed in `data/README.md` the scripts reproduce the base `results/repro/section16_stats.json`; the submitted corrected results require the additional audit step.

Scripts 01 to 05 are the cells of the analysis notebook in which this work was developed, in their original order and with the original section numbers in the comments, unchanged except for file paths. The notebook also contained analyses that were later withdrawn and are not reported in the manuscript. Those cells are not included, and the ledger keys they would have written are absent.

## Computing OIS in new individuals

```python
import json, pickle
d = 'model_output/repro_sbr_models/'
model  = pickle.load(open(d + 'sbr_all_Ridge.pkl', 'rb'))   # sklearn Ridge, alpha = 1.0
scaler = pickle.load(open(d + 'scaler.pkl', 'rb'))          # StandardScaler fitted on the training cohort
meta   = json.load(open(d + 'metadata.json'))
X = df[meta['feature_names']].values   # the 33 binding ratios in this order, then age (years), sex (1 male, 0 female), years of education
OIS = model.predict(scaler.transform(X))
OMI = df['upsit'].values - OIS
```

The binding ratios must be the PPMI imaging core laboratory quantities (region minus reference over reference, cerebellar white matter reference), with the earliest scan per participant. Scores are in UPSIT points and are not constrained to 0 to 40.

## Citation

To be added on publication.

## License

MIT, see `LICENSE`.
