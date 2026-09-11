# TRIPOD+AI checklist

Collins GS, Moons KGM, Dhiman P, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ* 2024;385:e078378. doi:10.1136/bmj-2023-078378 (PMID 38626948)

Manuscript: *Decomposition of the olfactory score by dopamine transporter imaging improves Parkinson's disease risk stratification in hyposmic individuals*

Two notes on scope. The model developed here predicts a continuous behavioural score (UPSIT) from imaging, and its fitted value (OIS) is then evaluated as a prognostic marker for conversion to Parkinson's disease. Items phrased for a model whose target is the clinical outcome are answered against that two-step structure. Items marked **Not done** are stated as such in the manuscript rather than left silent.

| Item | Requirement (abbreviated) | Where reported |
|---|---|---|
| 1 | Title identifies model, target population, outcome | Title. Names the target population (hyposmic individuals) and the outcome (Parkinson's disease risk stratification). The words "prediction model" are not used, the decomposition being the subject of the paper |
| 2 | Abstract per TRIPOD+AI for Abstracts | Abstract. Data source, training set and size, model, application set with events, discrimination with 95% CI, and the negative control |
| 3a | Healthcare context, prognostic or diagnostic, rationale, existing models | Introduction. Prognostic. Current two-step practice and the gap it leaves. Competing continuous imaging staging is discussed in the Discussion [PMID 41917515] |
| 3b | Target population, intended purpose, intended users | Introduction and Discussion, clinical implications. Hyposmic individuals who have already had both examinations, read by the clinician alongside current practice rather than replacing it |
| 3c | Known health inequalities between sociodemographic groups | **Partly.** No inequality analysis was possible. The cohort composition (92% to 96% White, mean education 16 years) and the resulting untested transportability are stated in the Discussion, limitations |
| 4 | Study objectives, development or validation | Introduction, final paragraph. Development in Parkinson's disease and controls, then application unchanged to a prodromal cohort |
| 5a | Sources of data for development and evaluation, rationale, representativeness | Methods, Data and measurements. Methods, Model, for why the training set is Parkinson's disease plus controls. Representativeness in Discussion, limitations |
| 5b | Dates of participant data, accrual, end of follow-up | Methods, Data and measurements. Baseline visits July 2010 to December 2025, prodromal accrual from July 2013, latest visit December 2025, data cut 23 February 2026 |
| 6a | Setting, number and location of centres | Methods, Data and measurements. 44 PPMI centres in North America and Europe |
| 6b | Eligibility criteria | Methods, Participant flow and follow-up, with Supplementary Fig. 1 |
| 6c | Treatments received and how handled | Table 1 reports levodopa equivalent daily dose in the training cohort. Treatment was not entered as a covariate and no adjustment was made |
| 7 | Data pre-processing and quality checking | Methods, Data and measurements (scan selection, tracer, missing binding ratios, complete-case rule) and Methods, Model (standardisation) |
| 8a | Outcome definition, time horizon, how and when assessed, rationale | Methods, Participant flow and follow-up. First recorded primary diagnosis of Parkinson's disease. Two years is the primary horizon, with the follow-up distribution given as the reason |
| 8b | Qualifications of outcome assessors | Methods, Participant flow and follow-up. Site investigator within PPMI, not re-adjudicated here |
| 8c | Blinding of outcome assessment | Methods, Participant flow and follow-up. **Not possible.** Secondary analysis of recorded diagnoses |
| 9a | Choice of initial predictors and any pre-selection | Methods, Model. All 33 released binding ratios plus age, sex and years of education, with no pre-selection |
| 9b | Definition and measurement of predictors | Methods, Data and measurements |
| 9c | Qualifications of predictor assessors | Methods, Data and measurements. Central quantification by the PPMI imaging core laboratory [PMID 30564614] |
| 10 | How study size was arrived at, sample size calculation | Methods, Statistical analysis. **No a priori calculation.** Analysis sets fixed by the cohort and the stated exclusions. Events per analysis are given with every estimate |
| 11 | Handling of missing data, reasons for omissions | Methods, Data and measurements (complete case, UPSIT item imputation rule) and Participant flow and follow-up (stepwise exclusions with counts) |
| 12a | How data were used, partitioning | Methods, Model and Statistical analysis. Development in 1,398 participants, five-fold internal cross-validation, leave-one-centre-out over 44 centres, then unchanged application to 1,759 prodromal participants |
| 12b | Handling of predictors, rescaling, transformation | Methods, Model. Standardised on the training set, the fitted scaler applied unchanged elsewhere |
| 12c | Model type, rationale, building steps, hyperparameter tuning, internal validation | Methods, Model. Ridge, penalty fixed at 1.0 and not tuned, chosen among six families compared under the same five-fold cross-validation (Supplementary Table 1) |
| 12d | Heterogeneity across clusters | Methods, Statistical analysis. Leave-one-centre-out across all 44 centres (Fig. 2b, Supplementary Table 1) and a discovery against external centre split, whose inconsistency is reported |
| 12e | Performance measures and plots, and their rationale | Methods, Statistical analysis. Concordance index with bootstrap intervals, cumulative/dynamic AUC with inverse probability of censoring weighting, Kaplan-Meier with log-rank, calibration by decile (Supplementary Fig. 2) |
| 12f | Model updating or recalibration | Methods, Model. **None performed** in either cohort |
| 12g | How predictions were calculated for evaluation | Methods, Model, and Code availability. Fitted coefficients and standardisation parameters are released |
| 13 | Class imbalance methods | **Not applicable.** The model target is a continuous score, not a class |
| 14 | Approaches to model fairness | **Not done.** Stated in the Discussion, limitations, with the reason |
| 15 | Model output, classification and thresholds | Methods, Model, and Results, risk stratification. Output is a continuous score in UPSIT points. Cut points are reported as exploratory and are not proposed as clinical thresholds |
| 16 | Differences between development and evaluation data | Table 1 places the training cohort and the two prodromal groups side by side. Discussed in Results, study design |
| 17 | Ethics approval and consent | Methods, Ethics |
| 18a | Funding and role of funders | Funding |
| 18b | Conflicts of interest | Competing interests |
| 18c | Study protocol | Methods, Ethics. **No protocol was prepared in advance.** Analyses are specified in the released scripts |
| 18d | Registration | Methods, Ethics. **Not registered** |
| 18e | Availability of data | Data availability |
| 18f | Availability of analytical code | Code availability |
| 19 | Patient and public involvement | Methods, Ethics. **None** |
| 20a | Flow of participants, outcome counts, follow-up summary | Methods, Participant flow and follow-up, and Supplementary Fig. 1 |
| 20b | Characteristics by data source, key dates, predictors, events, follow-up, missing data | Table 1, with per-variable denominators in the released table. Follow-up by recruitment cohort in Supplementary Table 2 |
| 20c | Comparison of predictor distributions, development against evaluation | Table 1 |
| 21 | Participants and events in each analysis | Given with every estimate in Results, in Tables 2 to 4 and in the Supplementary tables |
| 22 | Full model available for prediction in new individuals | Code availability. Coefficients, standardisation parameters and the fitted object |
| 23a | Performance with confidence intervals, including subgroups | Tables 2 to 4, Fig. 3, Supplementary Table 3. Every concordance index carries a 95% CI |
| 23b | Heterogeneity in performance across clusters | Fig. 2b and Supplementary Table 1, per-centre correlations. The centre split result is reported as heterogeneity and not as support |
| 24 | Results of model updating | **Not applicable.** No updating was performed |
| 25 | Overall interpretation | Discussion |
| 26 | Limitations and their effects | Discussion, final two paragraphs. Proportional hazards checked by Schoenfeld residuals (Supplementary Table 14) |
| 27a | Handling of poor quality or unavailable input data at implementation | Discussion, clinical implications. Regional binding ratios are required and no imputation is proposed where a scan fails quality control or is quantified only as a single region |
| 27b | User interaction and expertise required | Discussion, clinical implications. The score is computed from an existing quantified scan and three routinely recorded covariates, and is read in the points of the olfactory test |
| 27c | Next steps for research, applicability and generalisability | Discussion, limitations. External validation in an independent cohort is the stated next step |
