# ForamMgCa_PSM
R scripts and files associated with Saenger and Evans submission to Paleoceanography.

formMgCa_cal.val.R: constructs the sequential univariate and multivariate calibrations summarized in Figure 3. Calibrations are constructed from 75% of Mg/Ca data for each species and validated against withheld data. This requires the cal.val.R and cal.val.nls.R functions. Witheld predicted vs. observed values are used to calculate CE, which is summarized in Figure 3 of the manuscript.

formMgCa_all.R: constructs calibrations from all Mg/Ca data using the optimal independent variables identified in formMgCa.cal.val.R

formMgCaPSMfigs.R: generates manuscript figures.

cal.val.R: function for calibrating, and independently validating Mg/Ca relationships to multiple variables.

cal.val.nls.R: function for using non-linear least squares regression to calibrate, and independently validate Mg/Ca relationships to multiple variables.

ForamMg_culture.csv: compiled Mg/Ca data from foraminifera culture experiments.

ForamMg_Gbulloides, Ginflata, Gruber, Npachys: compiled Mg/Ca data used in core top calibrations. A separate file exists for each species (Gbulloides, Ginflata, Gruber, Npachys). Columns are core name, latitude, longitude, depth of core, measured Mg/Ca, oxygen isotope composition, oxygen isotope "calcification temperature", flag for reductive cleaning (Y=reductive cleaning was used), the mean size fraction analyzed (e.g. 300 = 250-350 um or 200-400 um, and original reference.

Gbulloides, Ginflata, Gruber, Npachys_summary: Summary of r2, Bayesian Information Criterion (BIC), RMSE, CE and their 95% confidence intervals for univariate/multivariate calibrations constructed using the variables in the final column. CE data from these files is summarized in Figure 3 of the manuscript.
