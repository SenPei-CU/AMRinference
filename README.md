# AMRinference

Code for running a synthetic test in Pei S. et al. Synergistic use of multimodal observations improves identification of asymptomatic carriers of antimicrobial-resistant organisms (2025)

## Data

contactnetwork.mat: The daily time-varying contact network. (1) patlist: the list of patient IDs on each day. (2) adm: the first admission date for each patient. (3). readmittedpat: patient IDs of readmitted patients on each day. (4) nl (neighbor list) and (5) part (partition) jointly specify the contact network structure on each day. Specifically, nl{t} is the list of contact IDs on day t. part(:,t) records the rows of contacts for each patient in nl{t}. The neighbors of patient I on day t are nl{t}(part(i):part(i+1)-1,t).

syntheticdata.mat: Data for a synthetic outbreak. (1) Y: clinical culture results. Three columns - sample collection time, patient ID, test result (1 positive; 0 negative)
(2) M: transmission matrix. Mij=1 if WGS data suggest transmission from i to j; Mij=0 otherwise. (3) xrec: the actual status of patients on each day. (4) gamma_truth: average importation rate in model simulation. (5) beta_truth: transmission rate in model simulation. (6) alpha_truth: decolonization rate in model simulation. (7) rho: effective sensitivity in model simulation.

admdis.mat: The range of hospitalization history for each patient. The earliest admission date and last discharge date for each patient. Three columns - patient ID, earliest admission data, last
%discharge date. Note that each patient may have multiple hospitalizations.

neighborinfo.mat: All contacts for each patient during the study period. nlall and partall jointly specify the information about contacts for each patient. nlall has three columns - neighbor ID, contact start date, contact end date. partall records the rows in nlall corresponding to each patient. For instance, the contacts of patient i are nlall(partall(i):partall(i+1)-1,1)

hospduration.mat: The period of hospitalization for each patient. Three columns - patient ID, admission date, discharge date. Note that each row is for one single hospitalization and one patient may have multiple rows for several hospitalizations.

medicalinfo.mat: Individual patient characteristics for an example population. Variables for each patient are shown as the column names.

regressionmodel_select.mat: Regression model to estimate initial colonization probability upon first admission. (1) bs: coefficients for covariates. (2) variables: the selected 14 variables in the regression model.

## Functions

harmonizeM.m - The function to translate the transmission matrix M informed by WGS to colonization probabilities of patients on paths connecting genetically related patients.

inference.m - The function to estimate individual colonization probabilities using information from patient characteristics, culture test results, and transmission events identified using WGS data.

## Run codes

To run the code, first compile the C++ programs in MATLAB using "mex ABMsimulation_chain_p_beta.c" and "mex masterequation_p_beta.c"

Run harmonizeM() to estimate colonization probabilities for patients on paths connecting genetically related patients.

Run inference() to estimate colonization probabilities for all patients using information from multimodal observations.
