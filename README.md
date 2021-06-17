# Data Analysis Code For MRT Paper Submitted To Psychological Methdos

Tianchen Qian

analysis.R is the file conducts the analysis.

xgeepack.R is needed as it contains the function that implements the WCLS estimator.

synthetic_data/ contains a synthetic data set (mimicking some features of HeartSteps V1 data), as well as example analysis code which implements the WCLS estimator. The synthetic data set is "synthetic_data_37subject_210time.csv". See Section 1 of "analysis_synthetic_data_with_results.pdf" for explanation of each variable in the synthetic data set.

(Note: For SAS users, please download "synthetic_data_37subject_210time_SAS.csv" instead, where the variables are renamed to conform with SAS's variable naming requirement.)
