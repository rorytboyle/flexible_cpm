# flexible_cpm
Modified connectome-based predictive model

1) Enables inclusion of covariates at feature selection or model-building stage or both stages
2) Permits different choices of k-fold cross-validation schemes
3) Returns weights (i.e. regression coefficients) from trained model for external validation on other datasets
4) Creates masks for visualisation in bioimagesuite and enables creation of masks when restricted connectivity matrices are used (i.e. when nodes are removed from the original 268 x 268 connectivity matrix) 

Functions with _prefix_ "_CPM___" are adapted/generalised functions of code written by Xilin Shen and Emily Finn and is uploaded here with their permission (codeshare_behavioral_prediction.m as obtained here: 
https://www.nitrc.org/frs/download.php/8071/nn_code.zip). This code is therefore copyrighted to Xilin Shen and Emily Finn as detailed below:

Copyright 2015 Xilin Shen and Emily Finn 
This code is released under the terms of the GNU GPL v2. This code is not FDA approved for clinical use; it is provided freely for research purposes. If using this in a publication please reference this properly as: 
Finn ES, Shen X, Scheinost D, Rosenberg MD, Huang, Chun MM,Papademetris X & Constable RT. (2015). Functional connectome fingerprinting: Identifying individuals using patterns of brain connectivity. Nature Neuroscience 18, 1664-1671.

All other functions and scripts (e.g. functions with _suffix_ "__CPM_") are original and this code is also released under the terms of the GNU GPL v2

Code is written and tested on MATLAB R2020a
