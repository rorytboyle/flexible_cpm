function [R_pos, P_pos, R_neg, P_neg, R_combined, P_combined,...
    rsq_pos, rsq_neg, rsq_combined, mae_pos, mae_neg, mae_combined]...
    = evaluate_CPM(all_behav, behav_pred_pos, behav_pred_neg,...
    behav_pred_combined)
% Evaluates accuracy/performance of CPM model. Calculates Pearson's
% correlation, R-squared/coefficient of determination, and MAE for each
% network strength model.
%
% Author: Rory Boyle
% Contact: rorytboyle@gmail.com
% Date: 24/01/2021
% Updated: 31/05/2021 to remove rsq_p_pos, rsq_p_neg, rsq_p_combined as
% values equivalent to p-value from Pearson's correlation
% Updated 15/09/2022 to correct MAE formula so that MAE is divided by sample size
%
%% 1)Pearson's correlations
[R_pos, P_pos] = corr(behav_pred_pos,all_behav);
[R_neg, P_neg] = corr(behav_pred_neg,all_behav);
[R_combined, P_combined] = corr(behav_pred_combined,all_behav);

%% 2) R^2 coefficient of determination
% recommended by Poldrack et al 2020 JAMA Psychiatry
pos_mdl = fitlm(behav_pred_pos, all_behav);
rsq_pos = pos_mdl.Rsquared.Ordinary;

neg_mdl = fitlm(behav_pred_neg,all_behav);
rsq_neg = neg_mdl.Rsquared.Ordinary;

combined_mdl = fitlm(behav_pred_combined,all_behav);
rsq_combined = combined_mdl.Rsquared.Ordinary;

%% 3) Mean Absolute Error (MAE)
mae_pos = (mean(abs(all_behav - behav_pred_pos)))/length(all_behav);
mae_neg = (mean(abs(all_behav - behav_pred_neg)))/length(all_behav);
mae_combined = (mean(abs(all_behav - behav_pred_combined)))/length(all_behav);
end
