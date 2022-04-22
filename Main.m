close all; clear; clc;
 
root_path = 'C:\Users\Cecilia\Desktop\Paper_LFPvsManifold\Código';
data_path = fullfile(root_path,'Data');
addpath(genpath(data_path));
% Add trial_data repo to path
addpath(genpath(fullfile(root_path,"TrialData")));

% Import trial_data structure
load(fullfile(data_path,'\filenames.mat'));
file = 10;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 

% Change original bin size (10 ms) to 30 ms
trial_data = binTD(trial_data,3);

% Trim data to window of interest
epoch = 'exec'; % exec(execution), feed(feedback processing), prep(preparation), rest
trial_data = trim_data(trial_data,epoch);


%%%%%%%%%%%%%%%
% MAIN ANALYSIS

% Canonical Correlation Analysis LFP vs Latent dynamics
pca_dims = 10; % Number of latent dynamics from PCA (M1: 10, PMd: 15, S1: 8)
[cca_coef,cca_coef_surr,p] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',pca_dims,'surrogate_iter',1,'doPlot',true));

% Velocity decoder
bins_to_past = 3; % Negative number means bins to future (for feedback control analysis)
[decoder_result] = vel_decoder(trial_data,struct('array',filenames{file,2},'folds',5,'bins_to_past',bins_to_past,'pca_dims',pca_dims,'doPlot',true,'eval',{{'r2','vaf'}}));

% Target classifier
[classifier_result] = target_classifier(trial_data,struct('array',filenames{file,2},'folds',5,'pca_dims',pca_dims,'doPlot',true));

% Unit activity correlation with LFPs 
[within_elect_corr, across_elect_corr] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',pca_dims,'doPlot',true));


%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL ANALYSIS

% LFP variance
[lfp_var] = lfp_variance(trial_data,struct('array',filenames{file,2},'doPlot',true));

% LFP-LFP correlation
[corr_data] = lfp_lfp_corr(trial_data,struct('array',filenames{file,2},'doPlot',true));

% LFP manifold vs latent dyanmics
lfp_dims = 10;
[cca_coef_supl,cca_coef_surr_supl,p_supl] = fCCA_LFPmanifold(trial_data,struct('array',filenames{file,2},'pca_dims',pca_dims,'lfp_dims',lfp_dims,'surrogate_iter',0,'doPlot',true));

% LFP power spectrum
% WARNING: This function works with the original dataset, not the trial_data structure
% Add ClassyDataAnalysis to working path
addpath(genpath(fullfile(root_path,"ClassyDataAnalysis")));
cds_path = fullfile(root_path,'CDS','Chewie_CO_CS_BL_10212016_001.mat');
[lfp_ps,ps_freqs] = lfp_power_spectrum(cds_path,struct('array',filenames{file,2},'epoch','exec','max_freq',400,'doPlot',true));

