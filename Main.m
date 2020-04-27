close all; clear; clc

array = 'M1';
filename =  ['Mihili_CO_20140303_' array '.mat'];


trial_data = loadTDfiles(filename,{@getTDidx,{'result','R'}});

switch array
    case 'PMd'
        pca_dims = 12; 
        trial_data = trimTD(trial_data,{'idx_movement_on',-39},{'idx_movement_on',6});
    case 'M1'
        pca_dims = 8; BinToPast = 5;
        trial_data = trimTD(trial_data,{'idx_movement_on',-12},{'idx_movement_on',42});
    case 'LeftS1Area2'
        pca_dims = 6; BinToPast = -5;
        trial_data = trimTD(trial_data,{'idx_movement_on',0},{'idx_movement_on',60});
end



% Canonical Correlation Analysis
[cca_coef,cca_bcoef,p] = fCCA(trial_data,struct('array',array,'pca_dims',pca_dims,'boot_iter',3,'boot_type','TME','doPlot',true));


% Spikes and LFP correlation
[correlation,boot_correlation] = ArrayCorrelation(trial_data,struct('array',array,'doPlot',true));


% Velocity decoders
[result] = LFP_Decoder(trial_data,struct('array',array,'folds',50,'BinToPast',BinToPast,'doPlot',true,'eval',{{'vaf','r2'}}));


% Target classifier
[result_classif,result_boot] = TargetClassifier(trial_data,struct('array',array,'folds',50,'doPlot',true));
