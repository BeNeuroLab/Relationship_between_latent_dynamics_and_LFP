%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [total_result] = LFP_Decoder(trial_data,params)
%
%   Will generate velocity decoders based on individual frequency bands and
%   return the VAF of the prediction together with the actual and predicted
%   velocities
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) array location (e.g. 'M1')
%     .BinToPast  : (int) Number of bins of history for the model, after downsampling
%     .folds      : (int) How many decoders to compute
%     .doPlot     : (logical) Whether to visualize result
%     .eval       : (cell) Indicate evaluation method (VAF, R-square, or both
%
% Written by Cecilia Gallego. Updated April 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [total_result] = LFP_Decoder(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
BinToPast = 5;
folds = 20;
eval = {'vaf'};
doPlot = true;
pca_dims = 8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiking channel (we are only using LFPs from channels with spiking activity)
spk_ch = unique(trial_data(1).([array '_unit_guide'])(:,1));
lfp_guide = trial_data(1).([array '_lfp_guide']);
freqs = unique(lfp_guide(:,3));
lfp_ch = ismember(lfp_guide(:,1),spk_ch);

% Smooth the signals
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'calc_fr',true,'width',0.05));               
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'calc_fr',true,'width',0.05));               

% Compute neural modes
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));


% Single frquencies
load('bands_name.mat');
result = cell(folds,length(freqs)+1);
for b = 1:length(freqs)
    
    disp(' '); disp(['Band ' num2str(b) ' out of ' num2str(length(freqs))])
    aux_trial_data = trial_data;
    isFreq = lfp_guide(:,3) == freqs(b);
    lfp_idx = isFreq & lfp_ch;
    for n = 1:length(trial_data)
        one_band_lfp = trial_data(n).([array '_lfp'])(:,lfp_idx);
        aux_trial_data(n).([array '_lfp']) = one_band_lfp;
    end
    for n = 1:folds
        disp(['Fold ' num2str(n) ' out of ' num2str(folds)])
        result{n,b} = TD_ComputeModel(aux_trial_data,BinToPast,[array '_lfp'],eval);
        result{n,b}.band = b;
    end
end
% Latent varibles
disp(' '); disp('Latent var'); band_names{b+1} = 'Lat Var';
name = [array '_pca']; 
for n = 1:folds
    disp(['Fold ' num2str(n) ' out of ' num2str(folds)])
    result{n,b+1} = TD_ComputeModel(aux_trial_data,BinToPast,name,eval);
    result{n,b+1}.band = b+1;
end

total_result = [result{:}];
bands = [total_result.band];

if doPlot
    
    if length(eval) == 2
        vaf = []; r2 = [];
        for n = 1:length(unique(bands))
            idx = bands == n;
            vaf(:,n) = [total_result(idx).vaf]';
            r2(:,n) = [total_result(idx).r2]';
        end

        figure
        subplot(1,2,1)
        boxplot(vaf,'color','k','symbol','k+');
        c = flipud(parula(9)); c = cat(1,[0 0 0],c);
        h = findobj(gca,'Tag','Box'); 
        for j=1:length(h) 
        p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
        end 
        title(trial_data(1).monkey); %ylim([-0.2,1]);
        ylabel('Decoder performance (VAF)');
        set(gca, 'XTick', 1:10, 'XTickLabel',band_names);

        subplot(1,2,2)
        boxplot(r2,'color','k','symbol','k+');
        c = flipud(parula(9)); c = cat(1,[0 0 0],c);
        h = findobj(gca,'Tag','Box'); 
        for j=1:length(h) 
        p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
        end 
        title(trial_data(1).monkey); %ylim([-0.2,1]);
        ylabel('Decoder performance (r2)');
        set(gca, 'XTick', 1:10, 'XTickLabel',band_names);
        
        
    else
        vaf = [];
        for n = 1:length(unique(bands))
            idx = bands == n;
            vaf(:,n) = [total_result(idx).(eval{1})]';
        end
        
        figure
        boxplot(vaf,'color','k','symbol','k+');
        c = flipud(parula(9)); c = cat(1,[0 0 0],c);
        h = findobj(gca,'Tag','Box'); 
        for j=1:length(h) 
        p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
        end 
        title(trial_data(1).monkey); %ylim([-0.2,1]);
        ylabel(['Decoder performance (' eval{1} ')']);
        set(gca, 'XTick', 1:10, 'XTickLabel',band_names);
    end
end


 
end

