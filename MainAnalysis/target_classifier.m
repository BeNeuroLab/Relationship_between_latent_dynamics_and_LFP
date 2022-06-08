% function [result] = target_classifier(trial_data,params)
%
%   Will generate target classifiers and return the performance
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) 'M1' or 'PMd' depending on the desaired array
%     .folds      : (int) Number of classifiers to compute
%     .doPlot     : (logical) Whether to visualize result
%     .pca_dims   : (int) Number of PCA dimensions
%
% OUTPUTS:
%   result: (array) #folds x #bands+1 array where each entry is the result
%           of one classsifier, each column corresponds to a frequency band
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = target_classifier(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
folds   = 20;
pca_dims  = 8;
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiking channel (we are only using LFPs from channels with spiking activity)
spk_ch = unique(trial_data(1).([array '_unit_guide'])(:,1));
lfp_ch = unique(trial_data(1).([array '_lfp_guide'])(:,1));
spk_ch = spk_ch(ismember(spk_ch,lfp_ch)); % Check that there is an LFP in each spiking channel

% Smooth the signals
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));

% Do PCA to obtain the latent dynamics
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));

% Get average activity per trial (one time bin per trial)
td = binTD(trial_data,'average');

% Get bands max freq and number of bands
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
bands_number = length(max_fq);


% Extract signals from trial_data structure
lfp_data = cell2mat({td.([array '_lfp'])}');
pca_data = cell2mat({td.([array '_pca'])}');
targets = [td.tgtDir];


result = zeros(folds,bands_number+1);
for iband = 1:bands_number+1   

    disp(['Signal ' num2str(iband) ' out of ' num2str(bands_number)])
    
    if iband == bands_number+1
        signal = pca_data;
    else
        lfp_idx = find(trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(iband));
        lfp_idx = lfp_idx(spk_ch);
        signal = lfp_data(:,lfp_idx);
    end
    
    % z-transform to normalize the data
    signal_mean = repmat(mean(signal),size(signal,1),1);
    signal_std = repmat(std(signal),size(signal,1),1);
    norm_signal = (signal-signal_mean)./signal_std;
    
    % Eliminate nan data
    bool = sum(isnan(norm_signal));
    norm_signal(:,bool~=0) = [];
    
    for ifold = 1:folds
        % Randomly select test and train trials
        idx = randperm(length(td));
        test_idx = idx(1:round(length(td)*0.1));
        train_idx = idx(round(length(td)*0.1)+1:end);
        
        % Fit and test the target classifier
        O1 = fitcnb(norm_signal(train_idx,:),targets(:,train_idx)');
        prediction = O1.predict(norm_signal(test_idx,:));
        
        % Save the result as number of correct predictions / all test trials
        result(ifold,iband) = sum(targets(:,test_idx)' == prediction)/length(prediction);
    end
end

if doPlot    

    % If there are names for the frequency bands
    if exist('bands_name.mat','file') == 2
        load('bands_name.mat');
    else % Else, assign ordered numbers
        bands_name = num2cell(1:bands_number);
    end
    bands_name{bands_number+1} = 'Lat Var';
    

%     figure
%     boxplot(result,'color','k','symbol','k+');
%     c = flipud(parula(bands_number)); c = cat(1,[0 0 0],c);
%     h = findobj(gca,'Tag','Box'); 
%     for j=1:length(h) 
%     p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
%     end 
%     % Plot chance level line
%     hold on; line(get(gca,'xlim'),[1 1]/length(unique(targets)),'color','k')
%     title(trial_data(1).monkey);
%     ylabel('Classifier performance'); ylim([0 1])
%     set(gca, 'XTick', 1:length(bands_name), 'XTickLabel',bands_name);
    
    figure
    c = parula(bands_number); c(bands_number+1,:) = [0 0 0];
    for n = 1:bands_number+1
        hold on; scatter(ones(size(result,1),1)*n-.25,result(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
        hold on; scatter(n,median(result(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
        hold on; errorbar(n,median(result(:,n)),std(result(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
    end
    hold on; line([0 bands_number+1.5],[1 1]/length(unique(targets)),'linewidth',1.5,'color',[0 0 0]);
    
    set(gca, 'XTick', 1:length(bands_name), 'XTickLabel',bands_name);
    
    xlim([0.5 bands_number+1.5]); ylim([0 1])
    ylabel('Classifier performance');
    set(gca,'TickDir','out'); box off
end

end