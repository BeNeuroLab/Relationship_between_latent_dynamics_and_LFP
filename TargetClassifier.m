% function [result,result_boot] = TargetClassifier(trial_data,params)
%
%   Will generate target classifiers and return the performance
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) 'M1' or 'PMd' depending on the desaired array
%     .folds      : (int) Number of classifiers to compute
%     .doPlot     : (logical) Whether to visualize result
%
% Written by Cecilia Gallego. Updated April 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result,result_boot] = TargetClassifier(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
folds   = 20;
pca_dims  = 8;
doPlot = false;
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

% Smooth the signals and compute neural modes
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array 'lfp'],'width',0.05));
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));


% Downsample (one sample per feature)
td = binTD(trial_data,'average');

targets = [td.tgtDir];

lfp_data = cell2mat({td.([array '_lfp'])}');
pca_data = cell2mat({td.([array '_pca'])}');


result = zeros(folds,length(freqs)+1);
result_boot = zeros(folds,length(freqs)+1);
for b = 1:(length(freqs)+1)    
    
    disp(['Signal ' num2str(b) ' out of ' num2str(length(freqs)+1)])
    
    if b == length(freqs)+1
        signal = pca_data;
    else
        isBand = lfp_guide(:,3) == freqs(b);
        lfp_idx = isBand & lfp_ch;
        signal = lfp_data(:,lfp_idx);
    end
    
    % z-transform to normalize the data
    signal_mean = repmat(mean(signal),size(signal,1),1);
    signal_std = repmat(std(signal),size(signal,1),1);
    norm_signal = (signal-signal_mean)./signal_std;
    
    % Eliminate nan data
    bool = sum(isnan(norm_signal));
    norm_signal(:,bool~=0) = [];
    
    for n = 1:folds
%         disp(['Fold ' num2str(n) ' out of ' num2str(folds)]);
        idx = randperm(length(td));
        test_idx = idx(1:round(length(td)*0.1));
        train_idx = idx(round(length(td)*0.1)+1:end);
        
        O1 = fitcnb(norm_signal(train_idx,:),targets(:,train_idx)');
        prediction = O1.predict(norm_signal(test_idx,:));

        result(n,b) = sum(targets(:,test_idx)' == prediction)/length(prediction);
        
        boot_idx = randperm(length(train_idx));
        boot_train_idx = train_idx(boot_idx);
        
        O1 = fitcnb(norm_signal(boot_train_idx,:),targets(:,train_idx)');
        prediction = O1.predict(norm_signal(test_idx,:));
        
        result_boot(n,b) = sum(targets(:,test_idx)' == prediction)/length(prediction);
    end
end

if doPlot    
    x = reshape([result; result_boot],[],20);

    pos1 = 1:3:30; pos2 = 2:3:31; positions = sort([pos1 pos2]);
    label_positions = (pos1+pos2)./2; 
   

    figure
    h = boxplot(x, 'positions', positions, 'color','k','symbol','k+');  
    set(h,'linewidth',1) 
    
    load('bands_name.mat'); bands_name{end+1} = 'Lat var';
    set(gca,'xtick',label_positions,'xticklabel',bands_name) 
    
    c = parula(18); c([19,20],:) = [0 0 0; 0 0 0];
    c = flipud(c); c(1:2:18,:) = 0;
    h = findobj(gca,'Tag','Box'); 
    for j=1:length(h) 
    p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
    end
    
    ylim([0 1])
    legend([p(4) p(5)],{'Regular','Bootstraping'})

    ylabel('Classifier performance');

end

end