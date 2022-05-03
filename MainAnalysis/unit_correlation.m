%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [within_elect, across_elect] = unit_correlation(trial_data,params)
%
%   Will compute the Pearson's correlation between the unit spiking
%   form a channel and the LFP signal at different frequencies form the
%   same channel (within_elect) or different channel (across_elect)
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) array location (e.g. 'M1')
%     .signal     : (char) unit signal to use. It can be SUA (single unit activity), 
%                   MUA (multi unit activity), denoised (PCA-denoised SUA)
%     .pca_dims   : (int) Number of PCA dimensions
%     .doPlot     : (logical) Whether to visualize result
%
% OUTPUTS:
%   within_elect: (array) #channels x #bands array where each entry is the 
%           correlation of the unit activity in channel (row) with the same 
%           channel LFP activity in the frequency band (column)
%   across_elect: (array) #channels x #bands array where each entry is the 
%           correlation of the unit activity in channel (row) with another 
%           channel LFP activity in the frequency band (column)
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [within_elect, across_elect] = unit_correlation(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
signal = 'SUA';
pca_dims = 10;
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
if ~strcmp(signal,'SUA')&&~strcmp(signal,'MUA')&&~strcmp(signal,'denoised'), error('Signal to work on is not valid'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiking channel (we are only using LFPs from channels with spiking activity)
spk_ch = unique(trial_data(1).([array '_unit_guide'])(:,1));
lfp_ch = unique(trial_data(1).([array '_lfp_guide'])(:,1));
spk_ch = spk_ch(ismember(spk_ch,lfp_ch)); % Check that there is an LFP in each spiking channel

if strcmp(signal,'MUA')
    % Get multiunit activity
    for itrial = 1:length(trial_data)
        aux_array = zeros(size(trial_data(itrial).([array '_spikes']),1),length(spk_ch));
        for iCh = 1:length(spk_ch)
            isChan =  trial_data(1).([array '_unit_guide'])(:,1) == spk_ch(iCh);
            aux_array(:,iCh) = sum(trial_data(itrial).([array '_spikes'])(:,isChan),2);
        end
        trial_data(itrial).([array '_spikes']) = aux_array;
    end
end

% Smooth signals
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));

% Concatenate trials
unit_data = cell2mat({trial_data.([array '_spikes'])}');
lfp_data = cell2mat({trial_data.([array '_lfp'])}');

if strcmp(signal,'denoised')
    % Get denoised firing rates
    [coeff, score,~,~,~,mu] = pca(unit_data);
    unit_data = score(:,1:pca_dims)*coeff(:,1:pca_dims)'+mu;
end


% Get bands max freq and number of bands
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
bands_number = length(max_fq);


% Get the array map if possible to get far apart electrodes for the across_electrode analysis (it will avoid surrounding channels)
try 
    load([trial_data(1).monkey '_' array '_ArrayMap.mat']);
catch 
    arr_map = reshape(1:100,10,10);
end


% Compute the correlation
% freqs = unique(lfp_guide(:,3));
within_elect = nan(length(spk_ch),bands_number);
across_elect = nan(length(spk_ch),bands_number);
for iCh = 1:length(spk_ch)
    % Get a random channel for across_elect excluding surrounding channels
    far_ch = arr_map((conv2(double(arr_map==spk_ch(iCh)),ones(3),'same')) == 0);
    far_ch = spk_ch(ismember(spk_ch,far_ch));
    random_idx = randi(length(far_ch));
    random_ch = far_ch(random_idx);
    
    isChan_within = trial_data(1).([array '_lfp_guide'])(:,1) == spk_ch(iCh);
    isChan_across = trial_data(1).([array '_lfp_guide'])(:,1) == random_ch;
    for iBand = 1:bands_number
        isBand = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(iBand);
        lfp_idx = isChan_within & isBand;
        within_elect(iCh,iBand) = corr(unit_data(:,iCh),lfp_data(:,lfp_idx));
        lfp_idx = isChan_across & isBand;
        across_elect(iCh,iBand) = corr(unit_data(:,iCh),lfp_data(:,lfp_idx));
    end
end

within_elect = abs(within_elect); 
across_elect = abs(across_elect);


if doPlot
    % If there are names for the frequency bands
    if exist('bands_name.mat','file') == 2
        load('bands_name.mat');
    else % Else, assign ordered numbers
        bands_name = num2cell(1:bands_number);
    end
    bands_name{bands_number+1} = 'Lat Var';

    % Plot
    x = reshape([within_elect; across_elect],[],bands_number*2);
    x(any(isnan(x), 2), :) = [];
    pos1 = 1:3:bands_number*3; pos2 = 2:3:bands_number*3+1;
    positions = repmat(sort([pos1 pos2]),size(x,1),1);
    label_positions = (pos1+pos2)./2; 

    c = parula(bands_number*2); c(2:2:bands_number*2,:) = 0;
    for n = 1:bands_number*2
        hold on; scatter(positions(:,n),x(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
        hold on; p(n) = scatter(positions(1,n)+0.5,median(x(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
        hold on; errorbar(positions(1,n)+0.5,median(x(:,n)),std(x(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
    end
    legend([p(1) p(2)],{'Within electrode','Across electrode'})
    set(gca,'xtick',label_positions,'xticklabel',bands_name) 
    ylabel('Pearson correlation coefficient');
    title(trial_data(1).monkey);

end


end
