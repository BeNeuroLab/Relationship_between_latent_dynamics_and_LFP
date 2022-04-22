%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [corr_data] = lfp_lfp_corr(trial_data,params)
%
%   Will compute the correlation between LFP frequency bands for each
%   channel individually
%
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array           : (char) array location (e.g. 'M1')
%     .doPlot          : (logical) Whether to visualize result
%
% OUTPUTS:
%   corr_data  : (matrix) #bands x #bands x #channels array where each
%      entry is the correlation between band1 (1st dim) and band2 (2 dim) 
%      for a given channel (3r dim)
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [corr_data] = lfp_lfp_corr(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array location'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiking channel (we are only using LFPs from channels with spiking activity)
spk_ch = unique(trial_data(1).([array '_unit_guide'])(:,1));
lfp_ch = unique(trial_data(1).([array '_lfp_guide'])(:,1));
spk_ch = spk_ch(ismember(spk_ch,lfp_ch)); % Check that there is an LFP in each spiking channel

% Smooth signals
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));

% LFP bands
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
bands_number = length(max_fq);
lfp_data = cell2mat({trial_data.([array '_lfp'])}');

% Reorder data
data = nan(size(lfp_data,1),length(spk_ch),bands_number);
for iBand = 1:bands_number
    band_idx = find(trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(iBand));
    data(:,:,iBand) = lfp_data(:,band_idx(spk_ch));
end
% Compute correlation
corr_data = nan(bands_number,bands_number,length(spk_ch));
for iCh = 1:size(data,2)
    corr_data(:,:,iCh) = corr(squeeze(data(:,iCh,:)));
end
corr_data = abs(corr_data);


if doPlot
    figure

    % If there are names for the frequency bands
    if exist('bands_name.mat','file') == 2
        load('bands_name.mat');
    else % Else, assign ordered numbers
        bands_name = num2cell(1:bands_number);
    end

    se = std(corr_data,0,3)./sqrt(size(corr_data,3));
    m = mean(corr_data,3);    

    matrix = nan(bands_number,bands_number);
    matrix(triu(ones(bands_number),1)==1) = m(triu(ones(bands_number),1)==1);
    matrix(tril(ones(bands_number),-1)==1) = se(tril(ones(bands_number),-1)==1);
    imagesc(matrix); colormap 'hot'; caxis([0,0.5]); colorbar
    set(gca,'yticklabel',bands_name(1:bands_number)); 
    set(gca,'xtick',1:bands_number,'xticklabel',bands_name(1:bands_number)); 
    set(gca,'TickDir','out'); box off
    xtickangle(45);

end

end