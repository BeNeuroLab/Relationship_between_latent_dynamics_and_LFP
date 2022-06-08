%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [variance] = lfp_variance(trial_data,params)
%
%   Will compute the variance of each normalized LFP frequency band 
%   activity during each trial. Then, it will plot the variance 
%   distribution of each freq band as a violin plot 
%
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array           : (char) array location (e.g. 'M1')
%     .doPlot          : (logical) Whether to visualize result
%
% OUTPUTS:
%   variance      : (matrix) #samples x #bands array where each column 
%      corresponds to the variance values that the corresponding freq band 
%      has for a all trials and channels
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variance] = lfp_variance(trial_data,params)
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

% LFP bands
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
bands_number = length(max_fq);

% Compute variance of each frequency band during individual trials 
variance = nan(length(trial_data)*length(spk_ch),bands_number);
for iBand = 1:bands_number
    count = 0;
    isFreq = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(iBand);
    for iCh = 1:length(spk_ch)
        isChan = trial_data(1).([array '_lfp_guide'])(:,1) == spk_ch(iCh);
        lfp_idx = isFreq & isChan;
        for iTrial = 1:length(trial_data)
            count = count +1;
            lfp_data = trial_data(iTrial).([array '_lfp'])(:,lfp_idx);
            norm_data = (lfp_data - min(lfp_data))./(max(lfp_data) - min(lfp_data));
            variance(count,iBand) = var(norm_data);
        end
    end
end

if doPlot
    
    figure
    col = parula(bands_number);
    % If there are names for the frequency bands
    if exist('bands_name.mat','file') == 2
        load('bands_name.mat');
    else % Else, assign ordered numbers
        bands_name = num2cell(1:bands_number);
    end

    for iBand = 1:bands_number
        
        % Violin plot      
        [N,~] = histcounts(variance(:,iBand),-0.005:0.01:1.005);
        idx = find(N);
        dis = smooth(N(idx)); dis = dis./max(dis); 
        edge_label = 0:0.01:1; newlabel = edge_label(idx);
        newlabel = cat(2,fliplr(newlabel),newlabel); 
        dis = cat(1,flipud(-dis),dis);
        hold on; h = fill(dis+(3*iBand),newlabel,'r');
        set(h,'FaceColor',col(iBand,:));
        m = mean(variance(:,iBand));
        [~,pos] = min(abs(newlabel-m)); x1 = dis(pos); x2 = -x1; 
        hold on; line([x1+(3*iBand) x2+(3*iBand)],[m m],'color','k','linewidth',1.5)
    end
    
    ylim([0,0.5]); ylabel('Variance')
    set(gca,'xtick',3:3:bands_number*3,'xticklabel',bands_name)
    set(gca,'TickDir','out'); box off;
    xlim([1,bands_number*3+2]);
end

end