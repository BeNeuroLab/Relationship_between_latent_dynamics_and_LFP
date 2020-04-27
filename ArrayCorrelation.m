%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cca_freq,cca_boot] = ArrayCorrelation(trial_data,params)
%
%   Will compute the Pearson's correlation between the multiunit spiking
%   form a channel and the LFP signal at different frequencies form the
%   same channel
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) array location (e.g. 'M1')
%     .doPlot     : (logical) Whether to visualize result
%
% Written by Cecilia Gallego. Updated April 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [correlation,boot_correlation] = ArrayCorrelation(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get multiunit activity
spikes_guide = trial_data(1).([array '_unit_guide']);
spikes_ch = unique(spikes_guide(:,1));
for ch = 1:length(spikes_ch)
    isChan =  spikes_guide(:,1) == spikes_ch(ch);
    for trial = 1:length(trial_data)
        trial_data(trial).([array '_MUA'])(:,ch) = sum(trial_data(trial).([array '_spikes'])(:,isChan),2);
    end
end

% Smooth signals and concatenate them across time
trial_data = smoothSignals(trial_data, struct('signals',[array '_MUA'],'width',0.05));
trial_data = smoothSignals(trial_data, struct('signals',[array '_lfp'],'width',0.05));
spike_data = cell2mat({trial_data.([array '_MUA'])}');
lfp_data = cell2mat({trial_data.([array '_lfp'])}');


% % Correlation of each lfp signal with the spikes in that same channel
lfp_guide = trial_data(1).([array '_lfp_guide']);
lfp_ch = unique(lfp_guide(:,1));


% Compute the correlation
freqs = unique(lfp_guide(:,3));
correlation = zeros(length(spikes_ch),length(freqs));
remove_idx = [];
for ch = 1:length(spikes_ch)
    % Check if the spiking channels have lfp signals
    if ~ismember(spikes_ch(ch),lfp_ch)
        remove_idx(end+1) = ch;
        continue
    end
    isChan = lfp_guide(:,1) == spikes_ch(ch);
    for b = 1:length(freqs)
        isBand = lfp_guide(:,3) == freqs(b);
        lfp_idx = isChan & isBand;
        correlation(ch,b) = corr(spike_data(:,ch),lfp_data(:,lfp_idx));
    end
end
% Eliminate the channels without lfp signal
if ~isempty(remove_idx)
    correlation(remove_idx,:) = [];
    spikes_ch(remove_idx) = [];
end


% Get the array map if possible for the bootstrapping (it will avoid surrounding channels)
try 
    load([trial_data(1).monkey '_' array '_ArrayMap.mat']);
catch 
    arr_map = reshape(1:100,10,10);
end


boot_correlation = zeros(size(correlation));
for ch = 1:length(spikes_ch)
    % Get a random index from bootstraping excluding surrounding channels
    boot_not_ch = arr_map((conv2(double(arr_map==spikes_ch(ch)),ones(3),'same')) == 1);
    boot_ch = spikes_ch(~ismember(spikes_ch,boot_not_ch));
    boot_idx = randperm(length(boot_ch)); boot_idx = boot_idx(1);
    
    % Compare the channel MUA with the LFP of a random channel (boot_idx)
    isChan = lfp_guide(:,1) == spikes_ch(boot_idx);
    for b = 1:length(freqs)
        isBand = lfp_guide(:,3) == freqs(b);
        lfp_idx = isChan & isBand;
        boot_correlation(ch,b) = corr(spike_data(:,ch),lfp_data(:,lfp_idx));
    end
end

correlation = abs(correlation); boot_correlation = abs(boot_correlation);

if doPlot
    % Final plot
    x = reshape([correlation; boot_correlation],[],18);

    pos1 = 1:3:27;
    pos2 = 2:3:28;
    positions = sort([pos1 pos2]);
    label_positions = (pos1+pos2)./2; 

    figure
    h = boxplot(x, 'positions', positions, 'color','k','symbol','k+');  
    set(h,'linewidth',1) 

    set(gca,'xtick',label_positions) 
    load('bands_name.mat');
    set(gca,'xticklabel',bands_name) 

    aux_corr = flipud(parula(18));
    aux_corr(1:2:18,:) = 0;
    h = findobj(gca,'Tag','Box'); 
    for j=1:length(h) 
    p(j) = patch(get(h(j),'XData'),get(h(j),'YData'),aux_corr(j,:),'FaceAlpha',.5); 
    end 

    legend([p(2) p(3)],{'Same-channel correlation','Random-channel correlation'})

    ylabel('Pearson correlation coefficient');
end


end
