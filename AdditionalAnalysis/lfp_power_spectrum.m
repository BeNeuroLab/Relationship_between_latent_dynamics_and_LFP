%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [freq_cont,freq] = lfp_power_spectrum(cds_path,params)
%
%   Will compute the power spectrum of the original LFP signal (before
%   filtering to specific frequency bands)
%
% INPUTS:
%   cds_path   : (char) path for the original file (the CDS file)
%   params     : struct with parameters
%     .array           : (char) array location (e.g. 'M1')
%     .epoch           : (char) behavioural epoch of interest (only 'exec' or 'prep')
%     .max_freq        : (int) maximum frequency to show in plot
%     .doPlot          : (logical) Whether to visualize result
%
% OUTPUTS:
%   freq_cont : (matrix) #samples x #freqs array where each column 
%      corresponds to the frequency content values that the corresponding 
%      freq has for a all trials and channels
%   freq      : 1D array containing the freq value for each column of freq_cont
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [freq_cont,freq] = lfp_power_spectrum(cds_path,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
epoch = 'exec';
max_freq = 400;
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array location'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Open origial dataset
load(cds_path)

fs = 1/(cds.lfp.t(2)-cds.lfp.t(1));
lfp_time = cds.lfp.t;
lfp = cds.lfp;

% Find LFP channels of the array of interest
names = cds.lfp.Properties.VariableNames;
lfp_idx = startsWith(names,array);    
if sum(lfp_idx) == 0 
    error('No lfp data for given array')
end
   
timing = cds.trials;
timing = timing(~isnan(timing.goCueTime),:);

epoch_idx = [];
switch epoch
    case 'exec'
        for iTrial = 1:size(timing,1)
            epoch_idx = cat(1,epoch_idx,find(lfp_time>timing.goCueTime(iTrial)+0.1,1):find(lfp_time>timing.goCueTime(iTrial)+0.1,1)+0.5*fs);     
        end
    case 'prep'
        for iTrial = 1:size(timing,1)
            epoch_idx = cat(1,epoch_idx,find(lfp_time>timing.goCueTime(iTrial)-0.3,1):find(lfp_time>timing.goCueTime(iTrial)-0.3,1)+0.5*fs);        
        end
    otherwise
        error('Not valid behavioural epoch for this analysis')
end

% Clear cds to save memory
clear 'cds'


lfp_single_aray = table2array(lfp(:,lfp_idx==1));

% filter the line noise
[B,A]             = butter(2,( 60 + [-2 2] )/(fs/2),'stop');
lfp_single_aray        = filtfilt(B,A,lfp_single_aray);
[B,A]             = butter(2,( 120 + [-2 2] )/(fs/2),'stop');
lfp_single_aray        = filtfilt(B,A,lfp_single_aray);
[B,A]             = butter(2,( 240 + [-2 2] )/(fs/2),'stop');
lfp_single_aray        = filtfilt(B,A,lfp_single_aray);
[B,A]             = butter(2,( 360 + [-2 2] )/(fs/2),'stop');
lfp_single_aray        = filtfilt(B,A,lfp_single_aray);


freq_cont = nan(size(epoch_idx,1)*size(lfp_single_aray,2),ceil(size(epoch_idx,2)/2));
pos = 0;
warning off
for iTrial = 1:size(epoch_idx,1)        
    for iCh = 1:size(lfp_single_aray,2)
        pos = pos +1;
        data = lfp_single_aray(epoch_idx(iTrial,:),iCh);
        N = length(data); 
        xdft = fft(data); 
        xdft = xdft(1:(N/2+1)); 
        psdx = (1/(fs*N)) * abs(xdft).^2; 
        psdx(2:end-1) = 2*psdx(2:end-1); 
        freq_cont(pos,:) = psdx;
    end
end
warning on
freq = linspace(0,fs/2,size(freq_cont,2)); 


if doPlot

    figure
    stdshade(10*log10(freq_cont),0.5,[0,0,0],freq);
    grid on; box off; set(gca,'TickDir','out'); 
    xlabel('Frequency (Hz)'); xlim([0,max_freq])
    ylabel('Power/Frequency (dB/Hz)')
end


end