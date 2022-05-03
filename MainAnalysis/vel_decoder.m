%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [total_result] = vel_decoder(trial_data,params)
%
%   Will generate velocity decoders based on individual frequency bands and
%   return the VAF of the prediction together with the actual and predicted
%   velocities
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array         : (char) array location (e.g. 'M1')
%     .bins_to_past  : (int) Number of bins of history for the model
%     .folds         : (int) How many decoders to compute
%     .doPlot        : (logical) Whether to visualize result
%     .eval          : (cell) Indicate evaluation method (VAF, R-square, or both)
%     .pca_dims      : (int) Number of PCA dimensions
%
% OUTPUTS:
%   total_result: (struct) result of all the decoderes containing the
%   following fields. It has as many entries as decoders computed ([number of LFP bands +1]*folds)
%     total_result.actual_vel    : (array) Tx2 containing the x and y actual velocities
%     total_result.actual_vel    : (array) Tx2 containing the x and y predicted velocities
%     total_result.vaf(optional) : (array) 1x2 containing VAF for x and y velocities
%     total_result.r2(optional)  : (array) 1x2 containing r2 for x and y velocities
%     total_result.band          : (int) Number indicating the frequency band for each decoder
%
% Written by Cecilia Gallego-Carracedo. Updated April 2022.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [total_result] = vel_decoder(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
bins_to_past = 5;
folds = 20;
eval = {'r2'};
doPlot = true;
pca_dims = 8;
fix_traintest = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get spiking channel (we are only using LFPs from channels with spiking activity)
spk_ch = unique(trial_data(1).([array '_unit_guide'])(:,1));
lfp_ch = unique(trial_data(1).([array '_lfp_guide'])(:,1));
spk_ch = spk_ch(ismember(spk_ch,lfp_ch)); % Check that there is an LFP in each spiking channel

% Smooth signals
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));               
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));               

% Do PCA to obtain the latent dynamics
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));

% Decoders based on LFP at single frquencies bands
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
bands_number = length(max_fq);

single_result = cell(folds,bands_number+1);
for iband = 1:bands_number
    disp(['Band ' num2str(iband) ' out of ' num2str(bands_number)])
    
    % Prepare data for decoding
    single_band_td = trial_data;
    lfp_idx = find(trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(iband));
    lfp_idx = lfp_idx(spk_ch);
    for itrial = 1:length(trial_data)
        one_band_lfp = trial_data(itrial).([array '_lfp'])(:,lfp_idx);
        single_band_td(itrial).([array '_lfp']) = one_band_lfp;
    end
    
    % Compute cross-validated decoder
    for ifold = 1:folds
        %disp(['Fold ' num2str(n) ' out of ' num2str(folds)])
        single_result{ifold,iband} = ComputeVelModel(single_band_td,[array '_lfp'],struct('bins_to_past',bins_to_past,'eval',{eval},'fix_traintest',fix_traintest));
        single_result{ifold,iband}.band = iband;
    end
end
% Decoders based on Latent varibles
disp('Latent variables'); 
for ifold = 1:folds
    single_result{ifold,iband+1} = ComputeVelModel(trial_data,[array '_pca'],struct('bins_to_past',bins_to_past,'eval',{eval},'fix_traintest',fix_traintest));
    single_result{ifold,iband+1}.band = iband+1;
end

% Combine all the decoder's results in a single structure
total_result = [single_result{:}];

% Plot the data if required
if doPlot    
    % If there are names for the frequency bands
    if exist('bands_name.mat','file') == 2
        load('bands_name.mat');
    else % Else, assign ordered numbers
        bands_name = num2cell(1:bands_number);
    end
    bands_name{bands_number+1} = 'Lat Var';
   
    bands = [total_result.band];
    if length(eval) == 2 % Two evaluation metrics
        vaf = nan(folds*2,length(bands_name)); r2 = nan(folds*2,length(bands_name));
        for iband = 1:length(bands_name)
            idx = bands == iband;
            vaf(:,iband) = [total_result(idx).vaf]';
            r2(:,iband) = [total_result(idx).r2]';
        end

        figure
        c = parula(bands_number); c(bands_number+1,:) = [0 0 0];
        
        subplot(1,2,1)
        for n = 1:bands_number+1
            hold on; scatter(ones(size(vaf,1),1)*n,vaf(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
            hold on; scatter(n+0.25,median(vaf(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
            hold on; errorbar(n+0.25,median(vaf(:,n)),std(vaf(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
        end        
        set(gca,'xtick',1:bands_number+1,'xticklabel',bands_name)         
        xlim([0.5 bands_number+1.5]);
        ylabel('Decoder accuracy (vaf)');
        set(gca,'TickDir','out'); box off


        subplot(1,2,2)
        for n = 1:bands_number+1
            hold on; scatter(ones(size(r2,1),1)*n,r2(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
            hold on; scatter(n+0.25,median(r2(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
            hold on; errorbar(n+0.25,median(r2(:,n)),std(r2(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
        end        
        set(gca,'xtick',1:bands_number+1,'xticklabel',bands_name)         
        xlim([0.5 bands_number+1.5]); ylim([0 1])
        ylabel('Decoder accuracy (r2)');
        set(gca,'TickDir','out'); box off

        
    else
        ev = nan(folds*2,length(bands_name));
        for itrial = 1:length(unique(bands))
            idx = bands == itrial;
            ev(:,itrial) = [total_result(idx).(eval{1})]';
        end
        
        figure
        c = parula(bands_number); c(bands_number+1,:) = [0 0 0];
        for n = 1:bands_number+1
            hold on; scatter(ones(size(ev,1),1)*n,ev(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
            hold on; scatter(n+0.25,median(ev(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
            hold on; errorbar(n+0.25,median(ev(:,n)),std(ev(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
        end        
        set(gca,'xtick',1:bands_number+1,'xticklabel',bands_name)         
        xlim([0.5 bands_number+1.5]); ylim([0 1])
        ylabel(['Decoder accuracy (' eval{1} ')']);
        set(gca,'TickDir','out'); box off
    end
end


end