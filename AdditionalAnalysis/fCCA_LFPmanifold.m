%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cca_coef,cca_coef_surr,p] = fCCA_LFPmanifold(trial_data,params)
%
%   Will perform Canonical Correlation Analysis between actual latent 
%   dynamics and LFP at each frequency band individually. 
%   Then, will perform CCA between surrogate latent dynamics and actual 
%   LFP at each frequency band individually.
%   
%   
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array           : (char) array location (e.g. 'M1')
%     .pca_dims        : (int) Number of PCA dimensions
%     .surrogate_iter  : (int) Number of bootstraping iterations
%     .doPlot          : (logical) Whether to visualize result
%
% OUTPUTS:
%   cca_coef      : (matrix) contains correlation coefficients between 
%                       actual latent dynamics and LFP, row -> channel and 
%                       column -> frequency band
%   cca_coef_surr : (matrix) contains correlation coefficients between
%                       surrogate latent dynamics and LFP, row -> channel  
%                       and column -> frequency band
%   p             : (vector) contains P-value of wilcoxon ranksum test 
%                       between actual and surrogate correlation coefficients
%
% Written by Cecilia Gallego-Carracedo. Updated September 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cca_coef,cca_coef_surr,p] = fCCA_LFPmanifold(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
pca_dims = 8;
lfp_dims = 10;
cca_dims = 4;
surrogate_iter = 3;
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
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));

% Perform PCA to obtain the latent dynamics
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));

% Prepare data for CCA
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));

pca_data = cell2mat({trial_data.([array '_pca'])}');
lfp_data = cell2mat({trial_data.([array '_lfp'])}');

% Compute CCA
cca_coef = nan(cca_dims,length(max_fq)); 
for band = 1:length(max_fq)
    lfp_idx = find(trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band));
    lfp_idx = lfp_idx(spk_ch);
    indiv_band_lfp = lfp_data(:,lfp_idx);

    % Apply PCA to single band LFP
    [~,lfp_score] = pca(indiv_band_lfp);  
    lfp_score = lfp_score(:,1:lfp_dims);
    
    % Compute CCA
    [~,~,r,~,~] = canoncorr(pca_data,lfp_score); 
    cca_coef(:,band) = r(1:cca_dims);
end

% CONTINUE HERE

% Surrogate data analysis
cca_coef_surr = zeros(cca_dims,length(max_fq),surrogate_iter);

if surrogate_iter > 0   
    
    % Generate surrogate data
    dataTensor = cat(3,trial_data.([array '_spikes']));
    for iter = 1:surrogate_iter
        disp(['Iteration for surrogate number: ' num2str(iter) ' out of ' num2str(surrogate_iter)]);

        surrData = computeTME(dataTensor,'surrogate-N');
        
        % Order surrogate data into trial_data
        surr_trial_data = trial_data;
        for n = 1:length(surr_trial_data)
            surr_trial_data(n).([array '_spikes']) = surrData(:,:,n);
        end
        
        % Smooth signals
        surr_trial_data = sqrtTransform(surr_trial_data,[array '_spikes']);
        surr_trial_data = smoothSignals(surr_trial_data,struct('signals',[array '_spikes'],'width',0.05));
        
        % Compute surrogate latent dynamics
        [surr_trial_data,~] = dimReduce(surr_trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));
        surr_pca_data = cell2mat({surr_trial_data.([array '_pca'])}');

        % CCA of surrogate latent dynamcs and actual LFP
        temp_cca_coef = nan(cca_dims,length(max_fq)); 
        for band = 1:length(max_fq)
            lfp_idx = find(trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band));
            lfp_idx = lfp_idx(spk_ch);
            indiv_band_lfp = lfp_data(:,lfp_idx);
        
            % Apply PCA to single band LFP
            [~,lfp_score] = pca(indiv_band_lfp);  
            lfp_score = lfp_score(:,1:lfp_dims);
            
            % Compute CCA
            [~,~,r,~,~] = canoncorr(surr_pca_data,lfp_score); 
            temp_cca_coef(:,band) = r(1:cca_dims);
        end

        cca_coef_surr(:,:,iter) = temp_cca_coef;
    end
    % Actual surrogate result is the mean of all surrogate analysis
    cca_coef_surr = mean(cca_coef_surr,3);

    % P-value of wilcoxon ranksum test between actual and surrogate CCA coefficients
    p = []; for n = 1:size(cca_coef,2); p(end+1) = ranksum(cca_coef_surr(:,n),cca_coef(:,n)); end
else
    disp('Skipping bootstrapping')
    cca_coef_surr = nan(size(cca_coef)); p = [];
end



% Final plot
if doPlot
    data = reshape([cca_coef; cca_coef_surr],[],size(cca_coef,2)*2);

    pos1 = 1:3:size(cca_coef,2)*3;
    pos2 = 2:3:size(cca_coef,2)*3+1;
    positions = repmat(sort([pos1 pos2]),size(data,1),1);
    label_positions = (pos1+pos2)./2; 

    figure
    c = parula(size(cca_coef,2)*2); c(2:2:end,:) = 0;
    for n = 1:size(cca_coef,2)*2
        hold on; scatter(positions(:,n),data(:,n),20,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:),'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
        hold on; scatter(positions(1,n)+0.5,median(data(:,n)),50,'MarkerEdgeColor',c(n,:),'MarkerFaceColor',c(n,:));
        hold on; errorbar(positions(1,n)+0.5,median(data(:,n)),std(data(:,n)),'LineStyle','none','color',c(n,:),'linewidth',2);
    end

    set(gca,'xtick',label_positions) 
    lfp_guide = trial_data(1).([array '_lfp_guide']);
    bands_name = num2cell((lfp_guide(1:band,2)+lfp_guide(1:band,3))/2);
    set(gca,'xticklabel',bands_name) 

    ylabel('CCA coefficient'); xlabel('Mid-band frequency (Hz)');
    title([trial_data(1).monkey,'  ',trial_data(1).date_time])
end



end