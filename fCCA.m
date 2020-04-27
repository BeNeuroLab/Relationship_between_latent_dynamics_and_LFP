%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [cca_coef,cca_bcoef,p] = fCCA(trial_data,params)
%
%   Will perform Canonical Correlation Analysis between PCA data and LFP at
%   each band individually
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   params     : struct with parameters
%     .array      : (char) array location (e.g. 'M1')
%     .pca_dims   : (int) Number of PCA dimensions
%     .boot_iter  : (int) Number of bootstraping iterations
%     .boot_type  : (char)
%           target_shift -> Rearrange target order
%           time_shuff -> Shuffle time bins in each trial
%           TME -> Tensor maximum entropy preserving Neurons covariance
%     .doPlot     : (logical) Whether to visualize result
%
% Written by Cecilia Gallego. Updated April 2020.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cca_coef,cca_bcoef,p] = fCCA(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default parameters
array = {};
pca_dims   = 8;
boot_iter = 5;
boot_type = 'TME';
doPlot = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_data = check_td_quality(trial_data);
if ~isempty(params), assignParams(who,params); end
if isempty(array), error('Need to provide a working array'); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Smooth spikes
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));

% Do PCA
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));

% Prepare data and perform CCA
channels = unique(trial_data(1).([array '_lfp_guide'])(:,1));
max_fq = unique(trial_data(1).([array '_lfp_guide'])(:,3));
min_fq = unique(trial_data(1).([array '_lfp_guide'])(:,2));
bands_name = cell(1,length(max_fq));
for n = 1:length(max_fq)
    if max_fq(n)== 0 && min_fq(n) == 0
        bands_name{n} = 'LMP';
    else
        bands_name{n} = [num2str(min_fq(n)) '-' num2str(max_fq(n))];
    end
end

lfp_data = cell2mat({trial_data.([array '_lfp'])}');
pca_data = cell2mat({trial_data.([array '_pca'])}');


% Compute CCA
cca_coef = cell(length(channels),length(max_fq)); 
U = cell(size(cca_coef)); V = cell(size(U));
A = cell(size(U)); B = cell(size(U));
for band = 1:length(max_fq)
    for ch = 1:length(channels)
        isFreq = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band);
        isChan = trial_data(1).([array '_lfp_guide'])(:,1) == channels(ch);
        idx = isFreq & isChan;
        indiv_lfp = lfp_data(:,idx);
        [A{ch,band},B{ch,band},cca_coef{ch,band},U{ch,band},V{ch,band}] = canoncorr(pca_data,indiv_lfp); 
    end
end
cca_coef = cell2mat(cca_coef);



% Boostrapping
cca_bcoef = zeros(length(channels),length(max_fq),boot_iter);
switch boot_type
    
    
    case 'target_shift'
        for i = 1:boot_iter
            disp(['Iteration number: ' num2str(i) ' out of ' num2str(boot_iter)]);
            idx = randperm(length(trial_data));
            boot_pca_data = cell2mat({trial_data(idx).([array '_pca'])}');

            temp_cca_coef = zeros(length(channels),length(max_fq)); 
            for band = 1:length(max_fq)
                for ch = 1:length(channels)
                    isFreq = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band);
                    isChan = trial_data(1).([array '_lfp_guide'])(:,1) == channels(ch);
                    idx = isFreq & isChan;
                    indiv_lfp = lfp_data(:,idx);
                    [~,~,temp_cca_coef(ch,band)] = canoncorr(boot_pca_data,indiv_lfp); 
                end
            end
            cca_bcoef(:,:,i) = temp_cca_coef;
        end
        cca_bcoef = mean(cca_bcoef,3);

        
        
    case 'time_shuff'
        for i = 1:boot_iter
            disp(['Iteration number: ' num2str(i) ' out of ' num2str(boot_iter)]);
            
            boot_trial_data = trial_data;
            for n = 1:length(boot_trial_data)
                single_trial = trial_data(n).([array '_spikes']);
                shuf_trial = single_trial(randperm(size(single_trial,1)),:);
                boot_trial_data(n).([array '_spikes']) = shuf_trial;
            end
            
            boot_trial_data = smoothSignals(boot_trial_data,struct('signals',[array '_spikes'],'width',0.05));
            [boot_trial_data,~] = dimReduce(boot_trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));
            boot_pca_data = cell2mat({boot_trial_data.([array '_pca'])}');


            temp_cca_coef = zeros(length(channels),length(max_fq)); 
            for band = 1:length(max_fq)
                for ch = 1:length(channels)
                    isFreq = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band);
                    isChan = trial_data(1).([array '_lfp_guide'])(:,1) == channels(ch);
                    idx = isFreq & isChan;
                    indiv_lfp = lfp_data(:,idx);
                    [~,~,temp_cca_coef(ch,band)] = canoncorr(boot_pca_data,indiv_lfp); 
                end
            end
            cca_bcoef(:,:,i) = temp_cca_coef;
        end
        cca_bcoef = mean(cca_bcoef,3);
        
        
        
    case 'TME'
        dataTensor = cat(3,trial_data.([array '_spikes']));
        for i = 1:boot_iter
            disp(['Iteration number: ' num2str(i) ' out of ' num2str(boot_iter)]);

            surrData = computeTME(dataTensor,'surrogate-TN');

            boot_trial_data = trial_data;
            for n = 1:length(boot_trial_data)
                boot_trial_data(n).([array '_spikes']) = surrData(:,:,n);
            end

            boot_trial_data = smoothSignals(boot_trial_data,struct('signals',[array '_spikes'],'width',0.05));
            [boot_trial_data,~] = dimReduce(boot_trial_data, struct('signals',[array '_spikes'],'num_dims',pca_dims));
            boot_pca_data = cell2mat({boot_trial_data.([array '_pca'])}');


            temp_cca_coef = zeros(length(channels),length(max_fq)); 
            for band = 1:length(max_fq)
                for ch = 1:length(channels)
                    isFreq = trial_data(1).([array '_lfp_guide'])(:,3) == max_fq(band);
                    isChan = trial_data(1).([array '_lfp_guide'])(:,1) == channels(ch);
                    idx = isFreq & isChan;
                    indiv_lfp = lfp_data(:,idx);
                    [~,~,temp_cca_coef(ch,band)] = canoncorr(boot_pca_data,indiv_lfp); 
                end
            end
            cca_bcoef(:,:,i) = temp_cca_coef;
        end
        cca_bcoef = mean(cca_bcoef,3);
    otherwise
        disp('Skipping bootstrapping')
end

p = []; for n = 1:9; p(end+1) = ranksum(cca_bcoef(:,n),cca_coef(:,n)); end



% Final plot
if doPlot
    x = reshape([cca_coef; cca_bcoef],[],18);

    pos1 = 1:3:27;
    pos2 = 2:3:28;
    positions = sort([pos1 pos2]);
    label_positions = (pos1+pos2)./2; 

    figure
    h = boxplot(x, 'positions', positions, 'color','k','symbol','k+');  
    set(h,'linewidth',1) 

    set(gca,'xtick',label_positions) 
    set(gca,'xticklabel',bands_name) 

    c = flipud(parula(18));
    c(1:2:18,:) = 0;
    h = findobj(gca,'Tag','Box'); 
    for j=1:length(h) 
    pa(j) = patch(get(h(j),'XData'),get(h(j),'YData'),c(j,:),'FaceAlpha',.5); 
    end 

%     legend([pa(2) pa(3)],{'Original CCA','Surrugate CCA'})

    ylabel('CCA coefficient');
    title([trial_data(1).monkey trial_data(1).date_time]);
    set(gca,'TickDir','out'); box off
    disp(p)
end



end