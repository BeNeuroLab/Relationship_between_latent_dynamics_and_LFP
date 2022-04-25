close all; clear; clc;
 
root_path = 'C:\Users\Cecilia\Documents\BeNeuroLab\Project_LFPvsMLatents';
data_path = fullfile(root_path,'Data');
addpath(genpath(data_path));
% Add trial_data repo to path
addpath(genpath(fullfile(root_path,"TrialData")));

% Import trial_data structure
load(fullfile(data_path,'filenames.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 2
load(fullfile(data_path,'raw_data_fig2')); % Reduced raw data form Chewie_CO_CS_BL_10212016_001
t = table2array(kin(:,1));
goCueTime = table2array(trials(:,6));
bin_size = 0.01; kernel_SD = 0.05;
% Behaviour
figure
subplot(4,1,1)
hold on; plot(t,table2array(kin(:,2))); hold on; plot(t,table2array(kin(:,3)))
for cue = goCueTime'
    hold on; line([cue cue],[-40 40],'Color','black');
    hold on; line([cue+0.24 cue+0.24],[-40 40],'Color','black');
    hold on; line([cue+0.7 cue+0.7],[-40 40],'Color','black');
end
set(gca,'TickDir','out'); box off;
% Firing rate
subplot(4,1,2)
fr=[]; t_plot = [];
for n = 1:3:size(spikes,1)-1
    fr(end+1,:) = sum(spikes(n:n+2,:));
    t_plot(end+1) = t(n);
end
fr = smooth_data(fr,bin_size,kernel_SD);
fr_plot = fr(:,mean(fr)>0.1); fr_plot = fr_plot./max(fr_plot);
surf(t_plot,1:size(fr_plot,2),fr_plot','EdgeColor', 'none'); colormap(hot); view(0,90);
for cue = goCueTime'
    hold on; line([cue cue],[75 80],'Color','black');
    hold on; line([cue+0.24 cue+0.24],[75 80],'Color','black');
    hold on; line([cue+0.7 cue+0.7],[75 80],'Color','black');
end
set(gca,'TickDir','out'); box off;
% Latent dynamics
subplot(4,1,3)
[~,pca_data] = pca(fr);
c = repmat([0:0.125:1]',1,3);
for n = 1:8
    norm_data = (pca_data(:,n)-min(pca_data(:,n)))./(max(pca_data(:,n))-min(pca_data(:,n)));
    hold on; plot(t_plot,norm_data+8-n,'color',c(n,:));
end
for cue = goCueTime'
    hold on; line([cue cue],[0 9],'Color','black');
    hold on; line([cue+0.24 cue+0.24],[0 9],'Color','black');
    hold on; line([cue+0.7 cue+0.7],[0 9],'Color','black');
end
set(gca,'TickDir','out'); box off;
% LFP
subplot(4,1,4)
c = parula(9);
samprate = 2000;
window_size = 0.05*samprate; % 50 ms windows
ch = 34;
t = table2array(lfp(:,1));   
ch_idx = lfp_guide(:,1) == ch;
for n = 1:9
    band_idx = lfp_guide(:,3) == lfp_guide(n,3);
    lfp_idx = ch_idx & band_idx;
    data = lfp_filt(:,lfp_idx);
    norm_data = (data-min(data))./(max(data)-min(data));
    hold on; plot(t,norm_data+n,'color',[0.5 0.5 0.5]);
    
    if lfp_guide(lfp_idx,2) == 0
        data = movmean(data,window_size);
    else
        temp_window_size = max(window_size,samprate/lfp_guide(lfp_idx,2));
        data = movmean(data.^2,temp_window_size);
    end
    data = smooth_data(data,bin_size,kernel_SD);
    norm_data = (data-min(data))./(max(data)-min(data));
    hold on; plot(t,norm_data+n,'color',c(n,:),'linewidth',1.5);
end
for cue = goCueTime'
    hold on; line([cue cue],[1 10],'Color','black');
    hold on; line([cue+0.24 cue+0.24],[1 10],'Color','black');
    hold on; line([cue+0.7 cue+0.7],[1 10],'Color','black');
end
set(gca,'TickDir','out'); box off;
xlabel('seconds')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3
file = 12;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = binTD(trial_data,3);
trial_data = trim_data(trial_data,'exec');
[cca_coef,cca_coef_surr,p,A,B,U,V] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));

trials = [2,5,7,18];
trial_len = 19; spacing = 5;
% Smooth signals
trial_data = sqrtTransform(trial_data,[array '_spikes']);
trial_data = smoothSignals(trial_data,struct('signals',[array '_spikes'],'width',0.05));
trial_data = smoothSignals(trial_data,struct('signals',[array '_lfp'],'width',0.05));
[trial_data,~] = dimReduce(trial_data, struct('signals',[array '_spikes'],'num_dims',10));
lfp = {trial_data.([array '_lfp'])};
pca = {trial_data.([array '_pca'])};
% Only Latent dynamics
figure
subplot(1,3,1)
pos = 1;
for t = 1:length(trials)
    for n = 1:10
        hold on; plot(pos:pos+trial_len-1,pca{trials(t)}(:,n)+45-5*(n-1),'color',[0 0 0]+(n-1)/10,'linewidth',1.5);
        name{10-(n-1)} = ['Lat var ' num2str(n)];
    end
    pos = pos+trial_len+spacing;
end
set(gca,'TickDir','out'); title('Neural modes before');
xticks([10,35,60,82]); xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4'});
yticks([0:5:45]); yticklabels(name); ylim([-5,50])
% Only LFPs
subplot(1,3,2)
c = parula(9);
ch = spk_ch(45);
ch_idx = trial_data(1).M1_lfp_guide(:,1) == ch; 
bands = unique(trial_data(1).M1_lfp_guide(:,3));
pos = 1;
for t = 1:length(trials)
    for b = 1:9
        band_idx = trial_data(1).M1_lfp_guide(:,3) == bands(b);
        lfp_idx = ch_idx & band_idx;
        data = lfp{trials(t)}(:,lfp_idx); data = (data-min(data))/max(data-min(data));
        hold on; plot(pos:pos+trial_len-1,data+3*(b-1),'color',c(b,:),'linewidth',1.5);
    end
    pos = pos+trial_len+spacing;
end
set(gca,'TickDir','out'); title('LFPs'); 
xticks([10,35,60,82]); xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4'});
load('bands_name.mat'); yticklabels(bands_name); yticks([0.5:3:24.5]); 
% LFP-Lantent dynamics aligned
subplot(1,3,3)
pos = 1; align = []; ch = 45;
for t = 1:length(trials)
    aux_pca = []; aux_lfp = [];    
    trial_idx = 1+(trial_len*(trials(t)-1)):1+(trial_len*(trials(t)-1))+trial_len-1;  
    for b = 1:9        
        if b == 1 ||  b == 5 ||  b == 6; data = V{ch,b}(trial_idx).*-1; % Signals inverted by CCA   
        else; data = V{ch,b}(trial_idx); end
        data = (data-min(data))/max(data-min(data));
        hold on; plot(pos:pos+trial_len-1,data+3*(b-1),'color',c(b,:),'linewidth',1.5);
        aux_lfp = cat(2,aux_lfp,V{ch,b}(trial_idx));
    end    
    for b = 1:9
        if b == 1 ||  b == 5 ||  b == 6; data = U{ch,b}(trial_idx).*-1;
        else; data = U{ch,b}(trial_idx); end
        data = (data-min(data))/max(data-min(data));
        hold on; plot(pos:pos+trial_len-1,data+3*(b-1),'color','k','linewidth',1.5);
        aux_pca = cat(2,aux_pca,U{ch,b}(trial_idx));
    end
    pos = pos+trial_len+spacing;
    aux_cor = corr(aux_pca,aux_lfp);
    align(t,1:9) = diag(aux_cor);    
end
set(gca,'TickDir','out'); yticks([0.5:3:24.5]); load('bands_name.mat'); yticklabels(bands_name)
xticks([10,35,60,82]); xticklabels({'Trial 1','Trial 2','Trial 3','Trial 4'});
title('Neural modes after');
targets = [trial_data.tgtDir];
disp(['Tragets: ' num2str(targets(trials))]);
disp(['Correlations: ' num2str(mean(align,1))]);



