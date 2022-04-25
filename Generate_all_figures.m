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
ch = 34;
t = table2array(lfp(:,1));   
ch_idx = lfp_guide(:,1) == ch;
for n = 1:9
    band_idx = lfp_guide(:,3) == lfp_guide(n,3);
    lfp_idx = ch_idx & band_idx;
    data = lfp_filt(:,lfp_idx);
    norm_data = (data-min(data))./(max(data)-min(data));
    hold on; plot(t,norm_data+n,'color',[0.5 0.5 0.5]);

    data = smooth_data(lfp_data(:,ch),bin_size,kernel_SD);
    norm_data = (data-min(data))./(max(data)-min(data));
    hold on; plot(t,norm_data+n,'color',c(n,:),'linewidth',1.5);
end
for cue = goCueTime'
    hold on; line([cue cue],[0 9],'Color','black');
    hold on; line([cue+0.24 cue+0.24],[0 9],'Color','black');
    hold on; line([cue+0.7 cue+0.7],[0 9],'Color','black');
end
set(gca,'TickDir','out'); box off;
xlabel('seconds')

