close all; clear; clc;
 
root_path = 'C:\Users\Cecilia\Documents\Projects\Project_LFPvsMLatents';
data_path = fullfile(root_path,'Data');
addpath(genpath(data_path));
% Add trial_data repo to path
addpath(genpath(fullfile(root_path,"TrialData")));

% Import filenames and bands_name
load(fullfile(data_path,'filenames.mat'));
load(fullfile(data_path,'bands_name.mat'));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 3a
file = 12;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[cca_coef,cca_coef_surr,p,A,B,U,V] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
trials = [2,5,7,18];
trial_len = 19; spacing = 5;
% Smooth signals
trial_data = sqrtTransform(trial_data,'M1_spikes');
trial_data = smoothSignals(trial_data,struct('signals','M1_spikes','width',0.05));
trial_data = smoothSignals(trial_data,struct('signals','M1_lfp','width',0.05));
[trial_data,~] = dimReduce(trial_data, struct('signals','M1_spikes','num_dims',10));
lfp = {trial_data.('M1_lfp')};
pca = {trial_data.('M1_pca')};
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
ch = 91;
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
% Figure 3b
figure
subplot(1,3,1)
file = 3;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',1,'doPlot',true));
subplot(1,3,2)
file = 12;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',1,'doPlot',true));
subplot(1,3,3)
file = 15;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',1,'doPlot',true));
% Figure 3c
figure
total_cca = {};
for file = 1:16
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data = trim_data(trial_data,'exec');
    [cca_coeff,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca{file} = cca_coeff;
end
data = cell2mat(total_cca');
edge_label = 0:0.02:1; c = parula(9);
for b = 1:9
    single_dis = histcounts(data(:,b),-0.01:0.02:1.01);
    idx = find(single_dis);
    dis = smooth(single_dis(idx)); dis = dis./max(dis); newlabel = edge_label(idx);
    newlabel = cat(2,fliplr(newlabel),newlabel); 
    dis = cat(1,flipud(-dis),dis);
    hold on; h(b) = fill(dis+(3*b),newlabel,'r');
    set(h(b),'FaceColor',c(b,:));
    single_data = data(:,b);
    m = mean(single_data); [~,pos] = min(abs(newlabel-m)); x1 = dis(pos); x = -x1; 
    hold on; line([x1+(3*b) x+(3*b)],[m m],'color','k','linewidth',1.5)
end
xlim([1,32]); ylim([0 1]); set(gca,'xtick',3:3:27,'xticklabel',bands_name)
ylabel('CCA coefficient'); set(gca,'TickDir','out'); box off;
title('All CCA M1')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4a
file = 12;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[decoder_result] = vel_decoder(trial_data,struct('array',filenames{file,2},'folds',1,'bins_to_past',3,'pca_dims',10,'doPlot',false,'fix_traintest',true));
figure
c = parula(9); c(10,:) = [0.5 0.5 0.5];
pos = 1; spacing = 5; trial_len = 16;
for t = 1:10
    idx = 1+(trial_len*(t-1));
    hold on; plot(pos:pos+trial_len-1,decoder_result(1).actual_vel(idx:idx+trial_len-1,1),'k','linewidth',2)
    for b = [5,9,10]
        hold on; plot(pos:pos+trial_len-1,decoder_result(b).predic_vel(idx:idx+trial_len-1,1),'color',c(b,:),'linewidth',1.5);
    end
    pos = pos+trial_len+spacing;
end
% Figure 4b
figure
total_r2 = {};
for file = 1:16
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data = trim_data(trial_data,'exec');
    [decoder_result] = vel_decoder(trial_data,struct('array',filenames{file,2},'folds',50,'bins_to_past',3,'pca_dims',10,'doPlot',false,'eval',{{'r2'}}));
    total_r2{file} = reshape([decoder_result.r2],[],10);
end
data = cell2mat(total_r2');
edge_label = 0:0.02:1; c = parula(9); c(10,:) = [0 0 0];
for b = 1:10
    single_dis = histcounts(data(:,b),-0.01:0.02:1.01);
    idx = find(single_dis);
    dis = smooth(single_dis(idx)); dis = dis./max(dis); newlabel = edge_label(idx);
    try newlabel = cat(2,edge_label(idx(1)-1),newlabel,edge_label(idx(end)+1)); dis = cat(1,0,dis,0); end
    newlabel = cat(2,fliplr(newlabel),newlabel); 
    dis = cat(1,flipud(-dis),dis);
    hold on; h(b) = fill(dis+(3*b),newlabel,'r');
    set(h(b),'FaceColor',c(b,:));
    single_data = data(:,b);
    m = mean(single_data); [~,pos] = min(abs(newlabel-m)); x1 = dis(pos); x = -x1; 
    if b == 10
        hold on; line([x1+(3*b) x+(3*b)],[m m],'color',[1 1 1],'linewidth',1.5)
    else
        hold on; line([x1+(3*b) x+(3*b)],[m m],'color','k','linewidth',1.5)
    end
end
xlim([1,32]); ylim([0 1]); set(gca,'xtick',3:3:30,'xticklabel',bands_name)
ylabel('Decoder Accuracy'); set(gca,'TickDir','out'); box off;
title('All decoders M1')
% Figure 4c
figure
c = parula(9); total_x = []; total_y = [];
for j = 1:16
    x = mean(total_cca{j});
    y = mean(total_r2{j}); y(:,10) = [];
    for b = 1:9
        if ismember(j,1:6)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',c(b,:),'MarkerFaceColor',c(b,:));
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',c(b,:),'LineWidth',2);
        elseif ismember(j,13:16)
            hold on; e(b) = scatter(x(b),y(b),'+','MarkerEdgeColor',c(b,:),'LineWidth',2);
        end
    end
    total_y = cat(1,total_y,y');
    total_x = cat(1,total_x,x');
end
mdl = fitlm(total_x,total_y);
hold on; h = plot(mdl);
delete(h(1)); delete(h(3)); delete(h(4)); set(h(2),'color','k')
legend([d(1) a(1) e(1)],{'Monkey M','Monkey CL','Monkey CR'});
ylabel('Decoder performance'); xlabel('CCA correlation coefficient');
set(gca,'TickDir','out'); box off; title('All M1');
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 5a
figure
subplot(1,3,1)
total_cca_exec = {}; total_cca_prep = {}; total_cca_rest = {};
for file = 1:16
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data_temp = trim_data(trial_data,'exec');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'prep');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_prep{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); hold on; stdshade(data_exec,0.5,[1,0,0])
data_prep = cell2mat(total_cca_prep'); hold on; stdshade(data_prep,0.5,[0,0,1])
data_rest = cell2mat(total_cca_rest'); hold on; stdshade(data_rest,0.5,[0,0,0])
set(gca,'xticklabel',bands_name); ylabel('CCA coefficient'); ylim([0,0.8])
set(gca,'TickDir','out'); box off; 
subplot(1,3,2)
for j = 1:16
    x = data_prep(j,:);
    y = data_exec(j,:);
    for b = 1:9
        if ismember(j,1:6)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        elseif ismember(j,13:16)
            hold on; e(b) = scatter(x(b),y(b),'+','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        end
    end
end
legend([d(1) a(1) e(1)],{'Monkey M','Monkey CL','Monkey CR'});
xlabel('Preparation'); ylabel('Execution'); xlim([0,0.8]); ylim([0,0.8])
set(gca,'TickDir','out'); box off; 
subplot(1,3,3)
for j = 1:16
    x = data_rest(j,:);
    y = data_exec(j,:);
    for b = 1:9
        if ismember(j,1:6)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        elseif ismember(j,13:16)
            hold on; e(b) = scatter(x(b),y(b),'+','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        end
    end
end
xlabel('Rest'); ylabel('Execution'); xlim([0,0.8]); ylim([0,0.8])
set(gca,'TickDir','out'); box off; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6b
figure
subplot(1,2,1)
file = 17;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'prep');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',1,'doPlot',true));
subplot(1,2,2)
file = 26;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'prep');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',1,'doPlot',true));
% Figure 6c&d
figure
total_cca_exec = {}; total_cca_prep = {}; total_cca_rest = {};
for file = 17:28
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data_temp = trim_data(trial_data,'exec');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'prep');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
    total_cca_prep{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); data_prep = cell2mat(total_cca_prep'); data_rest = cell2mat(total_cca_rest'); 
subplot(1,2,1)
for j = 1:12
    x = data_prep(j,:);
    y = data_exec(j,:);
    for b = 1:9
        if ismember(j,1:6)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        end
    end
end
legend([d(1) a(1)],{'Monkey M','Monkey CL'});
xlabel('Preparation'); ylabel('Execution'); xlim([0,0.8]); ylim([0,0.8])
set(gca,'TickDir','out'); box off; 
subplot(1,2,2)
for j = 1:12
    x = data_rest(j,:);
    y = data_exec(j,:);
    for b = 1:9
        if ismember(j,1:6)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        end
    end
end
xlabel('Rest'); ylabel('Execution'); xlim([0,0.8]); ylim([0,0.8])
set(gca,'TickDir','out'); box off; 
% 6e 
figure
total_clasif = {};
for file = 17:28
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data = trim_data(trial_data,'prep');
    [total_clasif{file}] = target_classifier(trial_data,struct('array',filenames{file,2},'folds',100,'pca_dims',15,'doPlot',false));
end
data = cell2mat(total_clasif');
edge_label = 0:0.02:1; c = parula(9); c(10,:) = [0 0 0];
for b = 1:10
    single_dis = histcounts(data(:,b),-0.01:0.02:1.01);
    idx = find(single_dis);
    dis = smooth(single_dis(idx)); dis = dis./max(dis); newlabel = edge_label(idx);
    try newlabel = cat(2,edge_label(idx(1)-1),newlabel,edge_label(idx(end)+1)); dis = cat(1,0,dis,0); end
    newlabel = cat(2,fliplr(newlabel),newlabel); 
    dis = cat(1,flipud(-dis),dis);
    hold on; h(b) = fill(dis+(3*b),newlabel,'r');
    set(h(b),'FaceColor',c(b,:));
    single_data = data(:,b);
    m = mean(single_data); [~,pos] = min(abs(newlabel-m)); x1 = dis(pos); x = -x1; 
    if b == 10
        hold on; line([x1+(3*b) x+(3*b)],[m m],'color',[1 1 1],'linewidth',1.5)
    else
        hold on; line([x1+(3*b) x+(3*b)],[m m],'color','k','linewidth',1.5)
    end
end
xlim([1,32]); ylim([0 1]); set(gca,'xtick',3:3:30,'xticklabel',bands_name)
ylabel('Accuracy'); set(gca,'TickDir','out'); box off;
title('All PMd Classifiers')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7b
figure
subplot(1,3,1)
file = 32;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'feed');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',8,'surrogate_iter',1,'doPlot',true));
subplot(1,3,2)
file = 35;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'feed');
[~,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',8,'surrogate_iter',1,'doPlot',true));
% Figure 7c
subplot(1,3,3)
total_cca_exec = {}; total_cca_rest = {};
for file = 29:36
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    trial_data_temp = trim_data(trial_data,'feed');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); data_rest = cell2mat(total_cca_rest'); 
for j = 1:8
    x = data_rest(j,:);
    y = data_exec(j,:);
    for b = 1:9
        if ismember(j,1:5)
            hold on; d(b) = scatter(x(b),y(b),'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        elseif ismember(j,7:12)
            hold on; a(b) = scatter(x(b),y(b),'square','MarkerEdgeColor',[0,0,0],'LineWidth',2);
        end
    end
end
xlabel('Rest'); ylabel('Execution'); xlim([0,0.7]); ylim([0,0.7])
set(gca,'TickDir','out'); box off; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 8a
figure
subplot(1,3,1)
file = 12;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'exec');
[~, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',10,'doPlot',true));
ylim([0,0.8]); set(gca,'TickDir','out'); box off; 
subplot(1,3,2)
file = 26;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'prep');
[~, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',15,'doPlot',true));
ylim([0,0.8]); set(gca,'TickDir','out'); box off; 
subplot(1,3,3)
file = 32;
trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
trial_data = trim_data(trial_data,'feed');
[~, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',8,'doPlot',true));
ylim([0,0.8]); set(gca,'TickDir','out'); box off; 
% Figure 8b
figure
col = parula(9);
for file = 1:36
    trial_data = loadTDfiles(filenames{file,1},{@getTDidx,{'result','R'}}); 
    if ismember(file,1:16)
        trial_data = trim_data(trial_data,'exec');
        [pears_coef, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',10,'doPlot',false));
        [cca_coef,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
        for band = 1:9
            x = median(cca_coef(:,band)); y = median(pears_coef(:,band));
            hold on; scatter(x,y,'MarkerEdgeColor',col(band,:),'MarkerFaceColor',col(band,:));
        end
    elseif ismember(file,17:28)
        trial_data = trim_data(trial_data,'prep');
        [pears_coef, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',15,'doPlot',false));
        [cca_coef,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',15,'surrogate_iter',0,'doPlot',false));
        for band = 1:9
            x = median(cca_coef(:,band)); y = median(pears_coef(:,band));
            hold on; scatter(x,y,'+','MarkerEdgeColor',col(band,:),'LineWidth',2);
        end
    else
        trial_data = trim_data(trial_data,'feed');
        [pears_coef, ~] = unit_correlation(trial_data,struct('array',filenames{file,2},'signal','SUA','pca_dims',8,'doPlot',false));
        [cca_coef,~,~] = fCCA(trial_data,struct('array',filenames{file,2},'pca_dims',8,'surrogate_iter',0,'doPlot',false));
        for band = 1:9
            x = median(cca_coef(:,band)); y = median(pears_coef(:,band));
            hold on; scatter(x,y,'^','MarkerEdgeColor',col(band,:),'MarkerFaceColor',col(band,:));
        end
    end
end
xlim([0,0.7]); ylim([0,0.7]); set(gca,'TickDir','out'); box off; 
ylabel('SU correlation'); xlabel('Latent dynamics correlation')

%% Supplementary figures
% Figure 7 - suppl 2
figure
% subplot(1,3,2)
total_cca_exec = {}; total_cca_prep = {}; total_cca_rest = {};
for file = 1:16
    load(filenames{file,1}); 
    trial_data_temp = trim_data(trial_data,'exec');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'prep');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_prep{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); hold on; stdshade(data_exec,0.5,[0.19,0.71,1])
data_prep = cell2mat(total_cca_prep'); hold on; stdshade(data_prep,0.5,[0.12,0.74,0.45])
data_rest = cell2mat(total_cca_rest'); hold on; stdshade(data_rest,0.5,[0.7,0.7,0.7])
set(gca,'xticklabel',bands_name); xtickangle(45); ylabel('LFP-latent dynamics correlation'); ylim([0,0.7])
set(gca,'TickDir','out'); box off; title('M1')

subplot(1,3,1)
total_cca_exec = {}; total_cca_prep = {}; total_cca_rest = {};
for file = 17:28
    load(filenames{file,1}); 
    trial_data_temp = trim_data(trial_data,'exec');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'prep');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_prep{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); hold on; stdshade(data_exec,0.5,[0.19,0.71,1])
data_prep = cell2mat(total_cca_prep'); hold on; stdshade(data_prep,0.5,[0.12,0.74,0.45])
data_rest = cell2mat(total_cca_rest'); hold on; stdshade(data_rest,0.5,[0.7,0.7,0.7])
set(gca,'xticklabel',bands_name); xtickangle(45); ylabel('LFP-latent dynamics correlation'); ylim([0,0.7])
set(gca,'TickDir','out'); box off; title('PMd')

subplot(1,3,3)
total_cca_exec = {}; total_cca_prep = {}; total_cca_rest = {};
for file = 29:36
    load(filenames{file,1}); 
    trial_data_temp = trim_data(trial_data,'feed');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_exec{file} = mean(cca_coeff);
    trial_data_temp = trim_data(trial_data,'rest');
    [cca_coeff,~,~] = fCCA(trial_data_temp,struct('array',filenames{file,2},'pca_dims',10,'surrogate_iter',0,'doPlot',false));
    total_cca_rest{file} = mean(cca_coeff);
end
data_exec = cell2mat(total_cca_exec'); hold on; stdshade(data_exec,0.5,[0.19,0.71,1])
data_rest = cell2mat(total_cca_rest'); hold on; stdshade(data_rest,0.5,[0.7,0.7,0.7])
set(gca,'xticklabel',bands_name); xtickangle(45); ylabel('LFP-latent dynamics correlation'); ylim([0,0.7])
set(gca,'TickDir','out'); box off; title('Area 2')