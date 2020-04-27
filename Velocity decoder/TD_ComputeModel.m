function [model] = TD_ComputeModel(trial_data,BinToPast,name,eval)

idx = randperm(length(trial_data));

test_idx = idx(1:round(length(trial_data)*0.1));
train_idx = idx(round(length(trial_data)*0.1)+1:length(trial_data));


% Dup and Shift the PCA projections
if BinToPast > 1
    trial_data = dupeAndShift(trial_data,name,-(1:BinToPast));  
    mod_params.in_signals = [name '_shift'];
elseif BinToPast < -1
    trial_data = dupeAndShift(trial_data,name,(1:abs(BinToPast)));  
    mod_params.in_signals = [name '_shift'];
else
    mod_params.in_signals = name;
end

% getModel will build the wiener cascade. prepare the inputs
mod_params.model_type = 'linmodel';
mod_params.out_signals = 'vel';
mod_params.train_idx = train_idx;
mod_params.polynomial = 3;

trial_data = getModel(trial_data,mod_params);

x_vel = getSig(trial_data(test_idx),{'vel',1});
y_vel = getSig(trial_data(test_idx),{'vel',2});
x_vel_pred = getSig(trial_data(test_idx),{'linmodel_default',1});
y_vel_pred = getSig(trial_data(test_idx),{'linmodel_default',2});


model.actual_vel = getSig(trial_data(test_idx),{'vel',[1,2]});
model.predic_vel = getSig(trial_data(test_idx),{'linmodel_default',[1,2]});


for n = 1:length(eval)
    switch eval{n}
        case 'vaf'
            vaf_x = compute_vaf(x_vel,x_vel_pred);
            vaf_y = compute_vaf(y_vel,y_vel_pred);
            model.vaf = [vaf_x vaf_y];
        case 'r2'
            r2_x = compute_r2(x_vel,x_vel_pred);
            r2_y = compute_r2(y_vel,y_vel_pred);
            model.r2 = [r2_x r2_y];
        otherwise
            error('Not a valid evaluation method: try ''vaf'' or ''r2''');
    end
end


end