%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [model] = ComputeVelModel(trial_data,bins_to_past,name,eval)
%
%  Will generate model using 90% of trials as training and the remaning 10%
%  as testing. The function will return the actual and predicted velicities
%  as well a the evaluation metric
%
% INPUTS:
%   trial_data  : (struct) trial_data struct
%   bins_to_past: (int) Number of bins of history for the model. 
%                       Negative number indicate "bins to the future" 
%   name        : (char) Name of the field from tial_data use as decoder input
%   eval        : (cell) Evaluation metric.
%                        Can be {'r2'}, {'vaf'} or {'r2','vaf'}
% 
% OUTPUTS:
%   model: (struct) result of the model containing the following fields
%      model.actual_vel    : (array) Tx2 containing the x and y actual velocities
%      model.predic_vel    : (array) Tx2 containing the x and y predicted velocities
%      model.vaf(optional) : (array) 1x2 containing VAF for x and y velocities
%      model.r2(optional)  : (array) 1x2 containing r2 for x and y velocities
%
% Written by Cecilia Gallego-Carracedo. Updated September 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [model] = ComputeVelModel(trial_data,bins_to_past,name,eval)

% Choose train and test trial indexes
idx = randperm(length(trial_data));

test_idx = idx(1:round(length(trial_data)*0.1));
train_idx = idx(round(length(trial_data)*0.1)+1:length(trial_data));


% Duplicate and Shift the data to add extra bins to the past or future
if bins_to_past > 1 % Bins to past
    trial_data = dupeAndShift(trial_data,name,-(1:bins_to_past));  
    mod_params.in_signals = [name '_shift'];
elseif bins_to_past < -1 % Bins to future
    trial_data = dupeAndShift(trial_data,name,(1:abs(bins_to_past)));  
    mod_params.in_signals = [name '_shift'];
else % No shift
    mod_params.in_signals = name;
end

% getModel will build the wiener cascade.
mod_params.model_type = 'linmodel';
mod_params.out_signals = 'vel';
mod_params.train_idx = train_idx;
mod_params.polynomial = 3;

trial_data = getModel(trial_data,mod_params);

% Extract actual and predicted velocities from trial_data
model.actual_vel = getSig(trial_data(test_idx),{'vel',[1,2]});
model.predic_vel = getSig(trial_data(test_idx),{'linmodel_default',[1,2]});

% Compare actual and predicted velocities using the evaluation metrics
for n = 1:length(eval)
    switch eval{n}
        case 'vaf'
            vaf_x = compute_vaf(model.actual_vel(:,1),model.predic_vel(:,1));
            vaf_y = compute_vaf(model.actual_vel(:,2),model.predic_vel(:,2));
            model.vaf = [vaf_x vaf_y];
        case 'r2'
            r2_x = compute_r2(model.actual_vel(:,1),model.predic_vel(:,1));
            r2_y = compute_r2(model.actual_vel(:,2),model.predic_vel(:,2));
            model.r2 = [r2_x r2_y];
        otherwise
            error('Not a valid evaluation method: try ''vaf'' or ''r2''');
    end
end


end