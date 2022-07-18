%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function trial_data = trim_data(trial_data,epoch)
%
%       This function trims the data into same-length trials at the event. 
%       Time points are selected for 30ms bin size
%   
%
% INPUTS:
%   trial_data : (struct) trial_data struct
%   epoch      : (char) behavioural epoch
%       This can be: exec(movement execution)
%                    feed(feedback processing)
%                    prep(movement preparation)
%                    rest
%
% OUTPUTS:
%   trial_data : (struct) trial_data after trimming
%
%
% Written by Cecilia Gallego-Carracedo. Updated September 2021.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function trial_data = trim_data_v2(trial_data,epoch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch epoch
        case 'exec'
            trial_data = trimTD(trial_data, {'idx_movement_on',0},{'idx_movement_on',15});
        case 'feed'
            trial_data = trimTD(trial_data, {'idx_goCueTime',5},{'idx_goCueTime',20});
        case 'prep'
            trial_data = trimTD(trial_data,{'idx_movement_on',-15},{'idx_movement_on',0});
        case 'rest'
            trial_data = trimTD(trial_data, {'idx_tgtOnTime',-15},{'idx_tgtOnTime',0});

            % Check than resting data does not include preparation time
            trial_data([trial_data.idx_tgtOnTime] ~= 16) = [];
            
        otherwise
            disp('Error! Trial epoch not found'); quit;
    end
    
    
end