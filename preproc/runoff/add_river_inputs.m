function river_info = add_river_inputs(river_info, river_time, varList, flux_type)
% Add river runoff at the source elements
%
%% Syntax
% river_info = add_river_inputs(river_info, river_time)
% river_info = add_river_inputs(river_info, river_time, varList)
% river_info = add_river_inputs(river_info, river_time, varList, flux_type)
%
%% Description
% river_info = add_river_runoff(river_info, model_time) adds runoff data
%       for the source elements. 
% river_info = add_river_inputs(river_info, river_time, varList) specifies
%       the variable list.
% river_info = add_river_inputs(river_info, river_time, varList, flux_type)
%       specifies the types (clim. or real-time) of river inputs.
%
%% Examples
% varList = {'runoff', 'temp', 'salt'};
% river_info = add_river_inputs(river_info, Mobj.time, varList, 'real_time');
%
%% Input Arguments
% river_info - river datastruct; datastruct.
%       the struct arrays containing river info, which typically comes from
%       the 'match_rivers' function.
% river_time  --- river time; datetime.
%       time vector (datetime) for the river inputs. 
%       Default: river_time = Mobj.time;
% varList - variable list; cell
%       the variables to be used for river inputs. 
%       Default: varList = {'runoff', 'temp', 'salt'};
% flux_type - flux type; char
%       the flux type of river inputs (clim or real_time). 
%       Default: flux_type = 'real_time'.
%
%% Output Arguments
% river_info - river datastruct
%       the updated struct arrays with riverine variables added.
%
%% Notes
% If the provided river data does not fully cover the simulation period,
% this function will automatically fill the uncovered time range using
% climatological data from the available period.
%
% When flux_type = 'clim', the function will return the monthly
% climatological data.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 6 May 2025.
% Email: wwu@vims.edu
%
% See also: match_rivers

%% Parse inputs
if nargin < 2; river_time = Mobj.time; end
if nargin < 3; varList = {'runoff', 'temp', 'salt'}; end
if nargin < 4; flux_type = 'real_time'; end
if ischar(varList); varList = {varList}; end

%% Define time vectors
unit_time = 'month';  % climatology monthly mean
time_shift = dateshift(river_time, 'start', unit_time);
time_clim = datetime(year(time_shift(1))-1,1,1):calmonths(1):datetime(year(time_shift(end))+1,12,1);

%% Extract runoff data.
nRivers = length(river_info);
nVars = numel(varList);

for iRiver = 1:nRivers
    river_info(iRiver).time = river_time(:);
    river_name = river_info(iRiver).River;
    for iVar = 1:nVars
        varName = varList{iVar};
        [time_raw, var_raw] = get_river_data(river_name, varName); 
        var_clim = accumarray(month(time_raw(:)), var_raw(:), [], @mean);
        var_full = repmat(var_clim(:)', 1, length(time_clim)/12);

        TT1 = array2timetable(var_full(:),'RowTimes', time_clim);
        TT2 = array2timetable(var_raw(:),'RowTimes', time_raw);
        switch flux_type
            case 'clim'
                TA = TT1;
            case {'real_time', 'real'}
                TA = synchronize(TT1, TT2); % fill missing data with clim. values
                ind_nan = ~isnan(TA.Var1_TT2);
                TA.Var1_TT1(ind_nan) = nan;
            otherwise
                error('Invalid flux type!')
        end
        var_new = mean(TA.Variables, 2, 'omitnan'); time_new = TA.Time;

        n_outlets = max(1, sum(strcmp({river_info.River}, river_name))*strcmpi(varName, 'runoff')); % Runoff distribution
        river_info(iRiver).(varName) = interp1(time_new, var_new, time_shift)/n_outlets;
    end
end

end

function [time, varData] = get_river_data(river, var)
% This function extracts runoff data for a specified river name.
% Note: The river dataset should be prepared based on your specific application.
% The returned variable 'time' must be in the datetime format.
% Also, ensure that the river name is used consistently across the entire toolbox.

RivData = load('example_river_data.mat');
D = RivData.(river);  time = D.time;  varData = D.(var);

end

















