function river_info = add_river_runoff(river_info, model_time, flux_type)
% Add river runoff on the source elements
% 
%% Syntax
% river_info = add_river_runoff(river_info, model_time)
% river_info = add_river_runoff(river_info, model_time, flux_type)
% 
%% Description 
% river_info = add_river_runoff(river_info, model_time) adds runoff data
% for the source elements.
% river_info = add_river_runoff(river_info, model_time, flux_type) decides
% to use climatological or real-time runoff.
%
%% Examples
% river_info = add_river_runoff(river_info, Mobj.time, 'real_time');
%
%% Input Arguments
% river_info --- the datastruct containing river info, which typically
% results from the 'match_rivers' function.
% model_time --- time vector (datetime) for the simulation.model_time =
% Mobj.time for most cases.
% flux_type --- a string used to decides the flux type (real_time/clim).
% Default: flux_type = 'real_time'.  
%
%% Output Arguments
% river_info --- the updated datastruct containing river info, with runoff added
% 
%% Notes
% If your provided runoff data doesn't fully cover the timespan of your
% simulation, this function will automatically pad the uncovered part with
% climatological data during the available runoff period.
%
% when flux_type = 'clim', you will get the monthly clim. river runoff.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: match_rivers and add_river_tracer

%% Parse inputs
if nargin < 3
    flux_type = 'real_time';
end

nRivers = numel(river_info.river_list);
river_info.river_time = model_time(:);

unit_time = 'month';  % Need imporvements

time_shift = dateshift(model_time, 'start', unit_time);
time_clim = datetime(year(time_shift(1))-1,1,1):calmonths(1):datetime(year(time_shift(end))+1,12,1);

for iRiver = 1:nRivers
    river_name = river_info.river_list{iRiver};
    outlets_num = sum(strcmp(river_info.river_list, river_name));

    [time_raw, runoff_raw] = get_river_runoff(river_name);  % NOTE: define this function on your own

    runoff_clim = accumarray(month(time_raw), runoff_raw, [], @mean); 
    runoff_full = repmat(runoff_clim(:)', 1, length(time_clim)/12);

    TT1 = array2timetable(runoff_full(:),'RowTimes', time_clim);
    TT2 = array2timetable(runoff_raw(:),'RowTimes', time_raw);
    switch flux_type
        case 'clim'
            TA = TT1;
        case {'real_time', 'real'}
            TA = synchronize(TT1, TT2);
            ind_nan = ~isnan(TA.Var1_TT2);
            TA.Var1_TT1(ind_nan) = nan;
        otherwise
            error('Invalid flux type!')
    end
    runoff_new = mean(TA.Variables, 2, 'omitnan');
    time_new = TA.Time;

    river_info.runoff(:,iRiver) = interp1(time_new, runoff_new, time_shift)/outlets_num;
end

end

function [time, runoff] = get_river_runoff(river_name)
% This function aims to extract the runoff data according to the river
% name. Note that the river data should be created based on your own
% applications. The returned 'time' should be in the format of datetime. In
% addiiton, make sure the river name is the consistent within the whole toolbox. 

RivData = load('example_river_data.mat');

D = RivData.(river_name);
time = D.time;
runoff = D.runoff;

end

















