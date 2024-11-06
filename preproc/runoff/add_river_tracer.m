function river_info = add_river_tracer(river_info, tracer_list, flux_type)
% 
% 
%% Syntax
% 
%
%% Description 
% 
%
%% Examples
%
%
%% Input Arguments
%
%
%% Output Arguments
% 
% 
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-24.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
nRivers = numel(river_info.river_list);
nTracers = numel(tracer_list);

unit_time = 'month';  % Need improvements

model_time = river_info.river_time;
time_shift = dateshift(model_time, 'start', unit_time);
time_clim = datetime(year(time_shift(1))-1,1,1):calmonths(1):datetime(year(time_shift(end))+1,12,1);

for iRiver = 1:nRivers
    river_name = river_info.river_list{iRiver};
    for iTracer = 1:nTracers
        tracer_name = tracer_list{iTracer};

        [time_raw, tracer_raw] = get_river_tracer(river_name, tracer_name);  % NOTE: define this function on your own
        tracer_clim = accumarray(month(time_raw), tracer_raw, [], @mean);
        tracer_full = repmat(tracer_clim(:)', 1, length(time_clim)/12);

        TT1 = array2timetable(tracer_full(:),'RowTimes', time_clim);
        TT2 = array2timetable(tracer_raw(:),'RowTimes', time_raw);
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
        tracer_new = mean(TA.Variables, 2, 'omitnan');
        time_new = TA.Time;

        river_info.(tracer_name)(:,iRiver) = interp1(time_new, tracer_new, time_shift);
    end
end

end
function [tracer_time, tracer_var] = get_river_tracer(river_name, tracer_name)
% This function aims to extract the tracer data according to the river name.

RivData = load('example_river_data.mat');

D = RivData.(river_name);
tracer_time = D.time;
tracer_var = D.(tracer_name);
end


