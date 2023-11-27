function DS = prep_river_source(river_info, tracer_list)
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
DS.source_elems = river_info.river_elems;
DS.sink_elems = [];

DS.vsource.runoff = river_info.runoff;
DS.vsource.time = river_info.river_time;
DS.vsource.dt = unique(seconds(diff(river_info.river_time)));

DS.vsink.runoff = 0*river_info.runoff;  % sink is set to be zero
DS.vsink.time = river_info.river_time;
DS.vsink.dt = unique(seconds(diff(river_info.river_time)));

 nTracers = numel(tracer_list);
for iTracer = 1:nTracers
    tracer_name = tracer_list{iTracer};
    DS.msource.(tracer_name)  = river_info.(tracer_name);
end
DS.msource.time = river_info.river_time;
DS.msource.dt = unique(seconds(diff(river_info.river_time)));

end
