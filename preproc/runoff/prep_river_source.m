function D = prep_river_source(river_info, tracer_list)
% Prepare river source data based on "river_info"
% 
%% Syntax
% D = prep_river_source(river_info, tracer_list)
%
%% Description 
% D = prep_river_source(river_info, tracer_list) prepares river source data
%       based on "river_info".
%
%% Examples
% varList = {'runoff', 'temp', 'salt'};
% river_info = add_river_inputs(river_info, Mobj.time, varList, 'real_time');
% D = prep_river_source(river_info, varList);
%
%% Input Arguments
% river_info - river datastruct
%       the updated struct arrays with riverine variables added.
% tracer_list - the tracer list; cell
%       the tracer variables to be used for river inputs;
%
%% Output Arguments
% D - variable datastruct; datastruct.
%       the variable datastruct containing source data.
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 6 May 2025.
% Email: wwu@vims.edu
% 
% See also: match_rivers and add_river_inputs

%% Parse inputs
river_time = river_info(1).time;  % time vector must be the same for all rivers.
tracer_list = lower(tracer_list); tracer_list(contains(tracer_list, 'runoff')) = [];

%% Begin to convert
D.source_elems = [river_info.Elements];
D.sink_elems = [];

D.vsource.runoff = [river_info.runoff];
D.vsource.time = river_time;
D.vsource.dt = unique(seconds(diff(river_time)));

D.vsink.runoff = zeros(size(river_time));  % sink is set to be zero
D.vsink.time = river_time;
D.vsink.dt = unique(seconds(diff(river_time)));

nVars = numel(tracer_list);
for iVar = 1:nVars
    varName = tracer_list{iVar};
    D.msource.(varName)  = [river_info.(varName)];
end
D.msource.time = river_time;
D.msource.dt = unique(seconds(diff(river_time)));

end
