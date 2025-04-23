function write_schism_nu_nc(Mobj, prefix_name, NdgCnd) 
% Write the *_nu.nc files for SCHISM.
% 
%% Syntax
% write_schism_nu_nc(Mobj, prefix_name, NdgCnd) 
%
%% Description 
% write_schism_th_nc(Mobj, prefix_name, BdryCnd) writes the
%       <prefix_name>_nu.nc file for SCHISM.
%
%% Examples
% bdry_time = Mobj.time(1):Mobj.time(end);
% write_schism_th_nc(Mobj, 'TEM', NdgCnd) 
% write_schism_th_nc(Mobj, 'SAL', NdgCnd) 
% write_schism_th_nc(Mobj, 'ICM', NdgCnd) 
% 
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       prefix name of the *_nu.nc file, e.g., TEM, SAL, COS;
%       the file extension will be omitted automatically.
% NdgCnd - the nudging data; datastruct
%       the datastruct array used to store nudging data.
% 
%% Output Arguments
% None
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 23 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: write_schism_th_nc

%% Parse inputs
prefix_name = regexprep(prefix_name, '\_nu.nc$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '_nu.nc']);

%% Merge variables
switch prefix_name
    case 'TEM'
        varList = {'temp'};
    case 'SAL'
        varList = {'salt'};
    case 'COS'
        varList = Mobj.tracer_sheet(3:end, 8);
    case 'ICM'
        varList = Mobj.tracer_sheet(3:end, 7);
end
% Merge boundary variables into "time_series"
D = merge_nudge_vars(NdgCnd, varList);

%% Begin to write
if exist(filepath,'file')==2; delete(filepath); end

[nTracers, nLevs, nNodes, nTimes] = size(D.tracer_concentration);
nccreate(filepath,'time','Dimensions',{'time', nTimes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', D.time);
nccreate(filepath,'map_to_global_node','Dimensions',{'node', nNodes},'Datatype','int64','Format','netcdf4')
ncwrite(filepath,'map_to_global_node', D.map_to_global_node);
nccreate(filepath,'tracer_concentration','Dimensions', ...
    {'ntracers', nTracers, 'nLevels', nLevs, 'node', nNodes, 'time', nTimes},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'tracer_concentration', D.tracer_concentration);

disp([prefix_name, '_nu.nc has been successfully created!'])
end

function D = merge_nudge_vars(NdgCnd, varList)
% Merge boundary variables into "time_series"
% All variables in "varList" must have the same size (nz*nps*nt)

if ischar(varList); varList = {varList}; end
nVars = numel(varList);  nudge_times = cell(nVars,1); 
for iVar = 1:nVars
    varName = lower(varList{iVar});
    ind_var = find(strcmp({NdgCnd.Variable}, varName), 1);
    % check the existence of variable
    if ~isempty(ind_var)
        nudge_times{iVar} = NdgCnd(ind_var).Time;
        varData = NdgCnd(ind_var).Data;
    else
        error('Variable %s not found in NdgCnd', varName);
    end
    D.tracer_concentration(iVar,:,:,:) = flip(varData, 1);  % the last row indicates the surface layer.
end
% Check the time vector
if all(cellfun(@(x) isequal(x, nudge_times{1}), nudge_times))
    nudge_time = nudge_times{1};
else
    error('"Time" fields must be identical in the same NetCDF file!')
end
D.time = seconds(nudge_time - nudge_time(1));
D.map_to_global_node = NdgCnd(ind_var).Nodes;

% Fill NaNs with junk values
n_nans = numel(find(isnan(D.tracer_concentration)));
if n_nans~=0
    warning on; warning([num2str(n_nums), ' NaNs have been filled with junk values'])
    D.tracer_concentration(isnan(D.tracer_concentration)) = -9999; % Junk values, will not be nudged by SCHISM.
end

end




