function write_schism_th_nc(Mobj, prefix_name, BdryCnd)
% Write *th.nc files for SCHISM.
% 
%% Syntax
% write_schism_th_nc(Mobj, prefix_name, BdryCnd)
%
%% Description 
% write_schism_th_nc(Mobj, prefix_name, BdryCnd) writes the
%       <prefix_name>.th.nc file for SCHISM.
%
%% Examples
% bdry_time = Mobj.time(1):Mobj.time(end);
% write_schism_th_nc(Mobj, 'elev2D', BdryCnd)
% write_schism_th_nc(Mobj, 'COS_3D', BdryCnd)
% 
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       prefix name of the *th.nc file, e.g., elev2D, uv3D, TEM_3D, COS_3D;
%       the file extension will be omitted automatically.
% BdryCnd - the boundary data; datastruct
%       the datastruct array used to store boundary inputs.
% 
%% Output Arguments
% None
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 21 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: write_schism_th

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.th.nc$', '');
filepath = [Mobj.aimpath, prefix_name, '.th.nc'];

%% Merge variables
switch prefix_name
    case 'elev2D'
        varList = {'ssh'};
    case 'TEM_3D'
        varList = {'temp'};
    case 'SAL_3D'
        varList = {'salt'};
    case 'uv3D'
        varList = {'uvel', 'vvel'};
    case 'COS_3D'
        varList = Mobj.tracer_sheet(3:end, 8);
    case 'ICM_3D'
        varList = Mobj.tracer_sheet(3:end, 7);
end
% Merge boundary variables into "time_series"
D = merge_bdry_vars(BdryCnd, varList);

%% Begin to write
if exist(filepath,'file')==2; delete(filepath); end  % over-write
if any(isnan(D.time_series(:))); error('NaNs were found in the input data!'); end

[nComponents, nLevels, nOpenBndNodes, nTimes] = size(D.time_series);

nccreate(filepath,'time_step','Dimensions',{'one', 1},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'time_step', D.time_step);
nccreate(filepath,'time','Dimensions',{'time', nTimes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', D.time);
nccreate(filepath,'time_series','Dimensions', ...
    {'nComponents', nComponents, 'nLevels', nLevels, 'nOpenBndNodes', nOpenBndNodes, 'time', nTimes},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'time_series', D.time_series);

disp([prefix_name, '.nc has been successfully created!'])
end

function D = merge_bdry_vars(BdryCnd, varList)
% Merge boundary variables into "time_series"
% All variables in "varList" must have the same size (nz*nps*nt)

if ischar(varList); varList = {varList}; end
nVars = numel(varList);  bdry_times = cell(nVars,1);
for iVar = 1:nVars
    varName = lower(varList{iVar});
    ind_var = find(strcmp({BdryCnd.Variable}, varName), 1);

    % check the existence of variable
    if ~isempty(ind_var)
        bdry_times{iVar} = BdryCnd(ind_var).Time;
        varData = BdryCnd(ind_var).Data;
    else
        error('Variable %s not found in BdryCnd', varName);
    end
    % clip to valid ranges
    varData = check_schism_var(varData, varName);
    
    D.time_series(iVar,:,:,:) = flip(varData, 1);  % the last row indicates the surface layer.
end

% check the time vector
if all(cellfun(@(x) isequal(x, bdry_times{1}), bdry_times))
    bdry_time = bdry_times{1};
else
    error('"Time" fields must be identical in the same NetCDF file!')
end
D.time_step = seconds(bdry_time(2) - bdry_time(1));
D.time = seconds(bdry_time - bdry_time(1));

end
