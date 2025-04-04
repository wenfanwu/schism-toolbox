function write_schism_th_nc(Mobj, prefix_name, BdryCnd)
% Write *.th_nc files for SCHISM
% 
%% Syntax
% write_schism_th_nc(Mobj, prefix_name, BdryCnd)
%
%% Description 
% write_schism_th_nc(Mobj, prefix_name, BdryCnd) creates prefix_name.th.nc
% file as boudary inputs.
%
%% Examples
% write_schism_th_nc(Mobj, 'elev2D', BdryCnd)
% write_schism_th_nc(Mobj, 'COS_3D', BdryCnd)
% 
%% Input Arguments
% Mobj --- the mesh object
% prefix_name --- prefix name for the th.nc file, e.g. elev2D, uv3D, TEM_3D, COS_3D.
% BdryCnd --- the datastruct that contains boundary data.
% 
%% Output Arguments
% None
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 3 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
BdryCnd = std_bdry_inputs(Mobj, BdryCnd);

switch prefix_name
    case 'elev2D'
        D = BdryCnd.ssh;
    case 'TEM_3D'
        D = BdryCnd.temp;
    case 'SAL_3D'
        D = BdryCnd.salt;
    case 'uv3D'
        D.time_step = BdryCnd.uvel.time_step;
        D.time = BdryCnd.uvel.time;
        D.time_series = cat(1,BdryCnd.uvel.time_series,BdryCnd.vvel.time_series);
    case 'COS_3D'
        D = merge_tracer_struct(Mobj, BdryCnd, 8);
    case 'ICM_3D'
        D = merge_tracer_struct(Mobj, BdryCnd, 7);
end

filepath = [Mobj.aimpath, prefix_name, '.th.nc'];
write_th_nc(filepath, D)
end

function BdryCnd = std_bdry_inputs(Mobj, BdryCnd)
% Standardized

time = Mobj.time;

varList = fieldnames(BdryCnd);
nVars = numel(varList);

for iVar = 1:nVars
    varName = varList{iVar};
    varData = BdryCnd.(varName);

    BdryCnd = rmfield(BdryCnd, varName);
    
    C.time_step = seconds(time(2)-time(1));
    C.time = seconds(time-time(1));
    varData = permute(varData, [2 1 3]);
    C.time_series(1,:,:,:) = flip(varData, 1);   % the bottom row in the depth dimension denotes the surface layer
    
    BdryCnd.(varName) = C;
    clear C
end
end

function D = merge_tracer_struct(Mobj, BdryCnd, N)

tracer_list = Mobj.tracer_sheet(3:end, N);
ntracers = find(~cellfun(@isempty, tracer_list), 1, 'last');

tmp_time = BdryCnd.(lower(tracer_list{1})).time;

D.time_step = BdryCnd.(lower(tracer_list{1})).time_step;
D.time = tmp_time;
varAll = [];
for ii = 1:ntracers
    tracer_name = lower(tracer_list{ii});
    varTmp = BdryCnd.(tracer_name).time_series;
    varAll = cat(1, varAll, varTmp);
end
D.time_series = varAll;
end

function  write_th_nc(filepath, D)

if exist(filepath,'file')==2
    delete(filepath)
end
if numel(find(isnan(D.time_series(:))))~=0
    error('NaN values exist in the input data!')
end

nComponents = size(D.time_series, 1);
nLevels = size(D.time_series,2);
nOpenBndNodes = size(D.time_series,3);
nTimes = length(D.time);
nccreate(filepath,'time_step','Dimensions',{'one', 1},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'time_step', D.time_step);
nccreate(filepath,'time','Dimensions',{'time', nTimes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', D.time);

nccreate(filepath,'time_series','Dimensions', ...
    {'nComponents', nComponents, 'nLevels', nLevels, 'nOpenBndNodes', nOpenBndNodes, 'time', nTimes},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'time_series', D.time_series);

[~, fileName] = fileparts(filepath);
disp([fileName, '.nc has been successfully created!'])
end


