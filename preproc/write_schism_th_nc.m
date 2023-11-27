function write_schism_th_nc(Mobj, BdryCnd, prefix_name)
% Write .th_nc file for SCHISM
% 
%% Syntax
% write_schism_th_nc(Mobj, BdryCnd, prefix_name)
%
%% Description 
% write_schism_th_nc(Mobj, BdryCnd, prefix_name) creates prefix_name.th.nc
% file as boudary inputs.
%
%% Examples
% write_schism_th_nc(Mobj, BdryCnd, 'elev2D')
% write_schism_th_nc(Mobj, BdryCnd, 'COS_3D')
% 
%% Input Arguments
% Mobj --- the mesh object
% BdryCnd --- the datastruct that contains boundary data.
% prefix_name --- prefix name for the th.nc file, e.g. elev2D, uv3D,
% TEM_3D, COS_3D.
% 
%% Output Arguments
% None
% 
%% Notes
% Make sure the provided BdryCnd contains required variables. 
% This function is still under development to accommodate more tracer 
% modules.
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022. 
% Last Updated on 2023-11-27.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
% sz_ref = [Mobj.maxLev, numel(Mobj.obc_nodes_tot), numel(Mobj.time)];

% varList = fieldnames(BdryCnd);
% if numel(contains(varList, 'time'))==0
%     error('the input datastruct should have the field named time!')
% end
% varList(contains(varList, 'time')) = [];
% for ii = 1:numel(varList)
%     varName = varList{ii};
%     sz_vars = size(BdryCnd.(varName));
%     if (sz_vars(1) ~= sz_ref(1) && sz_vars(1) ~= 1) || sz_vars(2) ~= sz_ref(2) || sz_vars(3) ~= sz_ref(3) 
%         error(['The size of ', varName, ' is not correct!'])
%     end
% end
%% Write NetCDF
BdryCnd = std_bdry_inputs(Mobj, BdryCnd);

if strcmp(prefix_name, 'elev2D')
    D = BdryCnd.ssh;
end
if strcmp(prefix_name, 'TEM_3D')
    D = BdryCnd.temp;
end
if strcmp(prefix_name, 'SAL_3D')
    D = BdryCnd.salt;
end
if strcmp(prefix_name, 'uv3D')
    D.time_step = BdryCnd.uvel.time_step;
    D.time = BdryCnd.uvel.time;
    D.time_series = cat(1,BdryCnd.uvel.time_series,BdryCnd.vvel.time_series);
end
if strcmp(prefix_name, 'ICM_3D')
    D = merge_tracer_struct(Mobj, BdryCnd, 7);
end
if strcmp(prefix_name, 'COS_3D')
    D = merge_tracer_struct(Mobj, BdryCnd, 8);
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
    
    S.time_step = seconds(time(2)-time(1));
    S.time = seconds(time-time(1));
    varData = permute(varData, [2 1 3]);
    S.time_series(1,:,:,:) = flip(varData, 1);   % the first row of the depth dimension denotes the bottom in SCHISM inputs
    
    BdryCnd.(varName) = S;
    clear S
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


