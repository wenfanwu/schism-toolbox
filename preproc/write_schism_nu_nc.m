function write_schism_nu_nc(Mobj, prefix_name, D)
% Write nudging files for SCHISM (TEM/SAL/ICM_nu.nc)
%
%% Syntax
% write_schism_nu_nc(Mobj, prefix_name, D)
%
%% Description
% write_schism_nu_nc(Mobj, prefix_name, D) writes nudging files (*_nu.nc)
%
%% Examples 
% 
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       this datastruct contains the mesh info.
% prefix_name - filename prefix; char
%       prefix name of the NetCDF files (TEM/SAL/ICM).
% D - input data; datastruct
%       a datastruct containing data to be written
%
%% Output Arguments
% None
%
%% Notes
% The bottom row in the vertical dimension indicates surface layer in the field "tracer_concentration"
%       squeeze(D.tracer_concentration(1,end,100,:)); ===> first tracer at surface layer for 100th node
%       squeeze(D.tracer_concentration(1,1,100,:));  ===> first tracer at bottom layer for 100th node
% 
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021. 
% Last Updated on 29 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
if ~isfield(D, 'time') || ~isfield(D, 'map_to_global_node') || ~isfield(D, 'tracer_concentration')
    error('missing fields in the input datastruct!')
end

if numel(size(D.tracer_concentration))~=4
    error('the dimension of tracer_concentration should be 4')
end

D.tracer_concentration(isnan(D.tracer_concentration)) = -9999; % Junk values, will not be nudged by SCHISM.

[nTracers, nLevs, nNodes, nTimes] = size(D.tracer_concentration);

if nNodes ~= length(D.map_to_global_node)
    error('the # of nudging nodes is inconsistent!')
end

if nTimes ~=  length(D.time)
    error('the nudging time is inconsistent!')
end

%% Begin to write
filepath = [Mobj.aimpath, prefix_name, '_nu.nc'];
if exist(filepath,'file')==2; delete(filepath); end

nccreate(filepath,'time','Dimensions',{'time', nTimes},'Datatype','double','Format','netcdf4')
ncwrite(filepath,'time', D.time);
nccreate(filepath,'map_to_global_node','Dimensions',{'node', nNodes},'Datatype','int64','Format','netcdf4')
ncwrite(filepath,'map_to_global_node', D.map_to_global_node);
nccreate(filepath,'tracer_concentration','Dimensions', ...
    {'ntracers', nTracers, 'nLevels', nLevs, 'node', nNodes, 'time', nTimes},'Datatype','single','Format','netcdf4')
ncwrite(filepath,'tracer_concentration', D.tracer_concentration);

[~, headName,~] = fileparts(filepath);
disp([headName, '.nc has been successfully created!'])
end
