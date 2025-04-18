function DS = get_hycom_init(Mobj, meta_data)
% Extract real-time hycom data as the initial field.
%  
%% Syntax
% DS = get_hycom_init(Mobj)
% DS = get_hycom_init(Mobj, meta_data)
%
%% Description 
% DS = get_hycom_init(Mobj) extracts hycom data as the initial field
% DS = get_hycom_init(Mobj, meta_data) specifies the medadata.
%
%% Example
% meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.mat';
% meta_data.raw_list = {'ssh','salt','temp'};      % original variable name in the file
% meta_data.new_list = {'ssh','salt','temp'};      % standard variable name for output
% meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
% meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
% DS = get_hycom_init(Mobj, meta_data);
%
%% Input Arguments
% Mobj - the mesh object; datastruct
%       the datastruct used to store the mesh info.
% meta_data -  the metadata; datastruct
%       the metadata of database.
%
%% Output Arguments
% DS - the data struct; datastruct
%       the datastruct used to store the initial data.
% 
%% Notes
% All hycom data are downloaded through "get_hycom_online" and they are
% organized as follows:
%       Wn79En72S33N41_19930101T0000Z.mat
%       Wn79En72S33N41_19930102T0000Z.mat
%       Wn79En72S33N41_19930103T0000Z.mat
%       ......
%       Wn79En72S33N41_20200101T0000Z.mat
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2022.
% Last Updated on 18 Apr 2025.
% Email: wwu@vims.edu
%
% See also: 

%% Define the meta data
if nargin < 2
    meta_data.src_file = 'E:\HYCOM\HYCOM_ChesBay\Wn79En72S33N41_****.mat';  % NEED TO BE CHANGED
    meta_data.raw_list = {'ssh', 'salt', 'temp'};      % original variable name in the file
    meta_data.new_list = {'ssh', 'salt', 'temp'};      % standard variable name for output
    meta_data.dim_vars = {'lon', 'lat', 'depth'};    % dimensional variables in the file
    meta_data.date_fmt = 'yyyymmddTHHMMZ';  % the date format of filename
end
%% Load data
filepath = strrep(meta_data.src_file, '****', datestr(Mobj.time(1), meta_data.date_fmt)); %#ok<DATST>

C = load(filepath);
nVars = numel(meta_data.raw_list);
for iVar = 1:nVars
    raw_name = meta_data.raw_list{iVar};
    new_name = meta_data.new_list{iVar};

    D = struct();
    D.lon = C.(meta_data.dim_vars{1});
    D.lat = C.(meta_data.dim_vars{2});
    D.depth = C.(meta_data.dim_vars{3});
    D.time = Mobj.time(1);
    D.var = squeeze(C.(raw_name));

    DS.(new_name) = D;
end

end
