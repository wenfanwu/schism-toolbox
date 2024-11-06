function [varData, lon, lat] = read_schism_gr3(filepath)
% Read the *gr3 or *ic file.
% 
%% Syntax
% [varData, lon, lat] = read_schism_gr3(filepath)
% 
%% Description 
% [varData, lon, lat] = read_schism_gr3(filepath) reads "depth" info from
% the gr3 file.
% 
%% Example
% filepath = 'watertype.gr3';
% varData = read_schism_gr3(filepath)
%
%% Input Arguments
% filepath - filepath; char
%       the absolute filepath of the *.gr3 file.
% 
%% Output Arguments
% varData - variable data; double
%       the variable data@nodes
% lon - longitude vector; double
%       longitude vector from the gr3 file.
% lat - latitude vector; double
%       latitude vector from the gr3 file.
%
%% Notes
% This function aims to read the "depth" info from generic gr3 files (e.g.,
% watertype.gr3, drag.gr3). As for hgrid.gr3, please use
% read_schism_hgrid.m to load its complete data instead.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 30 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: read_schism_hgrid

%% Parse inputs
if ~strcmp(filepath(end-2:end), 'gr3') & ~strcmp(filepath(end-2:end), '.ic')
    error('the input file must end with gr3 or ic!')
end

D = importdata(filepath, '%/s', inf);
D = cellfun(@(x) strtrim(x), D, 'UniformOutput',false);  % Adapt to earlier versions of MATLAB
%% Read "depth" info
head_info = strsplit(D{2});
nNodes = str2double(head_info(2));

node_part = double(split(string(D(3:3+nNodes-1))));
node_part(:, isnan(sum(node_part, 1))) = [ ];

lon = node_part(:,2);
lat = node_part(:,3);
varData = node_part(:, 4);

end