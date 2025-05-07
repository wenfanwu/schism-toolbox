function [c, x, y] = read_schism_gr3(filepath)
% Read the *.gr3 file.
% 
%% Syntax
% [c, x, y] = read_schism_gr3(filepath)
% 
%% Description 
% [c, x, y] = read_schism_gr3(filepath) reads "depth" info from the *.gr3 file.
% 
%% Example
% filepath = 'E:\watertype.gr3';
% c = read_schism_gr3(filepath);
%
%% Input Arguments
% filepath - filepath; char
%       the absolute filepath of the *.gr3 file.
% 
%% Output Arguments
% c - variable data; numeric
%       the variable data from the *.gr3 file.
% x - x-coordinate; numeric
%       x-coordinate vector from the *.gr3 file.
% y - y-coordinate; numeric
%       y-coordinate vector from the *.gr3 file.
%
%% Notes
% This function aims to read the "depth" info from generic gr3 files (e.g.,
% watertype.gr3, drag.gr3). As for hgrid.gr3, please use
% read_schism_hgrid.m instead. 
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 19 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_gr3

%% Parse inputs
if ~strcmp(filepath(end-2:end), 'gr3')
    error('the input file must end with gr3!')
end
%% Read "depth" info
% D = importdata(filepath, '%/s', inf);
% D = cellfun(@(x) strtrim(x), D, 'UniformOutput',false);  % Adapt to earlier versions of MATLAB
fid = fopen(filepath); D = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
D = strtrim(D{1});

head_info = strsplit(D{2});
nps = str2double(head_info(2));   % # of nodes

% vm = double(split(string(D(3:3+nps-1)))); vm(:, isnan(sum(vm, 1))) = [ ];
vm = reshape(sscanf(strjoin(D(3:3+nps-1)'), '%f'), [], nps)';
x = vm(:,2); y = vm(:,3); c = vm(:, 4);

end