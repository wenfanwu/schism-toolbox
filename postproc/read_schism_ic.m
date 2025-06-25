function [c, x, y] = read_schism_ic(filepath)
% Read the *.ic file (hvar or vvar).
% .
%% Syntax
% [c, x, y] = read_schism_ic(filepath)
% 
%% Description 
% [c, x, y] = read_schism_ic(filepath) reads data from the *.ic file.
% 
%% Example
% [c, z] = read_schism_ic('D:\ts.ic');  % vvar-type file
% [c, x, y] = read_schism_ic('D:\temp.ic');  % hvar-type file
%
%% Input Arguments
% filepath - filepath; char
%       the absolute filepath of the *.gr3 file.
% 
%% Output Arguments
% c - variable data; numeric
%       the variable data from the *.ic file.
% x - x- or z-coordinate; numeric
%       x- or z-coordinate vector from the *.ic file.
% y - y- or z-coordinate; numeric
%       y-or z-coordinate vector from the *.ic file.
%
%% Notes
% This function can read data from hvar- or vvar-type ic files.
% When reading data from the vvar-type ic file, x and y represent the
% z-coordinate.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 19 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_ic

%% Parse inputs
if ~strcmp(filepath(end-2:end), '.ic')
    error('the input file must end with ic!')
end

%% Read "depth" info
% D = importdata(filepath, '%/s', inf);
% D = cellfun(@(x) strtrim(x), D, 'UniformOutput',false);  % Adapt to earlier versions of MATLAB
fid = fopen(filepath); D = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
D = strtrim(D{1});

head_info = strsplit(D{2});
if numel(head_info)==2  % hvar
    nps = str2double(head_info(2));
    vm = double(split(string(D(3:3+nps-1))));
    vm(:, isnan(sum(vm, 1))) = [ ];
    x = vm(:,2); y = vm(:,3); c = vm(:, 4);
else    % vvar
    nz = str2double(strsplit(strtrim(D{1}(1:find(D{1} == '!', 1)-1))));
    vm = [];
    for iz = 2:nz+1
        if contains(D{iz}, '!')
            tline = double(split(string(strtrim(D{iz}(1:find(D{iz} == '!', 1)-1)))));  % omit the comments
        else
            tline = double(split(string(strtrim(D{iz}))));
        end
        vm = cat(1, vm, tline(2:end)');
    end
    c = vm(:,2:end); z = vm(:,1);
    x = z; y = z;
end

end