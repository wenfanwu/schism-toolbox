function [c, x, y] = read_schism_bp(filepath)
% Read the *.bp file.
%
%% Syntax
% [c, x, y] = read_schism_bp(filepath)
%
%% Description
% [c, x, y] = read_schism_bp(filepath) reads data from the *.bp file
%
%% Input Arguments
% filepath - filepath; char
%       the absolute filepath of the *.bp file
%
%% Output Arguments
% c - variable data; numeric
%       the variable data from the *.bp file.
% x - x-coordinate; numeric
%       x-coordinate vector from the *.bp file.
% y - y-coordinate; numeric
%       y-coordinate vector from the *.bp file.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 19 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: read_schism_prop

%% Parse inputs
if ~strcmp(filepath(end-2:end), '.bp')
    error('the input file must end with bp!')
end

%% Read data
% D = importdata(filepath, '%/s', inf);
fid = fopen(filepath); D = textscan(fid, '%s', 'Delimiter', '\n'); fclose(fid);
D = strtrim(D{1});

vm = double(split(string(D(3:end))));
x = vm(:,2); y = vm(:,3); c = vm(:,4);

end