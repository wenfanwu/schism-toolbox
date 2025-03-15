function bp = read_schism_bp(filepath)
% Read the *bp file.
%
%% Syntax
% varData = read_schism_bp(filepath)
%
%% Description
% varData = read_schism_bp(filepath) reads bp file
%
%% Input Arguments
% filepath - filepath; char
%       absolute filepath of the *.bp file
%
%% Output Arguments
% bp - variable datastruct; double
%       datastruct containing the bp info.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 21 Feb 2025. 
% Email: wwu@vims.edu
% 
% See also: read_schism_prop

%% Parse inputs
if ~strcmp(filepath(end-2:end), '.bp')
    error('the input file must end with bp!')
end

%% Read data
D = importdata(filepath, '%/s', inf);
vm = double(split(string(D(3:end))));

bp.x = vm(:,2);
bp.y = vm(:,3);
bp.c = vm(:,4);
end