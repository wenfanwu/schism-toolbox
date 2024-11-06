function varData = read_schism_prop(filepath)
% Read the *prop file.
%
%% Syntax
% varData = read_schism_prop(filepath)
%
%% Description
% varData = read_schism_prop(filepath)
%
%% Examples 
% filepath = 'tvd.prop';
% varData = read_schism_prop(filepath)
%
%% Input Arguments
% filepath - filepath; char
%       absolute filepath of the *.prop file
%
%% Output Arguments
% varData - variable data; double
%       variable vector from the prop file
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 31 Oct 2024. 
% Email: wwu@vims.edu
% 
% See also: read_schism_gr3

%% Parse inputs
if ~strcmp(filepath(end-4:end), '.prop')
    error('the input file must end with prop!')
end

%% Begin to read
D = importdata(filepath, '%/s', inf);

varTmp = double(split(string(D)));
varData = varTmp(:,2);
end