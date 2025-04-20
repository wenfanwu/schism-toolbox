function prop_flags = read_schism_prop(filepath)
% Read the *.prop file.
%
%% Syntax
% prop_flags = read_schism_prop(filepath)
%
%% Description
% prop_flags = read_schism_prop(filepath) reads data from the *.prop file.
%
%% Examples 
% filepath = 'E:\tvd.prop';
% prop_flags = read_schism_prop(filepath)
%
%% Input Arguments
% filepath - filepath; char
%       the absolute filepath of the *.prop file
%
%% Output Arguments
% prop_flags - prop flags; numeric
%       the variable data from the *.prop file.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2024. 
% Last Updated on 19 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_prop

%% Parse inputs
if ~strcmp(filepath(end-4:end), '.prop')
    error('the input file must end with prop!')
end

%% Begin to read
D = importdata(filepath, '%/s', inf);
vm = double(split(string(D)));
prop_flags = vm(:,2);

end