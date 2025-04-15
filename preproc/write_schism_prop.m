function write_schism_prop(Mobj, prefix_name, prop_flags)
% Write the *.prop file for SCHISM.
%
%% Syntax
% write_schism_prop(Mobj, prefix_name, prop_flags)
%
%% Description
% write_schism_prop(Mobj, prefix_name, prop_flags) writes the
%       <prefix_name>.prop file using 'prop_flags'. 
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       prefix name of the *.prop file, e.g. prefix_name = 'tvd', or
%       'tvd.prop'; the file extension will be omitted automatically.
% prop_flags - input data; numeric
%       the input data to be written. If 'prop_flags' is a scalar, it will
%       be expanded to a vector automatically , and the data at all nodes
%       are the same.in this case.
%
%% Output Arguments
% None
%
%% Notes
% This function aims to generate the input files ending with "prop" for
% SCHISM model, e.g. tvd.prop, fluxflag.prop etc. 
%
%% Author Info
% Created by Wenfan Wu,Virginia Institute of Marine Science in 2021. 
% Last Updated on 15 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_gr3

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.prop$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '.prop']);

%% Check inputs
if isscalar(prop_flags); prop_flags = prop_flags*ones(Mobj.nElems, 1); end
if any(isnan(prop_flags(:))); error('NaNs were found in the input data!'); end

%% Begin to write
fid = fopen(filepath,'w');
fprintf(fid, '%d   %d\n', [(1:Mobj.nElems)', prop_flags(:)]');
fclose(fid);

disp([prefix_name, '.prop has been created successfully!'])
end
