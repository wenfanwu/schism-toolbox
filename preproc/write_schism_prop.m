function write_schism_prop(Mobj, prefix_name, prop_flags)
% Create *.prop file for SCHISM.
%
%% Syntax 
% write_schism_prop(Mobj, prefix_name, prop_flags)
%
%% Description
% write_schism_prop(Mobj, prefix_name, prop_flags) writes a
% <prefix_name>.prop file for SCHISM with 'prop_flags'.
%
%% Input Arguments
% Mobj --- mesh object
% prefix_name --- the prefix name of .prop file, e.g. prefix_name = 'tvd'
% prop_flags --- the data to be written.
%
%% Output Arguments
% None
%
%% Notes
% This function aims to generate the input files with a prop suffix for
% SCHISM model, e.g. tvd.prop, fluxflag.prop etc. 
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 4 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: fprintf

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.prop$', '');

fileName = [Mobj.aimpath, prefix_name, '.prop'];
fid = fopen(fileName,'wt');

node_part = [(1:Mobj.nElems)', prop_flags(:)]';
node_fmt = repmat('%d   %d\n', 1, size(node_part,2));
fprintf(fid, node_fmt, node_part(:));

fclose(fid);
disp([prefix_name, '.prop has been created successfully!'])

end
