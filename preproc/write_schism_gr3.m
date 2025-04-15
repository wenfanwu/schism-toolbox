function write_schism_gr3(Mobj, prefix_name, varData)
% Write the *.gr3 file for SCHISM
%
%% Syntax 
% write_schism_gr3(Mobj, prefix_name, varData)
%
%% Description
% write_schism_gr3(Mobj, prefix_name, varData) writes a <prefix_name>.gr3
%       file for SCHISM using 'varData'.
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% prefix_name - prefix name; char
%       prefix name of the *gr3 file, e.g. prefix_name = 'drag', or
%       'drag.gr3'; the file extension will be omitted automatically.
% varData - input data; numeric
%       the input data to be written. If '' is a scalar, it will
%       be expanded to a vector automatically , and the data at all nodes
%       are the same.in this case.
% 
%% Output Arguments
% None
%
%% Notes
% This function aims to generate the input files ending with
% "gr3", e.g., watertype.gr3, albedo.gr3, windrot_geo2proj.gr3,
% diffmin.gr3, diffmax.gr3, bdef.gr3, windfactor.gr3, hdif.gr3, and so on.
%
%% Author Info
% Created by Wenfan Wu,Virginia Institute of Marine Science in 2021. 
% Last Updated on 15 Apr 2025. 
% Email: wwu@vims.edu
% 
% See also: write_schism_prop and write_schism_ic

%% Parse inputs
prefix_name = regexprep(prefix_name, '\.gr3$', '');
filepath = fullfile(Mobj.aimpath, [prefix_name, '.gr3']);
head_line = datestr(now, 'mmm/dd/yyyy HH:MM:SS'); %#ok<TNOW1,DATST>

%% Check inputs
if isscalar(varData); varData = varData*ones(Mobj.nNodes, 1); end
if any(isnan(varData(:))); error('NaNs were found in the input data!'); end

%% Begin to write
fid = fopen(filepath,'w');
fprintf(fid, [head_line, '\n']);
fprintf(fid,'%d %d\n', Mobj.nElems, Mobj.nNodes);

% Node and Elem
fprintf(fid, '%d   %14.6f   %14.6f    %13.7e\n', [(1:Mobj.nNodes)', Mobj.lon(:), Mobj.lat(:), varData(:)]');
fprintf(fid, '%d %d %d %d %d %d\n', [(1:Mobj.nElems)', Mobj.i34(:), Mobj.tri]');  % include NaN values
fclose(fid);

% Remove NaNs
fileText = fileread(filepath);
fid = fopen(filepath, 'w');
fwrite(fid, strrep(fileText, 'NaN', '')); % case-sensitive
fclose(fid);

disp([prefix_name, '.gr3 has been created successfully!'])
end

% Notes from Wenfan Wu
% ----------watertype.gr3
% Read in water type; the values for R, d_1, d_2 are given below 
% solar flux= R*exp(z/d_1)+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
% 1: 0.58 0.35 23 (Jerlov type I)
% 2: 0.62 0.60 20 (Jerlov type IA)
% 3: 0.67 1.00 17 (Jerlov type IB)
% 4: 0.77 1.50 14 (Jerlov type II)
% 5: 0.78 1.40 7.9 (Jerlov type III)
% 6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
% 7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
% ----------albedo.gr3
% it is generally set to be 0.06;
% ----------windrot_geo2proj.gr3
%  zero is ok for most cases;

% Elem part
% elem_part = [(1:Mobj.nElems)', Mobj.i34(:).*ones(Mobj.nElems,1), Mobj.tri]';
% elem_part_1d = elem_part(:);
% elem_part_1d(isnan(elem_part_1d)) = [];

% elem_fmt = repmat('%d %d %d %d %d %d\n', Mobj.nElems, 1);
% elem_fmt(Mobj.i34==3, :) = repmat('%d %d %d %d %d   \n', sum(Mobj.i34==3), 1);
% elem_fmt = elem_fmt'; elem_fmt_1d = elem_fmt(:);

% fprintf(fid, elem_fmt_1d, elem_part_1d);