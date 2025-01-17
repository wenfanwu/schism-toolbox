function write_schism_gr3(Mobj, prefix_name, varData)
% Write *.gr3 file for SCHISM
%
%% Syntax 
% write_schism_gr3(Mobj, prefix_name, varData)
%
%% Description
% write_schism_gr3(Mobj, prefix_name, varData) writes a <prefix_name>.gr3
% file for SCHISM with 'varData'.
%
%% Input Arguments
% Mobj --- the mesh object
% prefix_name --- prefix name of the gr3 file. e.g. prefix_name = 'drag'.
% varData --- a scalar or vector of data to be written. If the 'varData' is a
% scalar, it will be automatically expanded to a vector, and the data at
% all nodes will the same.
%
%% Notes
% this function aims to generate the input files with the suffix
% "gr3", eg. watertype.gr3, albedo.gr3, windrot_geo2proj.gr3,
% diffmin.gr3, diffmax.gr3, bdef.gr3, windfactor.gr3, hdif.gr3, and so on.
%
%% Author Info
% Created by Wenfan Wu,Virginia Institute of Marine Science in 2021. 
% Last Updated on 17 Sep 2024. 
% Email: wwu@vims.edu
% 
% See also: fprintf

%% Parse inputs
fileName = [Mobj.aimpath, prefix_name, '.gr3'];
headLine = datestr(now); %#ok<DATST>

if isscalar(varData)
    varData = varData*ones(1, Mobj.nNodes);
end
if numel(find(isnan(varData(:))))~=0
    error('There are NaNs in the input data!')
end

fid = fopen(fileName,'wt');
fprintf(fid, [headLine, '\n']);                                                              
fprintf(fid, [num2str(Mobj.nElems),' ',num2str(Mobj.nNodes), '\n']); 

% Node Info
node_part = [(1:Mobj.nNodes)', Mobj.lon(:), Mobj.lat(:), varData(:)]';
node_fmt = repmat('%d   %14.6f   %14.6f    %13.7e\n', 1, size(node_part,2));
fprintf(fid, node_fmt, node_part(:));

tri = Mobj.tri; 
if size(tri, 1)~=Mobj.nElems; tri = tri'; end

% Elem Info
elem_part = [(1:Mobj.nElems)', Mobj.i34(:).*ones(Mobj.nElems,1), tri]';
elem_part_1d = elem_part(:);
elem_part_1d(isnan(elem_part_1d)) = [];

elem_fmt = repmat('%d %d %d %d %d %d\n', Mobj.nElems, 1);
elem_fmt(Mobj.i34==3, :) = repmat('%d %d %d %d %d   \n', sum(Mobj.i34==3), 1);
elem_fmt = elem_fmt'; elem_fmt_1d = elem_fmt(:);

fprintf(fid, elem_fmt_1d, elem_part_1d);

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