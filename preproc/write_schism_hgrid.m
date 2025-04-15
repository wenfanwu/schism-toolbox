function write_schism_hgrid(Mobj)
% Write the hgrid.gr3 & hgrid.ll files for SCHISM.
%
%% Syntax
% write_schism_hgrid(Mobj)
% 
%% Description
% write_schism_hgrid(Mobj) writes the hgrid.gr3 & hgrid.ll files for SCHISM.
%
%% Example
% mesh_file = 'Exp1_BYS\inputs\BYES_34651.mat';
% Mobj = mesh2schism(mesh_file);
% write_schism_hgrid(Mobj)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 15 Apr 2025.
% Email: wwu@vims.edu
%
% See also: write_schism_gr3

%% Parse inputs
hgrid_gr3 = fullfile(Mobj.aimpath, 'hgrid.gr3');
hgrid_ll = fullfile(Mobj.aimpath, 'hgrid.ll');

%% Begin to write
if strncmpi(Mobj.coord, 'geographic', 3)
    write_hgrid_file(Mobj, hgrid_gr3)
    copyfile(hgrid_gr3, hgrid_ll)
else
    write_hgrid_file(Mobj, hgrid_gr3)
    write_hgrid_file(Mobj, hgrid_ll)
end

disp('hgrid.gr3 & hgrid.ll have been created successfully!')
end


function write_hgrid_file(Mobj, filepath)
% write hgrid files (hgrid.gr3 or hgrid.ll)
%
% Notes: Mobj.lon and Mobj.lat stores the coordinates for your given
% coordinate. In other words, they can also represent cartersian
% coordinates when your specified coordinate is 'cartersian".

switch filepath(end-2:end)
    case 'gr3'  % geographic (lon/lat)
        ux = Mobj.lon(:); uy = Mobj.lat(:);
    case '.ll' % cartersian
        ux = Mobj.x(:); uy = Mobj.y(:); % Add new fields "x" and "y" in this case
end

%%  Mesh info (Node & Elem)
fid = fopen(filepath,'w');
fprintf(fid, [Mobj.expname, ' ', datestr(now, 'mmm/dd/yyyy HH:MM:SS'), '\n']);          %#ok<TNOW1,DATST> % alphanumeric description; ignored by code
fprintf(fid,'%d %d\n', Mobj.nElems, Mobj.nNodes);       % # of elements and nodes in the horizontal grid

fprintf(fid, '%d   %14.6f   %14.6f    %13.7e\n', [(1:Mobj.nNodes)', ux(:), uy(:), Mobj.depth(:)]');
fprintf(fid, '%d %d %d %d %d %d\n', [(1:Mobj.nElems)', Mobj.i34(:), Mobj.tri]');  % include NaN values

%% List of open and land boundary segments (needed for hgrid.gr3 only; not needed for other *.gr3)
% ========= WRITE THE OPEN BOUNDARY PART ========= 
fprintf(fid, [num2str(Mobj.obc_counts), ' = Number of open boundaries\n']);
fprintf(fid, [num2str(numel(Mobj.obc_nodes_tot)), ' = Total number of open boundary nodes\n']);
for idx = 1:Mobj.obc_counts
    tmp_nodes = Mobj.obc_nodes(:,idx); tmp_nodes(tmp_nodes==0) = [];
    fprintf(fid, [num2str(length(tmp_nodes),'%d'), ' = Number of nodes for open boundary ', num2str(idx),'\n']);
    fprintf(fid, '%d\n', tmp_nodes(:)');
end

% ========= WRITE THE LAND PART ========= 
fprintf(fid, [num2str(Mobj.land_counts+Mobj.island_counts), ' = number of land boundaries\n']);   % including island numbers here
fprintf(fid, [num2str(numel(Mobj.land_nodes_tot)+numel(Mobj.island_nodes_tot)), ' = Total number of land boundary nodes\n']);
for idx = 1:Mobj.land_counts
    tmp_nodes = Mobj.land_nodes(:,idx); tmp_nodes(tmp_nodes==0) = [];
    fprintf(fid, [num2str(length(tmp_nodes),'%d'), ' 0 = Number of nodes for land boundary ',num2str(idx),'\n']);
    fprintf(fid, '%d\n', tmp_nodes(:)');
end

% ========= WRITE THE ISLAND PART ========= 
for idx = 1:Mobj.island_counts
    tmp_nodes = Mobj.island_nodes(:,idx); tmp_nodes(tmp_nodes==0) = [];
    fprintf(fid, [num2str(length(tmp_nodes),'%d'), ' 1 = Number of nodes for island boundary ',num2str(idx),'\n']);
    fprintf(fid, '%d\n', tmp_nodes(:)');
end
fclose(fid);

% Remove NaNs
fileText = fileread(filepath);
fid = fopen(filepath, 'w');
fwrite(fid, strrep(fileText, 'NaN', '')); % case-sensitive
fclose(fid);
end










