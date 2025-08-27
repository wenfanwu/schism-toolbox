function SS = def_schism_source(Mobj, ss_flags, load_flag, disp_flag)
% Define the source/sink elements on the map.
% 
%% Syntax
% SS = def_schism_source(Mobj, ss_flags)
% SS = def_schism_source(Mobj, ss_flags, load_flag)
% SS = def_schism_source(Mobj, ss_flags, load_flag, disp_flag)
%
%% Description 
% SS = def_schism_source(Mobj, ss_flags) defines source/sink elements on the map.
% SS = def_schism_source(Mobj, ss_flags, load_flag) loads or re-build the source/sink elements
% SS = def_schism_source(Mobj, ss_flags, load_flag, disp_flag) shows the results or not.
%
%% Example
% SS = def_schism_source(Mobj, [1 0], 1, 'load');
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store the mesh info.
% ss_flags - source/sink flags; numeric
%       a two-element vector, the first of which decides whether the source
%       is defined, while the latter is for the sink. Default: ss_flags = [1 0]; 
% load_flag - load flag; char
%       the flag used to rebuild, load, or add source/sink elements on the
%       map ('rebuild', 'load' or 'add'). Default: load_flag = 'rebuild'.
%       This function will save the defined source/sink as a MAT file, so
%       one can load or modify the old file without re-defining if there is
%       a "source_sink.mat" avaialble.  
% disp_flag - display flag; char
%       the flag used to show the results or not. (off/on). 
%       Default: disp_flag = 'off';
%
%% Output Arguments
% SS - source sink data; datastruct
%       the datastruct used to store source/sink data.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 6 May 2025.
% Email: wwu@vims.edu
% 
% See also: add_river_inputs and match_rivers

%% Parse inputs
if nargin < 2; ss_flags = [1 0]; end
if nargin < 3; load_flag = 'rebuild'; end  % add/rebuild/load
if nargin < 4; disp_flag = 'off'; end

%% Select
filepath = fullfile(Mobj.aimpath, 'source_sink.mat');

if strcmpi(load_flag, 'rebuild') || strcmpi(load_flag, 'add')
    if ss_flags(1) == 1
        title_str = 'Use datatips to select the source points, press ENTER to stop';
        dataTips = select_datatips(Mobj, title_str);
        [source.lonc, source.latc, source.elems] = get_datatips(Mobj, dataTips);
    else
        source.lonc = []; source.latc = []; source.elems = [];
    end
    if ss_flags(2) == 1
        title_str = 'please select sink points...';
        dataTips = select_datatips(Mobj, title_str);
        [sink.lonc, sink.latc, sink.elems] = get_datatips(Mobj, dataTips);
    else
        sink.lonc = []; sink.latc = []; sink.elems = [];
    end
    if strcmpi(load_flag, 'rebuild') % rebuild mode
        SS.source = source;
        SS.sink = sink;
    else  % add mode
        load(filepath) %#ok<*LOAD>
        SS.source.lonc = [SS.source.lonc; source.lonc]; %#ok<NODEF>
        SS.source.latc = [SS.source.latc; source.latc];
        SS.source.elems = [SS.source.elems; source.elems];
        SS.sink.lonc = [SS.sink.lonc; sink.lonc];
        SS.sink.latc = [SS.sink.latc; sink.latc];
        SS.sink.elems = [SS.sink.elems; sink.elems];
    end
    save(filepath, 'SS')
elseif strcmpi(load_flag, 'load') % load mode
    load(filepath)
else
    error('incorrect defining mode (valid modes: rebuild/load/add) !')
end
%% Display
if strcmpi(disp_flag, 'on')
    figure('Color', 'w')
    disp_schism_hgrid(Mobj, [1 0])
    hold on
    plot_schism_bnds(Mobj)
    colormap(white)
    colorbar off
    hold on
    scatter(SS.source.lonc, SS.source.latc, 50,'filled', 'r')
    scatter(SS.sink.lonc, SS.sink.latc, 50,'filled', 'b')
    title('source (red) and sink (blue) points')
    axis tight
    auto_center
end
end

function dataTips = select_datatips(Mobj, title_str)
% This function aims to select datatips on the map

figure('Color', 'w');
disp_schism_hgrid(Mobj, [0 0])
hold on
plot_schism_bnds(Mobj)
box on;
xlabel('Longitude (°E)', 'FontWeight','bold')
ylabel('Latitude (°N)', 'FontWeight','bold')
title(title_str)
hold on
scatter(Mobj.lonc, Mobj.latc,3,'filled', 'k')
colormap(white)
colorbar off
datacursormode
pause
dataTips = findobj(gcf,'Type','datatip');

end

function [lonc, latc, elems] = get_datatips(Mobj, dataTips)
% This function aims to get the coordinates from the datatips

if ~isempty(dataTips)
    coords = cell2array(arrayfun(@(x) [dataTips(x).X; dataTips(x).Y], 1:length(dataTips), 'UniformOutput', false));
    close all
    N = size(coords,1);
    lonc = coords(:,1);
    latc =  coords(:,2);
    elems = arrayfun(@(x) geomin(Mobj.lonc, Mobj.latc, lonc(x), latc(x)), 1:N);
    lonc = lonc(:);
    latc = latc(:);
    elems = elems(:);
else
    lonc = []; latc = []; elems = [];
end

end
