function sect_info = def_schism_transect(Mobj, mtype, N)
% Define a transect on the model grid
%
%% Syntax
% sect_info = def_schism_transect(Mobj, mtype)
% sect_info = def_schism_transect(Mobj, mtype, N)
%
%% Description
% sect_info = def_schism_transect(Mobj) defines a transect on the map.
% sect_info = def_schism_transect(Mobj, mtype) specifies the method to define transects.
% sect_info = def_schism_transect(Mobj, mtype, N) specifies the # of transect points.
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the mesh object with vertical grid info added.
% mtype - method type; numeric
%       Flag specifying the method to define transects. Six options are
%       available (see below). Default: mtype = 1;
%       Mode options:
%           0 : Select node centers manually using datatip; return selected nodes.
%           1 : Draw a straight line; return nearest nodes to the transect points.
%           2 : Draw polyline segments; return nearest nodes to the transect points.
%          -1 : Same as option 1, but points are not restricted to node centers.
%          -2 : Same as option 2, but points are not restricted to node centers.
%          -3 : Freehand curve (advanced option).
% N - # of transect points; numeric
%       when mtype < 0: N specifies the # of transect points;
%       when mtype = 0: N determines transect orientation
%                                  (>0: latitude-descending; <0: latitude-ascending);
%       when mtype > 0: N specifies the approximate # of transect points
%                                  (depending on the grid resolution near the transect).
%
%% Output Arguments
% sect_Info - transect info; datastruct
%       the datastruct containing transect data.
%
%% Examples
% figure('Color', 'w')
% disp_schism_hgrid(Mobj, [1 0], 'EdgeAlpha', 0.05, 'LineWidth', 0.5);
% sect_info = def_schism_transect(Mobj, -1, 100);
%
%% Notes
% If it is necessary to define transects strictly along curved channels or
% thalweg, the freehand curve approach (mtype = -3) would be very useful.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 23 Sep 2025.
% Email: wwu@vims.edu
%
% See also: disp_schism_vgrid and read_schism_transect

%% Parse inputs
if nargin < 2; mtype = 1; end
if nargin < 3; if mtype == 0; N = 1; else; N = 100; end; end

%% Define transect
hold on
switch mtype
    case 0  % scattered points (points are on the node center)
        disp('Click node centers on the map (press ENTER to end)')
        datacursormode on; pause;
        h = findobj(gcf,'Type','datatip');
        if ~isempty(h)
            xy_data = cell2array(arrayfun(@(x) [h(x).X; h(x).Y], 1:length(h), 'UniformOutput', false)); close gcf
            x_pts = xy_data(:,1); y_pts =  xy_data(:,2); nps = numel(x_pts);
        end
        if sign(N) == 1  % latitude-descending order
            [y_pts, ind] = sort(y_pts, 'descend'); x_pts = x_pts(ind);
        elseif sign(N) < 1 % latitude-ascending order
            [y_pts, ind] = sort(y_pts, 'ascend'); x_pts = x_pts(ind);
        end
        idx_pts = arrayfun(@(x) geomin(Mobj.lon, Mobj.lat, x_pts(x), y_pts(x)), 1:nps);

    case 1 %  straight line (points are on the node center)
        disp('Draw a line on the map (press ENTER to end)')
        h = drawline; pause;
        x_ept = h.Position(:,1)'; y_ept = h.Position(:,2)';  % endpoints
        x_roi = interp1(1:2, x_ept, linspace(1,2,N)); y_roi = interp1(1:2, y_ept, linspace(1,2,N));  % refine the transect points
        idx_pts = unique(geomin(Mobj.lon, Mobj.lat, x_roi, y_roi), 'stable');      % find the nearest nodes
        x_pts = Mobj.lon(idx_pts); y_pts = Mobj.lat(idx_pts);
        close gcf

    case -1 % straight line (points are not on the node center)
        disp('Draw a line on the map (press ENTER to continue)')
        h = drawline; pause;
        x_ept = h.Position(:,1)'; y_ept = h.Position(:,2)';
        x_pts = interp1(1:2, x_ept, linspace(1,2,N)); y_pts = interp1(1:2, y_ept, linspace(1,2,N));
        [x_pts, y_pts] = fine_tuning(h, x_pts, y_pts);  % fine-tuning
        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        idx_pts = minfind(Mobj.depth, F(x_pts, y_pts));  % use vertical layers of the depth-closest nodes
        close gcf

    case 2  % broken line (points are on the node center)
        disp('Draw broken lines on the map (press ENTER to end)')
        h = drawpolyline;
        x_ept = h.Position(:,1)'; y_ept = h.Position(:,2)';
        x_roi = []; y_roi = [];
        L_segs = sqrt(diff(x_ept).^2 + diff(y_ept).^2);  % segment lengh
        Ns = round(N.*L_segs./(sum(L_segs)));  % distribute # of points 
        for s = 1:numel(x_ept)-1
            x_seg = interp1(1:2, x_ept(s:s+1), linspace(1,2,Ns(s)));
            y_seg = interp1(1:2, y_ept(s:s+1), linspace(1,2,Ns(s)));
            x_roi = [x_roi; x_seg(:)]; y_roi = [y_roi; y_seg(:)]; %#ok<AGROW>
        end
        [~, ia, ~] = unique(hypot(x_roi, y_roi), 'stable');  % remove duplicate points
        x_roi = x_roi(ia); y_roi = y_roi(ia);
        idx_pts = unique(geomin(Mobj.lon, Mobj.lat, x_roi, y_roi), 'stable');
        x_pts = Mobj.lon(idx_pts); y_pts = Mobj.lat(idx_pts);  % find the nearest nodes
        close gcf; clear h

    case -2  % broken lines (points are not on the node center)
        disp('Please draw broken lines on the map (press ENTER to continue)')
        h = drawpolyline;
        x_ept = h.Position(:,1)'; y_ept = h.Position(:,2)';
        x_roi = []; y_roi = [];
        L_segs = sqrt(diff(x_ept).^2 + diff(y_ept).^2);  % segment lengh
        Ns = round(N.*L_segs./(sum(L_segs)));  % distribute # of points
        for s = 1:numel(x_ept)-1
            x_seg = interp1(1:2, x_ept(s:s+1), linspace(1,2,Ns(s)));
            y_seg = interp1(1:2, y_ept(s:s+1), linspace(1,2,Ns(s)));
            x_roi = [x_roi; x_seg(:)]; y_roi = [y_roi; y_seg(:)]; %#ok<AGROW>
        end
        [~, ia, ~] = unique(hypot(x_roi, y_roi), 'stable');

        x_pts = x_roi(ia); y_pts = y_roi(ia);
        [x_pts, y_pts] = fine_tuning(h, x_pts, y_pts);  % fine-tuning
        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        idx_pts = minfind(Mobj.depth, F(x_pts, y_pts));
        close gcf; clear h

    case -3  % freehand curve (advanced option)
        disp('Draw a curve freehand along the target feature (press ENTER to continue)')
        h = drawfreehand; pause
        d = [0; cumsum(sqrt(sum(diff(h.Position).^2,2)))];  % arc length parameter

        % Resample to uniform intervals
        window = 3;  % smoothing window (can be changed)
        d_new = linspace(0, d(end), N*3);   % Ensure enough points for smoothing
        x_new = smooth(interp1(d, h.Position(:,1), d_new, 'pchip'), window*3);
        y_new = smooth(interp1(d, h.Position(:,2), d_new, 'pchip'), window*3);
        x_pts = x_new(1:3:end); y_pts = y_new(1:3:end);

        [x_pts, y_pts] = fine_tuning(h, x_pts, y_pts);  % fine-tuning
        F = scatteredInterpolant(Mobj.lon, Mobj.lat, Mobj.depth);
        idx_pts = minfind(Mobj.depth, F(x_pts, y_pts));
        close gcf

end
%% Along-transect distance
switch lower(Mobj.coord(1:3))
    case 'geo'
        dist_pts = distance(y_pts(1), x_pts(1), y_pts, x_pts, [6378137 0.0818191910428158]); % meters
    case 'car'
        dist_pts = hypot(x_pts-x_pts(1), y_pts-y_pts(1));
end

%% Along-tranect Tangent & Normal vectors
[tvec, nvec] = transect_vector(x_pts, y_pts, 'left');

%% Store the transect info
sect_info.lon = x_pts(:);
sect_info.lat = y_pts(:);
sect_info.dist = dist_pts(:);

% Add depth layers.
if isfield(Mobj, 'depLayers')
    sect_info.depth = Mobj.depLayers(:, idx_pts(:));
end
sect_info.tvec = tvec;  % tangent vector
sect_info.nvec = nvec;  % normal vector

end

function [x_new, y_new] = fine_tuning(h, x_raw, y_raw)

% Manual fine-tuning
delete(h);    % clean existed object
disp('Perform manual fine-tuning if necessary (press ENTER to end)')
h = drawpolyline('Position', [x_raw(:) y_raw(:)], 'InteractionsAllowed','reshape'); pause
x_new = h.Position(:,1); y_new = h.Position(:,2);

end
