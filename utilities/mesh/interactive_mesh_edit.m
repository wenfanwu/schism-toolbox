function [x, y, tri] = interactive_mesh_edit(x0, y0, tri0)
% Interactive mesh editing tool for unstructured grids
%
%% Syntax
% [x, y, tri] = interactive_mesh_edit(x0, y0, tri0)
%
%% Description
% [x, y, tri] = interactive_mesh_edit(x0, y0, tri0) edits the unstructured grid.
%
%% Examples 
% 
%
%% Input Arguments
% x0 - x-axis coordinates; numeric
%       the x-axis coordinates defined at node centers.
% y0 - y-axis coordinates; numeric
%       the y-axis coordinates defined at node centers.
% tri0 - connectivity table; numeric
%       the connectivity table (N*4 or N*4) for elements.
%
%% Output Arguments
% x - x-axis coordinates; numeric
%       the updated x-axis coordinates defined at node centers.
% y - y-axis coordinates; numeric
%       the updated y-axis coordinates defined at node centers.
% tri - connectivity table; numeric
%       the updated connectivity table (N*4) for elements.
%
%% Tips
% Three modes are available:
%       1) MOVE Mode (shortcut key:1): hold down the mouse left to drag the selected node
%       2) ADD Mode (shortcut key:2): double-click the mouse left to add a new node
%       3) DELETE Mode (shortcut key:3): double-click the mouse left to delete the selected node
%
% It's recommened to activate a basemap first before editing, the basemap
% can be DEM or coastline map.
%
%% Notes
% This function was created with the help of ChatGPT. It is a lightweight
% unstructured grid optimization tool designed to manually refine some
% details on an existing grid.
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025. 
% Last Updated on 18 Mar 2025. 
% Email: wwu@vims.edu
% 
% See also: 

%% Parse inputs
disp([newline, ' MOVE Mode: hold down the mouse left to drag the selected node', newline, ...
        ' ADD Mode: double-click the mouse left to add a new node', newline, ...
        ' DELETE Mode: double-click the mouse left to delete the selected node'])
x = x0;
y = y0;
if size(tri0, 2) == 3
    tri = [tri0, nan(size(tri0,1),1)];
else
    tri = tri0;
end
mode = 'move';  
selected_node = [];
lastClickTime = 0;  
%% Begin to edit (basemap is recommented before editing)
fig = gcf;
set(fig, 'Color', 'w', 'Position', [100, 100, 800, 600]);
hold on;
box on;
axis image;
xlabel('X'); ylabel('Y');
xlim([min(x) max(x)])
ylim([min(y) max(y)])

% plot the unstructured grid
plt_elements = patch('Faces', tri(:, 1:end), 'Vertices', [x(:), y(:)], ...
    'FaceColor', 'none', 'EdgeColor', 'b');

% show all nodes (default)
plt_nodes = scatter([], [], 1, 'r', 'filled');

% UI button - select mode
move_btn = uicontrol('Style', 'pushbutton', 'String', 'Move', 'Position', [20, 520, 100, 30], 'Callback', @setMove);
add_btn = uicontrol('Style', 'pushbutton', 'String', 'Add', 'Position', [20, 480, 100, 30], 'Callback', @setAdd);
delete_btn = uicontrol('Style', 'pushbutton', 'String', 'Delete', 'Position', [20, 440, 100, 30], 'Callback', @setDelete);

% monitor keyboard key
set(fig, 'WindowKeyPressFcn', @keyPress);

% interactive tools
set(fig, 'WindowButtonDownFcn', @mouseDown);
set(fig, 'WindowButtonMotionFcn', @mouseMove);
set(fig, 'WindowButtonUpFcn', @mouseUp);

% "move" mode is activated by default
setMove();

    % shortcut keys to select mode
    function keyPress(~, event)
        switch event.Key
            case '1'
                setMove();
            case '2'
                setAdd();
            case '3'
                setDelete();
        end
    end

    function setMove(~, ~)
        mode = 'move';
        updateButtonColors();
    end

    function setAdd(~, ~)
        mode = 'add';
        updateButtonColors();
    end

    function setDelete(~, ~)
        mode = 'delete';
        updateButtonColors();
    end

    function updateButtonColors()
        move_btn.BackgroundColor = [0.94 0.94 0.94];  
        add_btn.BackgroundColor = [0.94 0.94 0.94];
        delete_btn.BackgroundColor = [0.94 0.94 0.94];

        switch mode
            case 'move'
                move_btn.BackgroundColor = [0.7 0.9 1];  
            case 'add'
                add_btn.BackgroundColor = [0.7 1 0.7];  
            case 'delete'
                delete_btn.BackgroundColor = [1 0.7 0.7];  
        end
    end

    function mouseDown(~, ~)
        cp = get(gca, 'CurrentPoint');
        pos = cp(1, 1:2);
        
        % check if the mouse position is within the axis limits
        x_limits = xlim(gca);  % get x-axis limits
        y_limits = ylim(gca);  % get y-axis limits
        
        if pos(1) < x_limits(1) || pos(1) > x_limits(2) || pos(2) < y_limits(1) || pos(2) > y_limits(2)
            % If the mouse is outside the axis limits, disable dragging
            return;
        end

        % calculate the closet points
        dists = sqrt((x - pos(1)).^2 + (y - pos(2)).^2);
        [minDist, idx] = min(dists);
    
        clickTime = now; %#ok<TNOW1>
    
        if strcmp(mode, 'move')
            if minDist < 0.1  
                selected_node = idx;
                updateNodeColors();
            else
                selected_node = [];
                updateNodeColors();
            end
        elseif strcmp(mode, 'add')
            if (clickTime - lastClickTime) * 24 * 60 * 60 < 0.3  
                addNode(pos);
            end
        elseif strcmp(mode, 'delete')
            if minDist < 0.1  
                % single-click to select (turns green)
                selected_node = idx;
                updateNodeColors();
                
                % double-click to delete node
                if (clickTime - lastClickTime) * 24 * 60 * 60 < 0.3  
                    deleteNode(idx);
                    selected_node = [];
                    updateNodeColors();
                end
            end
        end
        lastClickTime = clickTime;
    end

    function mouseMove(~, ~)
        if strcmp(mode, 'move') && ~isempty(selected_node)
            cp = get(gca, 'CurrentPoint');
            x(selected_node) = cp(1, 1);
            y(selected_node) = cp(1, 2);
            updatePlot();
        end
    end

    function mouseUp(~, ~)
        selected_node = [];
    end

    function deleteNode(idx)
        % delete the node "idx"
        x(idx) = [];
        y(idx) = [];
        
        % removes all triangles/quadrangles that contain this node
        tri(any(tri == idx, 2), :) = [];
       
        % update index
        tri(tri > idx) = tri(tri > idx) - 1;
        
        % remove suspended nodes (isolated nodes that does not form any elements)
        wild_nodes = setdiff(1:numel(x), unique(tri(:)));
        if ~isempty(wild_nodes)
            for ii = 1:numel(wild_nodes)
                idx = wild_nodes(ii);
                x(idx) = []; y(idx) = [];
                tri(any(tri == idx, 2), :) = [];
                tri(tri > idx) = tri(tri > idx) - 1;
            end
        end
        
        updatePlot();
    end
    
    function updateNodeColors()
        delete(plt_nodes); % delete the blank plt_nodes
        if ~isempty(selected_node) && selected_node <= length(x)
            plt_nodes = scatter(x(selected_node), y(selected_node), 40, 'g', 'filled'); % greenc color means the selected node
        else
            plt_nodes = scatter([], [], 5, 'r', 'filled');
        end
    end

    function addNode(pos)
        boundary_edges = findBoundaryEdges(tri);
        boundary_edges = boundary_edges(~any(isnan(boundary_edges), 2), :);
        boundary_points = unique(boundary_edges(:));
    
        if isempty(boundary_points) || length(boundary_points) < 2
            disp('Error: No boundary points found!');
            return;
        end
    
        % calculate the Euclidean distance from all boundary nodes ()
        dists = sqrt((x(boundary_points) - pos(1)).^2 + (y(boundary_points) - pos(2)).^2);
        [~, idx] = sort(dists);
        
        % selected the closet two nodes
        nearest_pts = boundary_points(idx(1:2));
    
        % check if the two nodes belong to the same boundary edge
        ind_e1 = find(boundary_edges(:,1)==nearest_pts(1) & boundary_edges(:,2)==nearest_pts(2));
        ind_e2 = find(boundary_edges(:,1)==nearest_pts(2) & boundary_edges(:,2)==nearest_pts(1));
    
        if ~isempty([ind_e1; ind_e2])
            % if so, add new nodes and elements
            x = [x; pos(1)];
            y = [y; pos(2)];
            new_idx = length(x);
            selected_node = new_idx;
            
            tri = [tri; new_idx, nearest_pts(1), nearest_pts(2), NaN];
        else
            % if not, directly connect the two nodes and add new elements
            ind_s1 = boundary_edges(:,1)==nearest_pts(1) | boundary_edges(:,2)==nearest_pts(1);
            ind_s2 = boundary_edges(:,1)==nearest_pts(2) | boundary_edges(:,2)==nearest_pts(2);

            s1 = boundary_edges(ind_s1,:); s2 = boundary_edges(ind_s2,:);
            tri = [tri; nearest_pts(1), nearest_pts(2), intersect(s1,s2), NaN];  % note the rotation direction
        end
    
        updatePlot();
    end

    function edges = findBoundaryEdges(tri)
        % calculate the edges of all elements
        edges = [tri(:, [1 2]); tri(:, [2 3]); tri(:, [3 1])];  
        if ~all(isnan(tri(:, 4)))
            edges = [edges; tri(:, [3 4]); tri(:, [4 1])];  
        end
        
        % count the # of occurrences of each edge
        edges = sort(edges, 2);
        [unique_edges, ~, ic] = unique(edges, 'rows');
        counts = accumarray(ic, 1);
        
        % find the edges that appears only once (outer boundary)
        edges = unique_edges(counts == 1, :);
    end

    function updatePlot()
        delete(plt_nodes);
        if ~isempty(selected_node) && selected_node <= length(x)
            plt_nodes = scatter(x(selected_node), y(selected_node), 40, 'g', 'filled'); % greenc color means the selected node
        else
            plt_nodes = scatter([],[], 5, 'r', 'filled');
        end
        plt_elements.Vertices = [x(:), y(:)];
        plt_elements.Faces = tri(:, 1:end);
        drawnow;
    end

uicontrol('Style', 'pushbutton', 'String', 'Close', 'Position', [20, 400, 100, 30], 'Callback', @closeFigure);

    function closeFigure(~, ~)
        assignin('base', 'x', x);
        assignin('base', 'y', y);
        assignin('base', 'tri', tri);
        uiresume(fig);
        close(fig);
    end

uiwait(fig);
end
