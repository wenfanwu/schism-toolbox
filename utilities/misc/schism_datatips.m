function output_txt = schism_datatips(obj,event_obj)
% Display data cursor position in a data tip
% obj          Currently not used
% event_obj    Handle to event object
% output_txt   Data tip text, returned as a character vector or a cell array of character vectors

pos = event_obj.Position;

nNodes = size(event_obj.Target.Vertices,1);
nElems = size(event_obj.Target.XData,2);
Cs = event_obj.Target.FaceVertexCData;

Xs = event_obj.Target.Vertices(:,1);
Ys = event_obj.Target.Vertices(:,2);

if nElems==1
    Cs = Cs(1:3:end); 
    Xs = (Xs(1:3:end)+Xs(2:3:end))/2; 
    Ys = (Ys(1:3:end)+Ys(2:3:end))/2; 
    
    Index = geomin(Xs(:), Ys(:), pos(1), pos(2));  % need to be changed

elseif numel(Cs)==nNodes
    Index = geomin(Xs(:), Ys(:), pos(1), pos(2));  

elseif numel(Cs)==nElems
    tri = event_obj.Target.Faces;
    i34 = sum(diff(sort(tri, 2), 1, 2) ~= 0, 2)+1;

    Xcs = nan(nElems, 1); Xcs(i34==3) = mean(Xs(tri(i34==3, 1:3)), 2);Xcs(i34==4) = mean(Xs(tri(i34==4, 1:4)), 2);
    Ycs = nan(nElems, 1); Ycs(i34==3) = mean(Ys(tri(i34==3, 1:3)), 2);Ycs(i34==4) = mean(Ys(tri(i34==4, 1:4)), 2);

    Index = geomin(Xcs(:), Ycs(:), pos(1), pos(2));
end

%********* Define the content of the data tip here *********%

% Display the x and y values:
output_txt = {['Index',formatValue(Index,event_obj)], ...
    ['X',formatValue(pos(1),event_obj)],...
    ['Y',formatValue(pos(2),event_obj)]};

% Display the C values:
C = Cs(Index);
output_txt{4} = ['C',formatValue(C,event_obj)];

%***********************************************************%

function formattedValue = formatValue(value, event_obj)
% Custom formatting function for displaying data tips
% Displays integers as whole numbers and non-integers with up to 4 decimal places.

% If you do not want TeX formatting in the data tip, uncomment the line below:
% event_obj.Interpreter = 'none';

% Check if TeX formatting is enabled
if strcmpi(event_obj.Interpreter, 'tex')
    % Define TeX formatting for color and bold text
    valueFormat = ' \color[rgb]{0 0.6 1}\bf';
    removeValueFormat = '\color[rgb]{.25 .25 .25}\rm';
else
    % Define non-TeX formatting
    valueFormat = ': ';
    removeValueFormat = '';
end

% Check if the value is an integer (considering floating-point precision)
if mod(value, 1) == 0
    % Format as an integer if the value is a whole number
    formattedValue = [valueFormat num2str(value, '%d') removeValueFormat];
else
    % Format as a floating-point number with 4 decimal places
    formattedValue = [valueFormat num2str(value, '%.4f') removeValueFormat];
end


