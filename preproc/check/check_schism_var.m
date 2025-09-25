function varData = check_schism_var(varData, varName)
% Clip SCHISM input variables to the valid range.

% Assign valid ranges
switch lower(varName)
    case {'elev', 'ssh'}
        vmin = -98; vmax = 98;
    case 'temp'
        vmin = -2; vmax = 40;
    case 'salt'
        vmin = 0; vmax = 42;
    case {'uvel', 'vvel'}
        vmin = -98; vmax = 98;
    otherwise  % for other tracers (ICM, CoSiNE, etc)
        vmin = 0; vmax = 999;
end
warning on

% # of NaN values
num_of_nans = numel(find(isnan(varData(:))));

% Check invalid values
if min(varData(:))<vmin || max(varData(:))>vmax
    warning([varName,' range: [',num2str(min(varData(:)), '%.3f'), ',', num2str(max(varData(:)), '%.3f') ,']; # of NaNs = ', num2str(num_of_nans)])
    disp(['• ', varName,' is clipped to [', num2str(vmin), ', ', num2str(vmax), ']'])
    varData = min(vmax, max(vmin, varData));
else
    disp(['• ', varName,' range: [',num2str(min(varData(:)), '%.3f'), ',', num2str(max(varData(:)), '%.3f') ,']; # of NaNs = ', num2str(num_of_nans)])
end

end


