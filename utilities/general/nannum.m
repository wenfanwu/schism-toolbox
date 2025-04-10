function num_of_nans = nannum(varData)
% Calculate the # of NaNs

num_of_nans = numel(find(isnan(varData(:))));
end