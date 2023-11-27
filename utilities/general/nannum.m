function num_of_nans = nannum(varData)
num_of_nans = numel(find(isnan(varData(:))));
end