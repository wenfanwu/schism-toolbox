function n_dims = dimnum(vd)
% Calculate the # of valid dimensions

if isscalar(vd)
    n_dims = 0;
else
    n_dims = numel(find(size(vd)~=1));
end
end