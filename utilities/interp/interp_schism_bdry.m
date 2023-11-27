function BdryCnd = interp_schism_bdry(Mobj, DS, varList)
% Interp the boundary data onto SCHISM vertical layers
%
%
%% Syntax
%
%
%% Description
%
%
%% Examples
%
%
%% Input Arguments
%
%
%% Output Arguments
%
%
%% Notes
%
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2022.
% Last Updated on 2022-10-25.
% Email: wenfanwu@stu.ouc.edu.cn
%
% See also:

%% Parse inputs
if nargin < 3
    varList = fieldnames(DS);
end
nVars = numel(varList);

%% Interpolate (vertical)
for iVar = 1:nVars
    varName = varList{iVar};
    D = DS.(varName);
    model_time = dateshift(Mobj.time, 'start', D.unit_time);  % expand the time
    [~, ind_exp] = find(model_time(:)==D.time(:)');
    if numel(ind_exp)~=numel(model_time)
        error('the time range of bdry inputs should cover the model time!')
    end
    depRaw = abs(D.depth);
    depRaw = depRaw-min(depRaw);
    depBnd = abs(Mobj.depLayers(:,D.ind));
    nTimes = numel(D.time);
    nNodes_obc = numel(D.ind);

    if strcmp(varName, 'ssh')
        varBnd = D.var;
    else
        varBnd = zeros(nNodes_obc, Mobj.maxLev, nTimes);
        for iNode = 1:nNodes_obc
            varRaw = squeeze(D.var(iNode,:,:));
            varRaw = fillmissing(varRaw, 'previous', 1);  % fill the potential missing values at deep layers
            depNew = depBnd(:,iNode);
            depNew = fillmissing(depNew, 'previous');
            varBnd(iNode,:,:) = multi_interp1(depRaw, varRaw, depNew, 1);
        end
    end

    % Kill nan values on each open boundaries to avoid wild values from interpolation
    obc_inds = [0, cumsum(Mobj.obc_lens)];
    for ii = 1:Mobj.obc_counts  
        ind_sect = obc_inds(ii)+1:obc_inds(ii+1);
        varSect = varBnd(ind_sect,:,:);
        varBnd(ind_sect,:,:) = kiil_nans(varSect);
    end
    if numel(find(isnan(varBnd(:)))) ~= 0
        error('there are still NaN values in the matrix')
    end

    varExp = varBnd(:,:, ind_exp);  % expand on the time dimension
    BdryCnd.(varName) = varExp;
end

end

function varBnd = kiil_nans(varBnd)
% Kill the NaNs values

if numel(find(isnan(varBnd(:)))) ~= 0
    warning on
    varBnd = fillmissing(varBnd, 'nearest', 1);
    varBnd = fillmissing(varBnd, 'nearest', 2);
    varBnd = fillmissing(varBnd, 'nearest', 3);
%     varBnd = fillmissing(varBnd, 'previous',1,'EndValues','previous');
%     varBnd = fillmissing(varBnd, 'previous',2,'EndValues','previous');
%     varBnd = fillmissing(varBnd, 'previous',3,'EndValues','previous');
    varBnd = fillmissing(varBnd, 'constant', mean(varBnd(:), 'omitnan'));
end

end















