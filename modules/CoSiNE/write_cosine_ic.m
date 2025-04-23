function write_cosine_ic(Mobj, InitCnd)
% Write the COS_hvar_*.ic files for the CoSiNE module
% 
%% Syntax
% write_cosine_ic(Mobj, InitCnd)
%
%% Description 
% write_cosine_ic(Mobj, InitCnd) writes the COS_hvar_*.ic files for the
% CoSiNE module.
%
%% Example
% write_cosine_ic(Mobj, InitCnd)
%
%% Input Arguments
% Mobj --- the mesh object
%
%% Output Arguments
% None
% 
%% Notes
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2022-05-29.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs 
varList = {'NO3','SiO4','NH4','S1'	'S2','Z1','Z2','DN','DSi','PO4','DOX','CO2','ALK'}; % DO NOT change the order

for iVar = 1:numel(varList)
    varName = varList{iVar};
    ind_var = find(strcmpi({InitCnd.Variable}, varName));
    if ~isempty(ind_var)
        varTmp = InitCnd(ind_var).Data(1,:);
        varTmp(varTmp<0) = 0;
        write_schism_ic(Mobj, ['COS_hvar_', num2str(iVar)], varTmp)
    end
end

end