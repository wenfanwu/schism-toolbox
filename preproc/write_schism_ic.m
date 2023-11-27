function write_schism_ic(Mobj, prefix_name, varData)
% Write the *ic files for SCHISM
% 
%% Syntax
% write_schism_ic(Mobj, prefix_name, varData)
%
%% Description 
% write_schism_ic(Mobj, prefix_name, varData) writes the [prefix_name].ic
% file using 'varData'.
% 
%% Example
% ssh = rand(Mobj.nNodes ,1);
% write_schism_ic(Mobj, 'elev', ssh)
%
% write_schism_ic(Mobj, 'salt', 30)
%
%% Input Arguments
% Mobj --- the mesh object
% prefix_name --- prefix name for the ic file.
% varData --- variable vector or scalar. If the input is a scalar, it will
% be expanded to a Mobj.nNodes*1 vector automatically.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 9 Nov. 2021. 
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: fprintf and write_schism_hotstart

%% Parse inputs
if numel(varData) == 1
    varData = varData*ones(1, Mobj.nNodes);
end
elemType = 3;
head_line = datestr(now);

%% Check
switch prefix_name
    case 'temp'
        if min(varData)<-2 || max(varData)>40
            warning on
            warning('Invalid temp. values have been removed!')
            varData = max(-2, varData);
            varData = min(40, varData);
        end
    case 'salt'
        if min(varData)<-2 || max(varData)>40
            warning on
            warning('Invalid salt. values have been removed!')
            varData = max(9, varData);
            varData = min(42, varData);
        end
end
%% Begin to Write
fileName = fullfile(Mobj.aimpath, [prefix_name, '.ic']);
fid = fopen(fileName,'wt');
fprintf(fid, [head_line, '\n']);
fprintf(fid, [num2str(Mobj.nElems),' ',num2str(Mobj.nNodes), '\n']);
for iNode = 1:Mobj.nNodes
    fprintf(fid, '%d  %14.6f  %14.6f  %13.7e\n', iNode, Mobj.lon(iNode), Mobj.lat(iNode), varData(iNode));
end
for iElem = 1:Mobj.nElems
    fprintf(fid, '%d %d %d %d %d\n',iElem, elemType, Mobj.tri(iElem,1), Mobj.tri(iElem,2), Mobj.tri(iElem,3));
end
fclose(fid);

disp([prefix_name, '.ic has been created successfully!'])

end