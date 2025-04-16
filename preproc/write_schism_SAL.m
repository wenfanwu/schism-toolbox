function write_schism_SAL(Mobj, SAL)
% Write loadtide_[FREQ].gr3 files for SCHISM
%
%% Syntax
% write_schism_SAL(Mobj, SAL)
%
%% Description
% write_schism_SAL(Mobj, SAL)
%
%% Examples
% tideList = {'S2','M2','N2','K2', 'K1','P1','O1','Q1'};
% SAL = get_fes2014_SAL(Mobj, tideList);
% write_schism_SAL(Mobj, SAL)
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
% SAL -  SAL data; datastruct
%       the datastruct containing SAL data (sal_amp, sal_pha).
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2025.
% Last Updated on 15 Apr 2025.
% Email: wwu@vims.edu
%
% See also: get_fes2014_SAL

%% Parse inputs
tide_list = SAL.tide_list; 
nTides = numel(tide_list);

%% Check inputs
if any(isnan(SAL.sal_amp(:))); error('NaNs were found in "sal_amp"'); end
if any(isnan(SAL.sal_pha(:))); error('NaNs were found in "sal_pha"'); end

%% Begin to write
for iTide = 1:nTides
    filename = ['loadtide_' upper(tide_list{iTide}) '.gr3'];
    filepath = fullfile(Mobj.aimpath, filename);

    fid = fopen(filepath, 'w');
    fprintf(fid,'%s\n', upper(tide_list{iTide}));

    amp = SAL.sal_amp(:, iTide);
    pha = SAL.sal_pha(:, iTide);

    % Node and Elem
    fprintf(fid,'%d %d\n', Mobj.nElems, Mobj.nNodes);
    fprintf(fid,'%d %f %f %f %f\n', [(1:Mobj.nNodes)' Mobj.lon(:) Mobj.lat(:) amp(:) pha(:)]');
    fprintf(fid,'%d %d %d %d %d %d\n', [(1:Mobj.nElems)', Mobj.i34(:), Mobj.tri]');  % include NaN values
    fclose(fid);

    % Remove NaNs
    fileText = fileread(filepath);
    fid = fopen(filepath, 'w');
    fwrite(fid, strrep(fileText, 'NaN', '')); % case-sensitive
    fclose(fid);

    disp([filename, ' has been created successfully!'])
end
fclose all;

end