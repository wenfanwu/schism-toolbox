function write_schism_vgrid(Mobj, version_num)
% Write the vgrid.in file for SCHISM.
% 
%% Syntax
% write_schism_vgrid(Mobj)
%
%% Description 
% write_schism_vgrid(Mobj) writes the vgrid.in file for SCHISM.
% 
%% Example
% Mobj = gen_schism_LSC2(Mobj, 3, [4 5 4 5], 0.35);
% write_schism_vgrid(Mobj);
%
% Mobj = gen_schism_SZ(Mobj, [10, 0.7, 5, 20], 20:2:150);
% write_schism_vgrid(Mobj);
%
%% Input Arguments
% Mobj - mesh object; datastruct
%       the datastruct used to store mesh info.
%
%% Output Arguments
% None
%
%% Author Info
% Created by Wenfan Wu, Virginia Institute of Marine Science in 2021.
% Last Updated on 23 Apr 2025.
% Email: wwu@vims.edu
% 
% See also: read_schism_vgrid and write_schism_hgrid

%% Parse inputs
% the format of vgrid.in has changed substantially since the version 5.10.
if nargin < 2
    version_num = 'v5.10';  % v5.10 or v5.9
else
    version_num = lower(version_num);
end

filepath = fullfile(Mobj.aimpath, 'vgrid.in');
switch upper(Mobj.vtype)
    case 'LSC2'
        if strcmp(version_num, 'v5.9') % old format of vgrid.in
            disp('vgird.in (LSC2) for the version v5.9 and below')
            fid = fopen(filepath,'wt');
            fprintf(fid, '%d !ivcor\n', 1);
            fprintf(fid, '%d !nvrt; below are node #, bottom level and sigma at each node (from bottom to surface)\n', Mobj.maxLev);
            fprintf(fid, ['%d %8d ', repmat('% 12.6f', 1, Mobj.maxLev), '\n'], [(1:Mobj.nNodes); Mobj.nLevs(:)'; Mobj.vgrids]);
            fclose(fid);
            % Remove NaNs
            fileText = fileread(filepath);
            fid = fopen(filepath, 'w');
            fwrite(fid, strrep(fileText, 'NaN', '')); % case-sensitive
            fclose(fid);
        end
        if strcmp(version_num, 'v5.10')   % new format of vgrid.in for the version v5.10 and above
            disp('vgird.in (LSC2) for the version v5.10 and above')
            vgrids = flipud(Mobj.vgrids); vgrids(isnan(vgrids)) = -9;
            fid = fopen(filepath,'w');
            fprintf(fid, '%d !ivcor\n', 1);
            fprintf(fid, '%d !nvrt (=Nz) \n', Mobj.maxLev);
            fprintf(fid,  [repmat('%10d', 1, Mobj.nNodes), '\n'], Mobj.maxLev+1-Mobj.nLevs);
            fprintf(fid,  ['%d ', repmat('%15.6f', 1, Mobj.nNodes), '\n'], [(1:Mobj.maxLev)', vgrids]');
            fclose(fid);
        end

    case {'SZ', 'S', 'Z'} 
        disp('vgird.in (SZ) for the version v5.10 and above')

        h_c = Mobj.s_conts(1);
        theta_b = Mobj.s_conts(2);
        theta_f = Mobj.s_conts(3);
        ks = Mobj.s_conts(4);
        h_s = min(abs(Mobj.zcors));
        kz = numel(Mobj.zcors);
        nvrt = kz+ks-1;
        
        fid = fopen(filepath,'w');
        fprintf(fid, '%d %s \n', 2, ' ! ivcor (1: LSC2; 2: SZ)');
        fprintf(fid, '%d %d %d. %s \n', nvrt, kz, h_s , '!nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)');
        fprintf(fid, '%s \n', 'Z levels');
        for ii = 1:kz
            if ii == kz
                fprintf(fid, '%d    %d. %s \n', ii, -h_s, '!last z-coord must match -h_s');
            else
                fprintf(fid, '%d    %d. \n', ii, Mobj.zcors(kz-ii+1));
            end
        end
        fprintf(fid, '%s \n', 'S levels');
        fprintf(fid, '%f %2.1f %d. %s\n', h_c, theta_b, theta_f, '!h_c, theta_b, theta_f');
        for jj = 1:ks
            switch jj
                case 1
                    fprintf(fid, '   %d   %d. %s\n', kz, -1, '!first S-level (S-coordinate must be -1)');
                case ks
                    fprintf(fid, '   %d    %d. %s \n', jj, 0, ' !last S-coord must be 0');
                otherwise
                    fprintf(fid, '   %d    %9.6f \n', jj, Mobj.sigma(ks-jj+1));
            end
        end
        fclose(fid);
end
disp(['vgrid.in (',Mobj.vtype,') has been created successfully!'])

end





