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
% Mobj --- the mesh object, with vertical information loaded. Mobj must be 
% outputed from 'gen_schism_LSC2' or 'gen_schism_SZ' function first.
%
%% Output Arguments
% None
%
%% Notes
%  sz coordinate is not fully supoorted so far
%
%% Author Info
% Created by Wenfan Wu, Ocean Univ. of China in 2021. 
% Last Updated on 2021-11-09
% Email: wenfanwu@stu.ouc.edu.cn
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
            for iNode = 1:Mobj.nNodes
                nLevs = max(Mobj.nLevs)+1-Mobj.nLevs(iNode);
                sigLevs = rot90(Mobj.vgrids(1:Mobj.nLevs(iNode), iNode), 3);
                format_str = ['%d %8d ', repmat('% 12.6f', 1, length(sigLevs)), '\n'];
                fprintf(fid, format_str, iNode, nLevs, sigLevs);
            end
            fclose(fid);
        end
        if strcmp(version_num, 'v5.10')   % new format of vgrid.in for the version v5.10 and above
            disp('vgird.in (LSC2) for the version v5.10 and above')
            fid = fopen(filepath,'w');
            fprintf(fid, '%d !ivcor\n', 1);
            fprintf(fid, '%d !nvrt (=Nz) \n', Mobj.maxLev);
            vgrids = flipud(Mobj.vgrids);
            vgrids(isnan(vgrids)) = -9;
            ind_btm = Mobj.maxLev+1-Mobj.nLevs;
            format_str = [repmat('%10d', 1, Mobj.nNodes), '\n'];
            fprintf(fid, format_str, ind_btm);
            for iLev = 1:Mobj.maxLev
                vgrid_lev = vgrids(iLev,:);
                format_str = ['%d ', repmat('%15.6f', 1, Mobj.nNodes), '\n'];
                fprintf(fid, format_str, iLev, vgrid_lev);
            end
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





