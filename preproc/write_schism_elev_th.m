function write_schism_elev_th(Mobj, ElevMatrix)
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
% Created by Wenfan Wu, Ocean Univ. of China in 2023. 
% Last Updated on 2023-11-26.
% Email: wenfanwu@stu.ouc.edu.cn
% 
% See also: 

%% Parse inputs
time_steps = Mobj.time(1):seconds(Mobj.dt):Mobj.time(end);
time_ticks = seconds(time_steps-time_steps(1));

fileName = [Mobj.aimpath, 'elev.th'];
fid = fopen(fileName,'wt');
formatStr = ['%d.', repmat('% 11.5f', 1, Mobj.nObcNodes), '\n'];
for iTime = 1:length(time_ticks)
    fprintf(fid, formatStr, time_ticks(iTime), ElevMatrix(iTime,:));
end
fclose(fid);

disp('elev.th has been created successfully!')
end

% ElevMatrix (timeSteps*nObcNodes)
% timeSteps = ceil(Mobj.rnday*24*3600/Mobj.dt)+1;
% timeTicks = Mobj.dt*(0:timeSteps);