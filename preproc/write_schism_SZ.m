function Mobj = write_schism_SZ(Mobj, tuning_coefs, zLayers, sigNums, hc)
% Need improvements Later
if isfield(Mobj, 'vcor_type')
    error(['vertical grid has already been designated as ',Mobj.vcor_type,' before!'])
else
    Mobj.vcor_type = 'SZ';
end
if nargin < 2
    tuning_coefs = [0.7 5];
    zLayers = [100 105 110 120 140 190 240 290 340 390 470 570 770 1000 1400 1800 2300 5000];
    sigNums = 20;
    hc = 30;
end
if nargin < 3
    zLayers = [100 105 110 120 140 190 240 290 340 390 470 570 770 1000 1400 1800 2300 5000];
    sigNums = 20;
    hc = 30;
end
if nargin < 4
    sigNums = 20;
    hc = 30;
end
if nargin < 5
    hc = 30;
end

zNums = length(zLayers);
hs = min(zLayers);

theta_b = tuning_coefs(1);
theta_f = tuning_coefs(2);
stretchFcn = @(x) (1-theta_b)*sinh(theta_f*x)/sinh(theta_f) + theta_b*(tanh(theta_f*(x+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5);

zLevs = 1:sigNums;
sigma = (zLevs-1.0)/(1.0-sigNums);
sigma_stret = stretchFcn(sigma);

zcors = sort(-zLayers);
scors = sort(sigma_stret);
nLayers = zNums+sigNums;

fileName = [Mobj.aimpath, 'vgrid_SZ.in'];
fid = fopen(fileName,'wt');
fprintf(fid, '%d %s \n', 2, '!ivcor (2: SZ; 1: VQS)');
fprintf(fid, '%d %d %d. %s \n', nLayers, zNums, hs , '!nvrt, kz (# of Z-levels); h_s (transition depth between S and Z)');
fprintf(fid, '%s \n', 'Z levels');
for ii = 1:zNums
    if ii == zNums
        fprintf(fid, '%d    %d. %s \n', ii, -hs, '!last z-coord. must be -h_s');
    else
        fprintf(fid, '%d    %d. \n', ii, zcors(ii));
    end
end
fprintf(fid, '%s \n', 'S levels');
fprintf(fid, '%d. %2.1f %d. %s\n', hc, theta_b, theta_f, '!h_c, theta_b, theta_f');
fprintf(fid, '   %d   %d. %s\n', zNums, -1, '!first S-level (S-coordinate must be -1)');
for jj = 1:sigNums
    if jj == sigNums
        fprintf(fid, '   %d    %d. %s \n', jj+1, 0, ' !last S-coodin. must be 0');
    elseif jj == 1
    else
        fprintf(fid, '   %d    %9.6f \n', jj, scors(jj));
    end
end
fclose(fid);
depLayers = ctranspose(Mobj.depth.*sigma_stret);
depLayers(depLayers<-hs) = nan;
depLayers = [depLayers; nan(zNums, size(depLayers,2))];


tmp = sum(~isnan(depLayers));
for ii = 1:size(depLayers,2)
    progressbar(ii/size(depLayers,2))
    if tmp(ii)~=sigNums
        ind = tmp(ii)+1:tmp(ii)+zNums;
        deps_cali = -zLayers;
        deps_cali(deps_cali<-Mobj.depth(ii)) = nan;
        deps_cali(find(isnan(deps_cali), 1)) = -Mobj.depth(ii);
        depLayers(ind, ii) = deps_cali;
    end
end
Mobj.depLayers = depLayers;
end











