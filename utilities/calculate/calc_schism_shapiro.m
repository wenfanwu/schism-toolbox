function shapiro_val = calc_schism_shapiro(Mobj, slope_list, disp_flag, hdif_max)

if nargin < 2
    slope_list = [0.015 0.025];
    disp_flag = 0;
    hdif_max = 0.5;
end
if nargin < 3
    disp_flag = 0;
    hdif_max = 0.5;
end
if nargin < 4
    hdif_max = 0.5;
end
slope_min = slope_list(1);
slope_ref = slope_list(2);

gam = Mobj.slope;
shapiro_val = hdif_max*tanh(2*gam/slope_ref);
shapiro_val(shapiro_val<slope_ref) = 1.e-3;
if disp_flag == 1
    disp_schism_var(Mobj, shapiro_val, 1)
    colormap(flipud(hot(35)))
end

end