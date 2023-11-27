function DS = get_cosine_init_test(Mobj)
% Extract BGC data from the fake dataset

region = Mobj.region;
varList = {'no3', 'sio4', 'nh4', 's1', 's2', 'z1', 'z2', 'dn', 'dsi', 'po4', 'dox', 'co2', 'alk'}; 
default_list = [10, 0.5, 8, 300, 1000, 0.1, 0.1, 0.1, 0.1, 3, 3, 1, 1024];

test_data.lon = (region(1):0.1:region(2))';
test_data.lat = (region(3):0.1:region(4))';
test_data.depth = (1:5:max(Mobj.depth)+10)';
test_data.time = Mobj.time(1);

nLons = numel(test_data.lon);
nLats = numel(test_data.lat);
nDeps = numel(test_data.depth);
test_data.var = ones(nLons, nLats, nDeps);

for iVar = 1:numel(varList)
    varName = varList{iVar};
    default_var = default_list(iVar);
    D = test_data;
    D.var = default_var*D.var;
    
    DS.(varName) = D;
    clear D
end
end