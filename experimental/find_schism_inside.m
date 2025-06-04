function is_in = find_schism_inside(Mobj, xq, yq)
% Return the index of points located within the model domain.

loop_nodes = find_loop_nodes(Mobj.tri, Mobj.edg);
idx = loop_nodes(:,1);  idx = [idx(:); idx(1)]; % index of the outer loop
is_in = inpolygon(xq(:), yq(:), Mobj.lon(idx), Mobj.lat(idx));

end