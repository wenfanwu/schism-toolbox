function com_array = cell2array(input_cell, s_dim)
% Convert the numeric vectors inside a cell to an array struct. 

if nargin == 1; s_dim = 0; end
input_cell = input_cell(:);

if s_dim == 0
    max_sz = size(input_cell{1});
    com_array = cell2mat(cellfun(@(x) reshape(x, [1 max_sz]), input_cell, 'UniformOutput', false));
else
    com_array = cell2mat(cellfun(@(x) permute(x, [s_dim setdiff(1:3, s_dim)]), input_cell, 'UniformOutput', false));
    re_seq = nan(1,3);
    re_seq(s_dim) = 1;
    re_seq(isnan(re_seq)) = setdiff(1:3, 1);
    com_array = permute(com_array, re_seq);
end
end