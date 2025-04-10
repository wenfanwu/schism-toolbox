function comArray = cell2array(inputCell, specDim)
% Convert the numeric vectors inside a cell to an array.
% struct. By default, it is assumed that all the components are of the same
% size, and then it will integrate them along a new dimension. For
% example, if the components are 3-D, the resulting intergrated array will
% be 4-D, and the length of the first dimension will be equal to the number
% of components. In addtiion, you can specify the intergrate dimension on
% your own, if one dimension of the components differs. For example, there
% are three matrixes in the cell struct, namely A(100*200*20), B(100*200*30),
% C(100*200*10), then, the resulting intergrated array will be D(100*200*60).
% ----------- Created by Wenfan Wu at 31 Dec. 2020, in Ocean University of China.
if nargin == 1
    specDim = 0;
end
inputCell = inputCell(:);
if specDim == 0
    matrixSize = size(inputCell{1});
    comArray = cell2mat(cellfun(@(x) reshape(x, [1 matrixSize]), inputCell, 'UniformOutput', false));
else
    comArray = cell2mat(cellfun(@(x) permute(x, [specDim setdiff(1:3, specDim)]), inputCell, 'UniformOutput', false));
    reSeq = nan(1,3);
    reSeq(specDim) = 1;
    reSeq(isnan(reSeq)) = setdiff(1:3, 1);
    comArray = permute(comArray, reSeq);
end
end