function [r_shift, c_shift] = VesselGenerator_CalculateMove(kernel)

size_kernel = size(kernel);

% Scores are random prob mult by kernel score. Kernel element with value 0
% means vessel cannot extend in that direction
score = kernel.* rand(size(kernel));
[~, match_ind] = max(score(:));

% convert to row col
[r, c] = ind2sub(size(kernel),match_ind);

% Compute what the shift is
r_shift = r - round(size_kernel(1)+1)/2;
c_shift = c - round(size_kernel(2)+1)/2;

% keyboard