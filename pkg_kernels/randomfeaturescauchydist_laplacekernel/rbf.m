
function [K] = rbf(X1, X2, sigma)
K = X1 * X2';
X1_row_sq = sum(X1.^2, 2) / 2;
X2_row_sq = sum(X2.^2, 2) / 2;
K = bsxfun(@minus, K, X1_row_sq);
K = bsxfun(@minus, K, X2_row_sq');
K = K / (sigma^2);
K = exp(K);
