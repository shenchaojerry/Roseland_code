function [U, S] = roseland(data, Dim, ref_set, denoise)

[n, ~] = size(data) ;
% form affinity matrix wrt the ref set
affinity_ext = pdist2(data, ref_set).^2;
% affinity matrix W

%sig = quantile(max(affinity_ext, [], 2), .8) / 30;
sig = median(median(affinity_ext, 2)) / 5;

if denoise
    W_ref = exp( - 2 * affinity_ext / sig );
else
    W_ref = exp( - affinity_ext / sig );
end


W_ref = sparse(W_ref);

% make W row stochastic
% form sparse D = D^{-.5}
D = W_ref * sum(W_ref, 1)';

D = D.^(-.5);
D = sparse(1:n, 1:n, D, n, n);
W_ref = D * W_ref;


% SVD
[U, S, ~] = svds(W_ref, Dim+1);
U = D * U;
S = diag(S);
S = S .^ 2;
end