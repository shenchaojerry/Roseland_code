function [U, S] = roseland(data, Dim, ref, denoise, sig)

[n, ~] = size(data) ;
% form affinity matrix wrt the ref set
affinity_ext = pdist2(data, ref.set).^2;
% affinity matrix W
if nargin < 5
    sig = quantile(max(affinity_ext, [], 2), .8) / 20;
end

W_ref = exp( - affinity_ext / sig );

if denoise
    % do a truncation
    trunc = quantile(W_ref, .8, 2);
    W_ref(W_ref <= trunc) = 0;
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