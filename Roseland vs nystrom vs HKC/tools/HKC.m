function [U] = HKC(data, Dim, ref, denoise, sig)

dist = pdist2(data, ref.set).^2;
if nargin < 5
    sig = quantile(max(dist, [], 2), .8) / 20;
end

W_ref = exp( - dist / sig );

if denoise
    % do a truncation
    trunc = quantile(W_ref, .8, 2);
    W_ref(W_ref <= trunc) = 0;
end

W_r = normr(W_ref);
W = W_r' * W_r;
[u_Haddad, ~] = eigs(W, Dim+1);
U = W_r * u_Haddad;
    
end