function [U] = HKC(data, Dim, ref_set, denoise)

dist = pdist2(data, ref_set).^2;

%sig = quantile(max(dist, [], 2), .8) / 30;
sig = median(median(dist, 2)) / 5;


if denoise
    W_ref = exp( - 2 * dist / sig );
else
    W_ref = exp( - dist / sig );
end

W_r = normr(W_ref);
W = W_r' * W_r;
[u_Haddad, ~] = eigs(W, Dim+1);
U = W_r * u_Haddad;
    
end