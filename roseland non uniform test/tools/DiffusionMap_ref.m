function [U, S] = DiffusionMap_ref(data, Dim, ref)

[n, ~] = size(data) ;

% form reference set
if ~isempty(ref.set)
    ref_set = ref.set;
else   
    if ref.ratio == 0 && sum(ref.idx) == 0
        p = randperm(n);
        ref.idx = sort(p(1:round(n*.01)));  % take 1% data pt to be ref pt
    elseif ref.ratio ~= 0 && sum(ref.idx) == 0
        p = randperm(n);
        ref.idx = p(1:round(n*ref.ratio));
    end   
    ref_set = data(ref.idx, :);
end


% form affinity matrix wrt the reference set
distance = pdist2(data, ref_set).^2;

sig = quantile(max(distance, [], 2), .8) ;
W_ref = exp( - 20 * distance / sig );
W_ref = sparse(W_ref);

%
% make W row stochastic
% form sparse D = D^{-.5}
D = W_ref * sum(W_ref, 1)';
D = D.^(-.5);
D = sparse(1:n, 1:n, D, n, n);
W_ref = D * W_ref;
%}


% SVD on D * W_ref
[U, S, ~] = svds(W_ref, Dim+1);
U = D * U;
S = S.^2;
U = U*S;


end