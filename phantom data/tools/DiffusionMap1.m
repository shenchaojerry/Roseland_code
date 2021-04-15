function [U, S] = DiffusionMap1(data, Dim, KNN, alpha)

[n, ~] = size(data) ;


% search for K nearnest neighbors
[index,distance]= knnsearch(data, data, 'k', KNN);
%[index,distance]= knnsearch(data, data, 'k', KNN, 'Distance', 'cosine');

% make sure d(x,x)=0
distance(:,1)=0;     


% find sigma
sig = quantile(distance(:,end).^2, .98) ;
% affinity matrix
ker = exp(-2*distance.^2/sig); 

ii = (1:n)'*ones(1,KNN);
W = sparse(ii, index, ker, n, n);
W = max(W, W');


if alpha
    % alpha-normalization
    D = sum(W, 2).^alpha; 
    W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
end   



% graph laplacian
D = sqrt(sum(W, 2)); 
W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
W = (W + W')./2;

% EVD
[U,S] = eigs(W, Dim + 1);
U = bsxfun(@rdivide, U(:, 2:end), D);
S = S(2:end, 2:end);
U = U * S;

end