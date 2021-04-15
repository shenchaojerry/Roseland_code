function [U, S] = DiffusionMap(data, denoise, Dim, KNN, alpha, sig)

[n, ~] = size(data) ;


% search for K nearnest neighbors
[index,distance]= knnsearch(data, data, 'k', KNN);

% make sure d(x,x)=0
distance(:,1)=0;     
  
if nargin < 6
    % find sigma
    sig = quantile(distance(:,end).^2, .98) ;
end

% affinity matrix
ker = exp(- distance.^2 / sig); 

ii = (1:n)'*ones(1,KNN);
W = sparse(ii, index, ker, n, n);
W = max(W, W');

if denoise
W(logical(eye(size(W)))) = 0;
end   

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
S = diag(S);
%U = U * S;

end