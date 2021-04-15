function [u, s] = randKNN(data, randknn, dim, denoise, sig)

% randknn: # of pts randomly selected from KNN

n = size(data,1);
KNN = 700;
[index, distance]= knnsearch(data, data, 'k', KNN+1);

if denoise
    distance(:,1) = 0;  %remove diagonal
end  

if nargin < 5
    sig = quantile(distance(:,end).^2, .5);
end

ind = [];
dist = [];
for i = 1:n
    temp = randperm(KNN,randknn)+1;
    ind(i,:) = [1 index(i,temp)];
    dist(i,:) = [1 distance(i,temp)];
end  

% affinity matrix
ker = exp(- dist.^2 / sig); 
ii = (1:n)'*ones(1,randknn+1);
W = sparse(ii, ind, ker, n, n);
W = max(W, W');
 

D = sum(W, 2);
D = D.^(-0.5);
D = diag(D);
W = D * W * D;
W = (W + W')./2;
[u, s] = eigs(W, dim + 1);
u = D * u;
s = diag(s);

end