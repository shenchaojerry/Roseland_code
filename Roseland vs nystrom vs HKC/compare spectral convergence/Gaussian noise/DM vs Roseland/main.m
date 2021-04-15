%% noisy S1 noisy landmark
%% data size=3000, landmark size=300

nystrom_err_inf = [];
ref_err_inf = [];
nystrom_err_l2 = [];
ref_err_l2 = [];
eigenval_nys = [];
eigenval_ref = [];
num_eigvec = 10;
% normalize kernel
fun = @(x) exp(- x.^2);    %kernel
q1 = integral(fun, 0, inf);    %mu_{1,0}^{0}
fun = @(x) (x.^2) .* exp(- x.^2) ./ q1; %mu_{1,2}^{0}
q2 = integral(fun, 0, inf);

iter = 10;

NN = 3000;

%average over many times
for K = 1:iter
    ref.size = 300;
    N = NN + ref.size;
    ref.idx = 0;
    %% generate uniform S^1 data set
    theta = rand(N,1)*2*pi; 
    theta = sort(theta);
    data = [cos(theta) sin(theta)];
    %% noisy data and noisy subset
    data = [data zeros(size(data,1), 100-size(data,2))];
    Noise = randn(size(data)) * .1;
    data = data + Noise ;
    %% get noisy subset
    refind = randperm(N);
    refind = refind(1:ref.size);
    ref.set = data(refind, :);   
    data(refind, :) = [];
    theta(refind) = [];  
    
    %% DM
    dist = pdist2(data, data);
    sig = quantile(dist(:,end).^2, .98) ;
    sig = sig / 10;
    W = exp( - dist.^2 / sig);
    D = sum(W, 2);
    D = D.^(-0.5);
    D = diag(D);
    W = D * W * D;
    W = (W + W') / 2;
    [u_dm, s] = eigs(W, num_eigvec+1);
    s = diag(s); 
    [s, ii] = sort(s, 'descend');
    s = s(2:num_eigvec+1);
    u_dm = u_dm(:, ii);
    u_dm = D * u_dm;
    S_nys = (1 - s) ./ sig * 2 ./ q2;
    eigenval_nys(:, K) = S_nys;
    
    %% roseland
    dist = pdist2(data, ref.set).^2;
    sig = quantile(max(dist, [], 2), .8) ;
    sig = sig / 20;
    [U2, s] = roseland(data,num_eigvec, ref, 1, sig);
    s = s(2:end);
    S_ref = (1 - s) ./ sig ./ q2;  %recover eigenvalues of the laplacian beltromi operator
    eigenval_ref(:, K) = S_ref;
    
    %% get eigenvectors and normalize
    U_nystrom = normc(u_dm(:,2:end));
    U_ref = normc(U2(:,2:end));
    %% alignment with ground truth using consecutive pairs
    n = size(data,1);

    % generate ground truth eigenvecors w.r.t. alignment
    nystrom_aligned = [];
    ref_aligned = [];
    for c = 1:num_eigvec
        
        if mod(c,2) == 1    %this one and next one form a pair
            true1 = sin(linspace(0,(c+1)*pi, n))';
            true2 = cos(linspace(0,(c+1)*pi, n))';
            Y1 = []; Y2 = [];
            for i = 1:length(true1)
                Y1(:, i) = circshift(true1, i);
                Y2(:, i) = circshift(true2, i);
            end
            Y1 = normc(Y1); Y2 = normc(Y2);
            True(:,c) = Y1(:,1);
            True(:,c+1) = Y2(:,1);
            
            %need to check signs (4 possibilities)     
            [m1, a1] = min(max(abs(Y1 - U_nystrom(:, c))) + max(abs(Y2 - U_nystrom(:, c+1))));
            [m2, a2] = min(max(abs(Y1 - U_nystrom(:, c))) + max(abs(Y2 + U_nystrom(:, c+1))));
            [m3, a3] = min(max(abs(Y1 + U_nystrom(:, c))) + max(abs(Y2 - U_nystrom(:, c+1))));
            [m4, a4] = min(max(abs(Y1 + U_nystrom(:, c))) + max(abs(Y2 + U_nystrom(:, c+1))));
            M = [m1 m2 m3 m4]; A = [a1 a2 a3 a4];
            [~, ind] = min(M);
            align_nystrom = A(ind);
            
                
            [m1, a1] = min(max(abs(Y1 - U_ref(:, c))) + max(abs(Y2 - U_ref(:, c+1))));
            [m2, a2] = min(max(abs(Y1 - U_ref(:, c))) + max(abs(Y2 + U_ref(:, c+1))));
            [m3, a3] = min(max(abs(Y1 + U_ref(:, c))) + max(abs(Y2 - U_ref(:, c+1))));
            [m4, a4] = min(max(abs(Y1 + U_ref(:, c))) + max(abs(Y2 + U_ref(:, c+1))));
            M = [m1 m2 m3 m4]; A = [a1 a2 a3 a4];
            [~, ind] = min(M);
            align_ref = A(ind);
           
            
            vec1 = U_nystrom(:,c); vec2 = U_nystrom(:,c+1); 
            nystrom_aligned(:, c) = circshift(vec1, -align_nystrom);
            nystrom_aligned(:, c+1) = circshift(vec2, -align_nystrom);
            
            vec1 = U_ref(:,c); vec2 = U_ref(:,c+1); 
            ref_aligned(:, c) = circshift(vec1, -align_ref);
            ref_aligned(:, c+1) = circshift(vec2, -align_ref);

        end   
    end   


    %% L2 & inf err

    for c = 1:num_eigvec
        inf_true = max(True(:,c));
        if max(abs(nystrom_aligned(:,c)-True(:,c))) < max(abs(nystrom_aligned(:,c)+True(:,c)))
            nystrom_err_inf(K, c) = max(abs(nystrom_aligned(:,c)-True(:,c))) / inf_true;
        else %flip sign
            nystrom_err_inf(K, c) = max(abs(nystrom_aligned(:,c)+True(:,c))) / inf_true;
            nystrom_aligned(:,c) = -nystrom_aligned(:,c);
        end  
        
        if max(abs(ref_aligned(:,c)-True(:,c))) < max(abs(ref_aligned(:,c)+True(:,c)))
            ref_err_inf(K, c) = max(abs(ref_aligned(:,c)-True(:,c))) / inf_true;
        else
            ref_err_inf(K, c) = max(abs(ref_aligned(:,c)+True(:,c))) / inf_true;
            ref_aligned(:, c) = -ref_aligned(:, c);
        end
        
        if sqrt(sum((nystrom_aligned(:,c)-True(:,c)).^2)) < sqrt(sum((nystrom_aligned(:,c)+True(:,c)).^2))
            nystrom_err_l2(K, c) = sqrt(sum((nystrom_aligned(:,c)-True(:,c)).^2));
        else
            nystrom_err_l2(K, c) = sqrt(sum((nystrom_aligned(:,c)+True(:,c)).^2));
        end
        
        if sqrt(sum((ref_aligned(:,c)-True(:,c)).^2)) < sqrt(sum((ref_aligned(:,c)+True(:,c)).^2))
            ref_err_l2(K, c) = sqrt(sum((ref_aligned(:,c)-True(:,c)).^2));
        else
            ref_err_l2(K, c) = sqrt(sum((ref_aligned(:,c)+True(:,c)).^2));
        end

    end
end

nystrom_err_inf = mean(nystrom_err_inf, 1);
nystrom_err_l2 = mean(nystrom_err_l2, 1);
ref_err_inf = mean(ref_err_inf, 1);
ref_err_l2 = mean(ref_err_l2, 1);
eigenval_nys = mean(eigenval_nys, 2);
eigenval_ref = mean(eigenval_ref, 2);

%% eigenvalues superimpose
num_eigenval = 8;
S = eigenval_nys(1:num_eigenval);
S2 = eigenval_ref(1:num_eigenval);

figure('Renderer', 'painters', 'Position', [10 10 900 800]); hold on;
A = [1:(num_eigenval/2); 1:(num_eigenval/2)]; A = A(:);
scatter(1:num_eigenval, A.^2, 200, 'filled');
scatter(1:num_eigenval, S, 200, 'filled');
scatter(1:num_eigenval, S2, 200, 'filled');
axis tight; grid on; 
xticks(1:1:num_eigenval)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
xlabel('The i th eigenvalue', 'fontsize', 35)
[l, hobj, hout, mout] = legend({'truth', 'DM', 'Roseland'}, 'fontsize', 38);
M = findobj(hobj,'type','patch');
set(M,'MarkerSize',20);

%% eigenvalues relative errors
figure('Renderer', 'painters', 'Position', [10 10 900 600]); hold on;
plot(1:num_eigenval, abs(A.^2 - S) ./ (A.^2), '--b s', 'MarkerSize', 20, 'linewidth', 1.5)
plot(1:num_eigenval, abs(A.^2 - S2) ./ (A.^2), '--r *', 'MarkerSize', 20, 'linewidth', 1.5)
axis tight; grid on;
xticks(1:1:num_eigenval)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
xlabel('The i th eigenvalue', 'fontsize', 35)
legend({'DM', 'Roseland'}, 'fontsize', 35);
%title('Relative error', 'fontsize', 20)

%% plot eigenvec errors
num_eigvec = 8;
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
hold on; ylim([.06 6]); xlim([1 8])
plot(1:num_eigvec, nystrom_err_inf(1:num_eigvec), '--b s', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, ref_err_inf(1:num_eigvec), '--b *', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, nystrom_err_l2(1:num_eigvec), '--r s', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, ref_err_l2(1:num_eigvec), '--r *', 'MarkerSize', 15, 'linewidth', 2)
legend({'DM inf-norm','Roseland inf-norm', 'DM 2-norm','Roseland 2-norm'}, 'fontsize', 30)
grid on; 
xticks(1:num_eigvec)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
%title('Relative error', 'fontsize', 20)
xlabel('The i th eigenvector', 'fontsize', 35)


%% compare embeddings
c = linspace(1, NN, NN);
fig_size = 500;
% Data
figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size]);hold on; 
scatter(data(:,1),data(:,2), 20, c, 'filled'); axis tight
scatter3(ref.set(:, 1), ref.set(:, 2), ref.set(:, 3), 90, 'r', 'filled');
axis off;

% DM 
figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size]);
scatter(U_nystrom(:, 1),U_nystrom(:, 2), 30, c, 'filled'); axis tight
axis off;

% ref DM
figure('Renderer', 'painters', 'Position', [10 10 fig_size fig_size]);
scatter(U_ref(:, 1), U_ref(:, 2), 30, c, 'filled'); axis tight
axis off;
