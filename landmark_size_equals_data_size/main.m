%% compare eigenvalue & eigenvec: DM vs Roseland

dm_err_inf = [];
ref_err_inf = [];
dm_err_l2 = [];
ref_err_l2 = [];
eigenval_dm = [];
eigenval_ref = [];
num_eigvec = 100;
beta = 1;

% normalize kernel
fun = @(x) exp(- x.^2);    %kernel
q1 = integral(fun, 0, inf);    %mu_{1,0}^{0}
fun = @(x) (x.^2) .* exp(- x.^2) ./ q1; %mu_{1,2}^{0}
q2 = integral(fun, 0, inf);

iter = 1;
%average over K many times
for K = 1:iter
    ref.size = round(2500^beta);
    N = 2500 + ref.size;
    ref.idx = 0;
    %% generate uniform S^1 data set
    theta = rand(N,1)*2*pi;
    theta = sort(theta);
    data = [cos(theta) sin(theta)];

    %% noisy data and noisy subset
    %{
    data = [data zeros(size(data,1), 100-size(data,2))];
    Noise = randn(size(data)) * .1;
    data = data + Noise ;
    %}

    %% get subset
    refind = randperm(N);
    refind = refind(1:ref.size);
    ref.set = data(refind, :);     %for Nystrom and ref
    data(refind, :) = [];

    %% noisy data, clean subset
    %{
    data = [data zeros(size(data,1), 100-size(data,2))];
    ref.set = [ref.set zeros(size(ref.set,1), 100-size(ref.set,2))];
    Noise = randn(size(data)) * .1;
    data = data + Noise ;
    %}

    %% DM
    h_dm = 0.0169;   %m=n=2500
    dist = pdist2(data, data);
    W = exp( - dist.^2 / h_dm);
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
    S_dm = (1 - s) ./ h_dm * 2 ./ q2;
    eigenval_dm(:, K) = S_dm;

    %% Roseland
    h_ref = 0.01; %m=n=2500
    dist = pdist2(data, ref.set);
    W_ref = exp( - dist.^2 / h_ref );

    D = W_ref * sum(W_ref, 1)';
    D = D.^(-.5);
    D = diag(D);
    W_ref = D * W_ref;
    [u_ref, s, ~] = svds(W_ref, num_eigvec+1);
    s = diag(s);
    s = s(2:num_eigvec+1);
    s = s.^2;
    u_ref = D * u_ref;
    S_ref = (1 - s) ./ h_ref ./ q2;
    eigenval_ref(:, K) = S_ref;



    %% get eigenvectors and normalize
    U_dm = u_dm(:, 2:end);
    U_ref = u_ref(:, 2:end);
    U_dm = normc(U_dm);
    U_ref = normc(U_ref);

    %% alignment with ground truth using consecutive pairs
    n = size(data,1);

    % generate ground truth eigenvecors w.r.t. alignment
    True_dm_aligned = [];
    True_ref_aligned = [];
    for c = 1:num_eigvec

        if mod(c,2) == 1    %this one and next one form a pair
            true1 = sin(linspace(0,(c+1)*pi, n))';
            true2 = cos(linspace(0,(c+1)*pi, n))';
            Y1 = []; Y2 = [];
            for i = 1:n
                Y1(:, i) = circshift(true1, i);
                Y2(:, i) = circshift(true2, i);
            end
            Y1 = normc(Y1); Y2 = normc(Y2);
            %need to check signs (4 possibilities)
            [m1, a1] = min(max(abs(Y1 - U_dm(:, c))) + max(abs(Y2 - U_dm(:, c+1))));
            [m2, a2] = min(max(abs(Y1 - U_dm(:, c))) + max(abs(Y2 + U_dm(:, c+1))));
            [m3, a3] = min(max(abs(Y1 + U_dm(:, c))) + max(abs(Y2 - U_dm(:, c+1))));
            [m4, a4] = min(max(abs(Y1 + U_dm(:, c))) + max(abs(Y2 + U_dm(:, c+1))));
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

            True_dm_aligned(:, c) = Y1(:, align_nystrom);
            True_dm_aligned(:, c+1) = Y2(:, align_nystrom);
            True_ref_aligned(:, c) = Y1(:, align_ref);
            True_ref_aligned(:, c+1) = Y2(:, align_ref);
        end
    end


    %% L2 & inf err

    for c = 1:num_eigvec
        inf_true = max(True_ref_aligned(:,c));
        dm_err_inf(K, c) = min(max(abs(True_dm_aligned(:,c)-U_dm(:,c))), max(abs(True_dm_aligned(:,c)+U_dm(:,c)))) / inf_true;
        ref_err_inf(K, c) = min(max(abs(True_ref_aligned(:,c)-U_ref(:,c))), max(abs(True_ref_aligned(:,c)+U_ref(:,c)))) / inf_true;

        dm_err_l2(K, c) = min(sqrt(sum((True_dm_aligned(:,c)-U_dm(:,c)).^2)), sqrt(sum((True_dm_aligned(:,c)+U_dm(:,c)).^2)));
        ref_err_l2(K, c) = min(sqrt(sum((True_ref_aligned(:,c)-U_ref(:,c)).^2)), sqrt(sum((True_ref_aligned(:,c)+U_ref(:,c)).^2)));
    end
end

dm_err_inf = mean(dm_err_inf, 1);
dm_err_l2 = mean(dm_err_l2, 1);
ref_err_inf = mean(ref_err_inf, 1);
ref_err_l2 = mean(ref_err_l2, 1);


%%
%{
cd('/home/grad/chaos/Dropbox/ref DM project/refDM vs DM vs Nystrom/compare DM vs refDM')
save('dm_err_inf.mat', 'dm_err_inf');
save('dm_err_l2.mat', 'dm_err_l2');
save('ref_err_inf.mat', 'ref_err_inf');
save('ref_err_l2.mat', 'ref_err_l2');
save('U_dm.mat', 'U_dm')
save('U_ref.mat', 'U_ref')
save('True_dm_aligned.mat', 'True_dm_aligned')
save('True_ref_aligned.mat', 'True_ref_aligned')
save('eigenval_dm.mat', 'eigenval_dm');
save('eigenval_ref.mat', 'eigenval_ref');
%}
