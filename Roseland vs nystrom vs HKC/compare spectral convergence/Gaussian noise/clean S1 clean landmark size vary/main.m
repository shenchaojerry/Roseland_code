%% clean data size = 90000
%% clean landmark size vary from 100 to 1000

nystrom_err_inf = [];
ref_err_inf = [];
nystrom_err_l2 = [];
ref_err_l2 = [];


count = 0;
num_eigvec = 12;
ref_size = 100:100:1000;

for K = ref_size
    nystrom_err_inf_temp = [];
    ref_err_inf_temp = [];
    nystrom_err_l2_temp = [];
    ref_err_l2_temp = [];

    count = count + 1;
    
    %average over many times
    for ii = 1:20
        ref.size = K;
        N = 90000 + ref.size;
        ref.idx = 0;
        %% generate uniform S^1 data set and subset
        theta = rand(N,1)*2*pi; 
        theta = sort(theta);
        data = [cos(theta) sin(theta)];

        refind = randperm(N);
        refind = refind(1:ref.size);
        ref.set = data(refind, :);  
        data(refind, :) = [];
        theta(refind) = [];
        
        %% Nystrom
        sample = ref.set;
        sample_size = ref.size; 
        [V, ~] = Nystrom3(data, sample, num_eigvec);
        V = V(sample_size+1:end, :);
        %% roseland
        [U2, ~] = roseland(data,num_eigvec, ref, 0);
     
        %% get eigenvectors and normalize
        pick = 1:15:90000;
        U_nystrom = V(pick, 2:end);
        U_nystrom = normc(U_nystrom);
        U_ref = U2(pick, 2:end);
        U_ref = normc(U_ref);
        %% alignment with ground truth using consecutive pairs
        n = size(data,1);
        
        % generate ground truth eigenvecors w.r.t. alignment
        True_nystrom_aligned = [];
        True_ref_aligned = [];
        for c = 1:num_eigvec  %this one and next one form a pair
            if mod(c,2) == 1    %this one and next one form a pair
                true1 = sin(linspace(0,(c+1)*pi, n))';
                true2 = cos(linspace(0,(c+1)*pi, n))';
                true1 = true1(pick);
                true2 = true2(pick);
                Y1 = []; Y2 = [];
                for i = 1:length(true1)
                    Y1(:, i) = circshift(true1, i);
                    Y2(:, i) = circshift(true2, i);
                end
                Y1 = normc(Y1); Y2 = normc(Y2);
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
                
                True_nystrom_aligned(:, c) = Y1(:, align_nystrom);
                True_nystrom_aligned(:, c+1) = Y2(:, align_nystrom);
                True_ref_aligned(:, c) = Y1(:, align_ref);
                True_ref_aligned(:, c+1) = Y2(:, align_ref);
            end   
        end   

        %% L2 & inf err
        for c = 1:num_eigvec
            inf_true = max(True_ref_aligned(:,c));
            nystrom_err_inf_temp(ii, c) = min(max(abs(True_nystrom_aligned(:,c)-U_nystrom(:,c))), max(abs(True_nystrom_aligned(:,c)+U_nystrom(:,c)))) / inf_true;
            ref_err_inf_temp(ii, c) = min(max(abs(True_ref_aligned(:,c)-U_ref(:,c))), max(abs(True_ref_aligned(:,c)+U_ref(:,c)))) / inf_true;

            nystrom_err_l2_temp(ii, c) = min(sqrt(sum((True_nystrom_aligned(:,c)-U_nystrom(:,c)).^2)), sqrt(sum((True_nystrom_aligned(:,c)+U_nystrom(:,c)).^2)));
            ref_err_l2_temp(ii, c) = min(sqrt(sum((True_ref_aligned(:,c)-U_ref(:,c)).^2)), sqrt(sum((True_ref_aligned(:,c)+U_ref(:,c)).^2)));

        end
    end
    
    nystrom_err_inf(count, 1:num_eigvec) = mean(nystrom_err_inf_temp, 1);
    nystrom_err_l2(count, 1:num_eigvec) = mean(nystrom_err_l2_temp, 1);
    ref_err_inf(count, 1:num_eigvec) = mean(ref_err_inf_temp, 1);
    ref_err_l2(count, 1:num_eigvec) = mean(ref_err_l2_temp, 1);
end




%% plots nystrom vs ref
ref_size = 100:100:1000;
A = 1;
B = 10;

%first eigenvector err
figure('Renderer', 'painters', 'Position', [10 10 1100 800]);
hold on; axis tight;
plot(ref_size, nystrom_err_inf(:,A), '--b s', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, ref_err_inf(:,A), '--b *', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, nystrom_err_l2(:,A), '--r s', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, ref_err_l2(:,A), '--r *', 'MarkerSize', 15, 'linewidth', 2)
legend({'Nystrom inf-norm','Roseland inf-norm', 'Nystrom 2-norm','Roseland 2-norm'}, 'fontsize', 28)
grid on;
xticks(100:100:1000)
xticklabels({'100','200','300','400','500','600','700','800','900','1000'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
xlabel('Size of subset', 'fontsize', 35)

%10 th eigenvector err
figure('Renderer', 'painters', 'Position', [10 10 1100 800]);
hold on; axis tight;
plot(ref_size, nystrom_err_inf(:,B), '--b s', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, ref_err_inf(:,B), '--b *', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, nystrom_err_l2(:,B), '--r s', 'MarkerSize', 15, 'linewidth', 2)
plot(ref_size, ref_err_l2(:,B), '--r *', 'MarkerSize', 15, 'linewidth', 2)
legend({'Nystrom inf-norm','Roseland inf-norm', 'Nystrom 2-norm','Roseland 2-norm'}, 'fontsize', 28)
grid on;
xticks(100:100:1000)
xticklabels({'100','200','300','400','500','600','700','800','900','1000'})
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
xlabel('Size of subset', 'fontsize', 35)
