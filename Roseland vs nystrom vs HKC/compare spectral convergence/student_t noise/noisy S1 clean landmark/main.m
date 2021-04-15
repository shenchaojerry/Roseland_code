%% noisy S1 clean landmark
%% data size=90000, landmark size=300

embed_dim = 1000;
nystrom_err_inf = [];
ref_err_inf = [];
nystrom_err_l2 = [];
ref_err_l2 = [];
num_eigvec = 26;
subset_size = 250;
roseland_subset_size = 200;

iter = 1;

%average over many times
for K = 1:iter
    N = 90000 + subset_size;
    %% generate uniform S^1 data set
    theta = rand(N,1)*2*pi; 
    theta = sort(theta);
    data = [cos(theta) sin(theta)];
    
    %% get subset
    refind = randperm(N);
    refind = refind(1:subset_size);
    subset = data(refind, :);   
    data(refind, :) = [];
    theta(refind) = []; 

    %% noisy data, clean subset
    data = [data zeros(size(data,1), embed_dim-size(data,2))];
    subset = [subset zeros(size(subset,1), embed_dim-size(subset,2))];  
    Noise = trnd(4, size(data)) / sqrt(2) * .1; % student t
    data = data + Noise ;
    %% Nystrom
    sample = subset;
    sample_size = subset_size; 
    [V, ~] = Nystrom(data, sample, num_eigvec, 1);
    V = V(sample_size+1:end, :);
    
    %% roseland
    ref_set = sample(1:roseland_subset_size, :);
    [U2, ~] = roseland(data, num_eigvec, ref_set, 1);

    %% HKC
    u_Haddad = HKC(data, num_eigvec, ref_set, 1);
    
    %% get eigenvectors and normalize
    % we down sample to speed up phase alignment & err calculation
    pick = 1:15:90000;
    U_nystrom = V(pick, 2:end);
    U_nystrom = normc(U_nystrom);
    U_ref = U2(pick, 2:end);
    U_ref = normc(U_ref);
    U_Haddad = u_Haddad(pick, 2:end);
    U_Haddad = normc(U_Haddad);
    %% alignment with ground truth using consecutive pairs
    n = size(data,1);
    % generate ground truth eigenvecors w.r.t. alignment
    nystrom_aligned = [];
    ref_aligned = [];
    Haddad_aligned = [];
    True = [];
    
    for c = 1:num_eigvec
        
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
            
            [m1, a1] = min(max(abs(Y1 - U_Haddad(:, c))) + max(abs(Y2 - U_Haddad(:, c+1))));
            [m2, a2] = min(max(abs(Y1 - U_Haddad(:, c))) + max(abs(Y2 + U_Haddad(:, c+1))));
            [m3, a3] = min(max(abs(Y1 + U_Haddad(:, c))) + max(abs(Y2 - U_Haddad(:, c+1))));
            [m4, a4] = min(max(abs(Y1 + U_Haddad(:, c))) + max(abs(Y2 + U_Haddad(:, c+1))));
            M = [m1 m2 m3 m4]; A = [a1 a2 a3 a4];
            [~, ind] = min(M);
            align_haddad = A(ind);
            
            vec1 = U_nystrom(:,c); vec2 = U_nystrom(:,c+1); 
            nystrom_aligned(:, c) = circshift(vec1, -align_nystrom);
            nystrom_aligned(:, c+1) = circshift(vec2, -align_nystrom);
            
            vec1 = U_ref(:,c); vec2 = U_ref(:,c+1); 
            ref_aligned(:, c) = circshift(vec1, -align_ref);
            ref_aligned(:, c+1) = circshift(vec2, -align_ref);
            
            vec1 = U_Haddad(:,c); vec2 = U_Haddad(:,c+1); 
            Haddad_aligned(:, c) = circshift(vec1, -align_haddad);
            Haddad_aligned(:, c+1) = circshift(vec2, -align_haddad);
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
        
        if max(abs(Haddad_aligned(:,c)-True(:,c))) > max(abs(Haddad_aligned(:,c)+True(:,c)))
            Haddad_aligned(:,c) = -Haddad_aligned(:,c);
        end  
    end
end

nystrom_err_inf = mean(nystrom_err_inf, 1);
nystrom_err_l2 = mean(nystrom_err_l2, 1);
ref_err_inf = mean(ref_err_inf, 1);
ref_err_l2 = mean(ref_err_l2, 1);

%% plot eigenvec errors
num_eigvec = 12;
figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
hold on; ylim([.06 2]); xlim([1 12])
plot(1:num_eigvec, nystrom_err_inf(1:num_eigvec), '--b s', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, ref_err_inf(1:num_eigvec), '--b *', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, nystrom_err_l2(1:num_eigvec), '--r s', 'MarkerSize', 15, 'linewidth', 2)
plot(1:num_eigvec, ref_err_l2(1:num_eigvec), '--r *', 'MarkerSize', 15, 'linewidth', 2)
legend({'Nystrom inf-norm','Roseland inf-norm', 'Nystrom 2-norm','Roseland 2-norm'}, 'fontsize', 30)
grid on; 
xticks(1:num_eigvec)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
%title('Relative error', 'fontsize', 20)
xlabel('The i th eigenvector', 'fontsize', 35)


%% compare first 6 eigenvectors

line1 = 1; line2 = 3;
figure('Renderer', 'painters', 'Position', [10 10 1700 660]);
subplot(3,6,1); hold on; axis tight
plot(nystrom_aligned(:,1), 'linewidth', line1)
plot(True(:,1), 'linewidth', line2)
ylabel('Nystrom', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
%legend({'nystrom', 'truth'}, 'fontsize', 10)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,2); hold on; axis tight
plot(nystrom_aligned(:,2), 'linewidth', line1)
plot(True(:,2), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,3); hold on; axis tight
plot(nystrom_aligned(:,3), 'linewidth', line1)
plot(True(:,3), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,4); hold on; axis tight
plot(nystrom_aligned(:,4), 'linewidth', line1)
plot(True(:,4), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,5); hold on; axis tight
plot(nystrom_aligned(:,5), 'linewidth', line1)
plot(True(:,5), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,6); hold on; axis tight
plot(nystrom_aligned(:,6), 'linewidth', line1)
plot(True(:,6), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


subplot(3,6,7); hold on; axis tight
plot(Haddad_aligned(:,1), 'linewidth', line1)
plot(True(:,1), 'linewidth', line2)
ylabel('HKC', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,8); hold on; axis tight
plot(Haddad_aligned(:,2), 'linewidth', line1)
plot(True(:,2), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,9); hold on; axis tight
plot(Haddad_aligned(:,3), 'linewidth', line1)
plot(True(:,3), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,10); hold on; axis tight
plot(Haddad_aligned(:,4), 'linewidth', line1)
plot(True(:,4), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,11); hold on; axis tight
plot(Haddad_aligned(:,5), 'linewidth', line1)
plot(True(:,5), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,12); hold on; axis tight
plot(Haddad_aligned(:,6), 'linewidth', line1)
plot(True(:,6), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


subplot(3,6,13); hold on; axis tight
plot(ref_aligned(:,1), 'linewidth', line1)
plot(True(:,1), 'linewidth', line2)
ylabel('Roseland', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
%legend({'ref', 'truth'}, 'fontsize', 10)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,14); hold on; axis tight
plot(ref_aligned(:,2), 'linewidth', line1)
plot(True(:,2), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,15); hold on; axis tight
plot(ref_aligned(:,3), 'linewidth', line1)
plot(True(:,3), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,16); hold on; axis tight
plot(ref_aligned(:,4), 'linewidth', line1)
plot(True(:,4), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,17); hold on; axis tight
plot(ref_aligned(:,5), 'linewidth', line1)
plot(True(:,5), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,18); hold on; axis tight
plot(ref_aligned(:,6), 'linewidth', line1)
plot(True(:,6), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


%% compare 7-12 eigenvectors

line1 = 1; line2 = 3;
figure('Renderer', 'painters', 'Position', [10 10 1700 660]);
subplot(3,6,1); hold on; axis tight
plot(nystrom_aligned(:,7), 'linewidth', line1)
plot(True(:,7), 'linewidth', line2)
ylabel('Nystrom', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
%legend({'nystrom', 'truth'}, 'fontsize', 10)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,2); hold on; axis tight
plot(nystrom_aligned(:,8), 'linewidth', line1)
plot(True(:,8), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,3); hold on; axis tight
plot(nystrom_aligned(:,9), 'linewidth', line1)
plot(True(:,9), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,4); hold on; axis tight
plot(nystrom_aligned(:,10), 'linewidth', line1)
plot(True(:,10), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,5); hold on; axis tight
plot(nystrom_aligned(:,11), 'linewidth', line1)
plot(True(:,11), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,6); hold on; axis tight
plot(nystrom_aligned(:,12), 'linewidth', line1)
plot(True(:,12), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


subplot(3,6,7); hold on; axis tight
plot(Haddad_aligned(:,7), 'linewidth', line1)
plot(True(:,7), 'linewidth', line2)
ylabel('HKC', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,8); hold on; axis tight
plot(Haddad_aligned(:,8), 'linewidth', line1)
plot(True(:,8), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,9); hold on; axis tight
plot(Haddad_aligned(:,9), 'linewidth', line1)
plot(True(:,9), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,10); hold on; axis tight
plot(Haddad_aligned(:,10), 'linewidth', line1)
plot(True(:,10), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,11); hold on; axis tight
plot(Haddad_aligned(:,11), 'linewidth', line1)
plot(True(:,11), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,12); hold on; axis tight
plot(Haddad_aligned(:,12), 'linewidth', line1)
plot(True(:,12), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


subplot(3,6,13); hold on; axis tight
plot(ref_aligned(:,7), 'linewidth', line1)
plot(True(:,7), 'linewidth', line2)
ylabel('Roseland', 'fontsize', 25)
%title('1 st eigenvector', 'fontsize', 15)
%legend({'ref', 'truth'}, 'fontsize', 10)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,14); hold on; axis tight
plot(ref_aligned(:,8), 'linewidth', line1)
plot(True(:,8), 'linewidth', line2)
%title('2 nd eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,15); hold on; axis tight
plot(ref_aligned(:,9), 'linewidth', line1)
plot(True(:,9), 'linewidth', line2)
%title('3 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,16); hold on; axis tight
plot(ref_aligned(:,10), 'linewidth', line1)
plot(True(:,10), 'linewidth', line2)
%title('4 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,17); hold on; axis tight
plot(ref_aligned(:,11), 'linewidth', line1)
plot(True(:,11), 'linewidth', line2)
%title('5 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])

subplot(3,6,18); hold on; axis tight
plot(ref_aligned(:,12), 'linewidth', line1)
plot(True(:,12), 'linewidth', line2)
%title('6 th eigenvector', 'fontsize', 15)
set(gca,'xtick',[]); set(gca,'ytick',[])


%% compare embeddings
c = linspace(1, 9000, 9000);
% Data
pick = 1:10:90000;
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);hold on; 
scatter(data(pick, 1),data(pick, 2), 30, c, 'filled'); axis tight
scatter(subset(:, 1), subset(:, 2), 90, 'r', 'filled');
axis off;

% Nystrom 
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
scatter(U_nystrom(:, 1),U_nystrom(:, 2), 30, c, 'filled'); axis tight
axis off;
% HKC
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
scatter(U_Haddad(:, 1), U_Haddad(:, 2), 30, c, 'filled'); axis tight
axis off;
% ref DM
figure('Renderer', 'painters', 'Position', [10 10 1000 1000]);
scatter(U_ref(:, 1), U_ref(:, 2), 30, c, 'filled'); axis tight
axis off;


