%% geodesic recovery by full DM, nys, Roseland

range = 100:2200;
subset_size = 50;
%subset_size = 300;
top = 40;


% parameters
t = 0.015; %diffusion time
h_dm = 0.1; 
h_nys = ( sqrt(log(subset_size)) / 0.3146 / sqrt(subset_size) )^(4/3);  %guided by the variance analysis
h_ref = h_nys / 1.5;


%1st column is h & t
%2nd is err for geodesic between x_{i} and x_{i+1}
%3rd is err for geodesic between x_{i} and x_{i+10}
dm_err_knn = cell(1, 3);
nys_err_knn = cell(1, 3);
ref_err_knn = cell(1, 3);

% normalize kernel
fun = @(x) exp(- x.^2);    %kernel
q1 = integral(fun, 0, inf);    %mu_{1,0}^{0}
fun = @(x) (x.^2) .* exp(- x.^2) ./ q1; %mu_{1,2}^{0}
q2 = integral(fun, 0, inf);

iter = 5;  

count = 0;

for i = 1:iter
    count = count + 1;

    % uniform S^1
    N = 2500;
    theta = rand(N,1)*2*pi;
    theta = sort(theta);
    data = [cos(theta) sin(theta)];
    
    % subset
    ref.size = subset_size;
    theta1 = rand(ref.size,1)*2*pi;
    theta1 = sort(theta1);
    ref.set = [cos(theta1) sin(theta1)];
    
    if count == 1
        dm_err_knn{count,1} = sprintf('h=%.3f, t=%.3f', h_dm,t);
        nys_err_knn{count,1} = sprintf('h=%.3f, t=%.3f', h_nys,t);
        ref_err_knn{count,1} = sprintf('h=%.3f, t=%.3f', h_ref,t);
    end
    
    %% DM
    [u, s] = DiffusionMap(data, 0, top, N, 0, h_dm);
    u = u * diag(1 ./ max(u)); 
    
    %A = [1:top/2; 1:top/2]; A = A(:); S_dm = A.^2;   %true eigenval

    %
    %get back eigenvalues of S1 laplacian
    S_dm = (1 - s) ./ h_dm * 2 ./ q2;
    S_dm = S_dm(1:top);
    %}
    
    % heat kernel coordinates
    CONST = ((2*t).^(0.75)) * sqrt(2) * ((4*pi).^(0.25));
    u_dm = (1/sqrt(pi)) * CONST * u * diag(exp( - t * S_dm));
    
    %% Nystrom
    
    [u_ext, s] = Nystrom(data, ref.set, top, 0, h_nys, 1);
    s = s(2:end);
    V = u_ext(ref.size+1:end,:);
    V = V(:,2:end);
    V = V * diag(1 ./ max(abs(V))); 
    
    
    %A = [1:top/2; 1:top/2]; A = A(:); S_nystrom = A.^2;   %true eigenval
    
    %
    % get back eigenvalues of S^1 laplacian
    S_nystrom = (1 - s) ./ h_nys * 2 ./ q2;
    S_nystrom = S_nystrom(1:top);
    %}
    
    % heat kernel coordinates
    CONST = ((2*t).^(0.75)) * sqrt(2) * ((4*pi).^(0.25));
    u_nys = (1/sqrt(pi)) * CONST * V * diag(exp( - t * S_nystrom));
    
    
    %% roseland
    
    [u, s] = roseland(data, top, ref, 0, h_ref);
    u = u(:, 2:end);
    s = s(2:end);
    u = u * diag(1 ./ max(abs(u))); 
    
    S_ref = (1 - s) ./ h_ref ./ q2;
    S_ref = S_ref(1:top);
    
    
    %A = [1:top/2; 1:top/2]; A = A(:); S_ref = A.^2;   %true eigenval
    
    % heat kernel coordinates
    CONST = ((2*t).^(0.75)) * sqrt(2) * ((4*pi).^(0.25));
    u_ref = (1/sqrt(pi)) * CONST * u * diag(exp( - t * S_ref));

    
    %% geodesic recovery
    ind = range; 
    true1 = theta(ind+1) - theta(ind); %true geodesic distances between x_{i} and x_{i+1}
    true2 = theta(ind+10) - theta(ind); %true geodesic distances between x_{i} and x_{i+10}
    
    for ii = 1:top  
        
        % dm
        coord1 = u_dm(ind, 1:ii);
        coord2 = u_dm(ind + 1, 1:ii);
        coord3 = u_dm(ind + 10, 1:ii);
        %knn=1 err
        dm_err_knn{count, 2}(ii) = mean(abs(sqrt(sum((coord2-coord1).^2, 2)) - true1) ./ true1);
        %knn=10 err
        dm_err_knn{count, 3}(ii) = mean(abs(sqrt(sum((coord3-coord1).^2, 2)) - true2) ./ true2);
        
        
        % nystrom
        coord1 = u_nys(ind, 1:ii);
        coord2 = u_nys(ind + 1, 1:ii);
        coord3 = u_nys(ind + 10, 1:ii);
        %knn=1 err
        nys_err_knn{count, 2}(ii) = mean(abs(sqrt(sum((coord2-coord1).^2, 2)) - true1) ./ true1);
        %knn=10 err
        nys_err_knn{count, 3}(ii) = mean(abs(sqrt(sum((coord3-coord1).^2, 2)) - true2) ./ true2);
        
        
        % ref dm
        coord1 = u_ref(ind, 1:ii);
        coord2 = u_ref(ind + 1, 1:ii);
        coord3 = u_ref(ind + 10, 1:ii);
        %knn=1 err
        ref_err_knn{count, 2}(ii) = mean(abs(sqrt(sum((coord2-coord1).^2, 2)) - true1) ./ true1);
        %knn=10 err
        ref_err_knn{count, 3}(ii) = mean(abs(sqrt(sum((coord3-coord1).^2, 2)) - true2) ./ true2);       
        
    end
end


%% plots
%% get errors
top = size(dm_err_knn{1,2},2);
% dm
dm1 = zeros(1,top);
for i = 1:size(dm_err_knn,1)
    dm1 = dm1 + dm_err_knn{i,2};
end  
dm1 = dm1 ./ size(dm_err_knn,1);
% nys 
Nys1 = zeros(1,top);
for i = 1:size(nys_err_knn,1)
    Nys1 = Nys1 + nys_err_knn{i,2};
end  
Nys1 = Nys1 ./ size(nys_err_knn,1);
%  Roseland
Ref1 = zeros(1,top);
for i = 1:size(ref_err_knn,1)
    Ref1 = Ref1 + ref_err_knn{i,2};
end  
Ref1 = Ref1 ./ size(ref_err_knn,1);

%% compare eigevalues
num_eigenval = 26;
figure('Renderer', 'painters', 'Position', [10 10 900 800]); hold on;
%subplot(1,2,1)
A = [1:(num_eigenval/2); 1:(num_eigenval/2)]; A = A(:);
hold on; axis tight
scatter(1:num_eigenval, A.^2, 100, 'filled');
scatter(1:num_eigenval, S_dm(1:num_eigenval), 100, 'filled');
scatter(1:num_eigenval, S_nystrom(1:num_eigenval), 100, 'filled');
scatter(1:num_eigenval, S_ref(1:num_eigenval), 100, 'filled');
grid on; 
xticks(1:3:num_eigenval)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
xlabel('The i th eigenvalue', 'fontsize', 40)
[l, hobj, hout, mout] = legend({'truth', 'DM','Nystrom', 'Roseland'}, 'fontsize', 35);
M = findobj(hobj,'type','patch');
set(M,'MarkerSize',20);

%% compare 1-NN errors 
%
dm_vec = 1:18;
nys_vec = 1:num_eigenval;
ref_vec = 1:15;
%}

%{
till = num_eigenval;
dm_vec = 1:3:till;
nys_vec = 1:3:till;
ref_vec = 1:3:till;
%}

figure('Renderer', 'painters', 'Position', [10 10 1200 800]); hold on;
%subplot(1,2,2); 
hold on; 
plot(dm_vec, dm1(dm_vec), '--^', 'MarkerSize', 25, 'linewidth', 3)
plot(nys_vec, Nys1(nys_vec), '-- s', 'MarkerSize', 25, 'linewidth', 3)
plot(ref_vec, Ref1(ref_vec), '-- *', 'MarkerSize', 25, 'linewidth', 3)
axis tight
legend({'DM 1-NN', 'Nystrom 1-NN','Roseland 1-NN'}, 'fontsize', 35)
grid on; 
xticks(1:2:num_eigenval)
%xticks(1:3:till)
yticks(0:.08:.95)
set(gca, 'FontSize', 30)
%ylabel('Relative errors', 'fontsize', 20)
xlabel('Number of eigenvectors used for embedding', 'fontsize', 35)
