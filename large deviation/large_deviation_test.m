M_range = 50:25:800;

%% test large deviation 1

%marginal X and Y are all uniform [0, 6]
fun = @(x,y) x .* y / 36; %density function
expect_true = integral2(fun, 0, 6, 0, 6, 'Method', 'iterated');
power = 2;   %M^power = N

% sample M x N samples from grid
exp_MN_err1 = [];

for i = 1:200
    exp_MN = [];
    for M = M_range
        N = round(M.^power);
        x = 6*rand(N,1);
        y = 6*rand(M,1);
        [X,Y] = meshgrid(x,y);
        A = X .* Y;

        expect_samp = sum(A(:)) / N / M;
        exp_MN = [exp_MN expect_samp];
    end

    exp_MN_err1(i, :) = abs(exp_MN - expect_true) / expect_true;  %relative error
end

% sample M iid samples
exp_M_err1 = [];

for i = 1:200
    exp_M = [];
    for M = M_range

        X = 6*rand(M,1);
        Y = 6*rand(M,1);
        A = X .* Y;

        expect_samp = sum(A) / M;
        exp_M = [exp_M expect_samp];
    end
    exp_M_err1(i,:) = abs(exp_M - expect_true) / expect_true;  %abs error

end

% sample N iid samples
exp_N_err1 = [];
for i = 1:200
    exp_N = [];
    for M = M_range
        N = round(M.^power);
        X = 6*rand(N,1);
        Y = 6*rand(N,1);
        A = X .* Y;

        expect_samp = sum(A) / N;
        exp_N = [exp_N expect_samp];
    end

    exp_N_err1(i,:) = abs(exp_N - expect_true) / expect_true;  %abs error
end



%% test large deviation 2

%marginal X and Y are all uniform [0, 1]
h = 0.01;
fun = @(x,y) exp( - ((y-.5).^2 + (x-y).^2)/ h) / h;  %density function
expect_true = integral2(fun, 0, 1, 0, 1, 'Method', 'iterated');
power = 2;   %M^power = N

% sample MN samples from grid
exp_MN_err2 = [];
for i = 1:200
    exp_MN = [];
    for M = M_range
        N = round(M.^power);

        x = rand(N,1);
        y = rand(M,1);
        [X,Y] = meshgrid(x,y);
        A = exp( - ((Y - .5).^2 + (X - Y).^2) / h) / h;

        expect_samp = sum(A(:)) / N / M;
        exp_MN = [exp_MN expect_samp];
    end

    exp_MN_err2(i,:) = abs(exp_MN - expect_true) / expect_true;  %relative error
end

% sample M iid samples
exp_M_err2 = [];
for i = 1:200
    exp_M = [];
    for M = M_range

        X = rand(M,1);
        Y = rand(M,1);
        A = exp( - ((Y - .5).^2 + (X - Y).^2) / h) / h;

        expect_samp = sum(A) / M;
        exp_M = [exp_M expect_samp];
    end

    exp_M_err2(i,:) = abs(exp_M - expect_true) / expect_true;  %abs error
end

% sample N iid samples
exp_N_err2 = [];
for i = 1:200
    exp_N = [];
    for M = M_range
        N = round(M.^power);

        X = rand(N,1);
        Y = rand(N,1);
        A = exp( - ((Y - .5).^2 + (X - Y).^2) / h) / h;

        expect_samp = sum(A) / N;
        exp_N = [exp_N expect_samp];
    end

    exp_N_err2(i,:) = abs(exp_N - expect_true) / expect_true;  %abs error
end



%% plots
M_range = 50:25:800;

% case 1
figure('Renderer', 'painters', 'Position', [10 10 1000 900]); hold on;
plot(log(M_range), log(mean(exp_MN_err1)), '--o', 'MarkerSize', 15, 'linewidth', 4)
plot(log(M_range), log(mean(exp_M_err1)), '--s', 'MarkerSize', 15, 'linewidth', 4)
plot(log(M_range), log(mean(exp_N_err1)), '--^', 'MarkerSize', 15, 'linewidth', 4)
grid on; grid minor
xt = get(gca, 'XTick');
set(gca, 'FontSize', 35)
axis tight
xlabel('log M', 'fontsize', 35)
ylabel('relative error', 'fontsize', 35)
legend({'M x N grid samples','M iid samples','N iid samples'}, 'fontsize', 28)


% case 2
figure('Renderer', 'painters', 'Position', [10 10 1000 900]); hold on;
plot(log(M_range), log(mean(exp_MN_err2)), '--o', 'MarkerSize', 15, 'linewidth', 4)
plot(log(M_range), log(mean(exp_M_err2)), '--s', 'MarkerSize', 15, 'linewidth', 4)
plot(log(M_range), log(mean(exp_N_err2)), '--^', 'MarkerSize', 15, 'linewidth', 4)
grid on; grid minor
xt = get(gca, 'XTick');
set(gca, 'FontSize', 35)
axis tight
xlabel('log M', 'fontsize', 35)
ylabel('relative error', 'fontsize', 35)
legend({'M x N grid samples','M iid samples','N iid samples'}, 'fontsize', 28)


