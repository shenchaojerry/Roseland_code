%% generate non-uniform S^1 data set via Rejection sampling
N = 6000;
%true pdf for angle theta
pdf = @(x) exp( - (x - pi).^2 ./ 5) + .2;
area = integral(pdf, 0, 2*pi); 
pdf = @(x) pdf(x) ./ area;  %normalize to a real pdf

subplot(2,2,1)
figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
x = linspace(0, 2*pi,1000);
plot(x, pdf(x), 'linewidth', 10); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
xlabel('\theta', 'fontsize',25)
title('pdf of data','fontsize',25)


%create many samples on interval [0, 2pi]
x_samples = 2 * pi * rand(2 * N, 1);
%evaluate for each sample
sample_value = pdf(x_samples);
%rejection step
max_value = max(sample_value);
accepted = rand(2 * N, 1) < (sample_value / max_value);
theta = x_samples(accepted);

subplot(2,2,2)
figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
hist(theta); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
xlabel('\theta', 'fontsize',25)
title('simulated theta','fontsize',25)

theta = sort(theta);
data = [cos(theta) sin(theta)];

data_forplot = randperm(size(data,1));
data_forplot = data_forplot(1:1000);


%compute reference pdf
f = @(x) 1 ./ x.^2;
q = @(x) f(pdf(x));
area = integral(q, 0, 2*pi); 
q = @(x) q(x) ./ area;  %normalize to a real pdf

subplot(2,2,3)
figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
plot(x, q(x),'linewidth', 10); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
xlabel('\theta', 'fontsize',25)
title('pdf of reference set','fontsize',25)
%export_fig('pdf_ref','-transparent','-eps')


%create many samples on interval [0, 2pi]
y_samples = 2 * pi * rand(round(2 * N), 1);
%evaluate for each sample
sample_value = q(y_samples);
%rejection step
max_value = max(sample_value);
accepted = rand(round(2 * N), 1) < (sample_value / max_value);
theta_ref = y_samples(accepted);
theta_ref = sort(theta_ref);

ref_set = [cos(theta_ref) sin(theta_ref)];

ref_forplot = randperm(size(ref_set,1));
ref_forplot = ref_forplot(1:200);

subplot(2,2,4)
figure('Renderer', 'painters', 'Position', [10 10 1200 800]);
hist(theta_ref); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
xlabel('\theta', 'fontsize',25)
title('simulated theta','fontsize',25)


%plot data and ref data
subplot(1,2,1)
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(data(data_forplot, 1), data(data_forplot, 2), 10, 'filled'); axis tight;
ax = gca; ax.FontSize = 35; 
%title(sprintf('data set, size = %d', size(data,1)), 'fontsize', 25)

subplot(1,2,2)
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(ref_set(ref_forplot, 1), ref_set(ref_forplot, 2), 10, 'filled');
ax = gca; ax.FontSize = 30; 
%title(sprintf('reference set, size = %d', size(ref_set,1)), 'fontsize', 25)


%% randomly choose subsets as ref set

ref.ratio = size(ref_set,1) / size(data,1);
ref.idx = 0;
ref.set = [];
[U1, ~] = DiffusionMap_ref(data, 2, ref);
%subplot(1,2,1)
%scatter(U1(:, 2), U1(:, 3), 5, c1, 'filled');     axis tight;
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(U1(data_forplot, 2), U1(data_forplot, 3), 30, 'filled'); axis tight;
axis off;
%title('embedding via random reference set', 'fontsize', 25)

%% ref set by pdf
ref.ratio = 0;
ref.idx = 0;
ref.set = ref_set;
[U2, S2] = DiffusionMap_ref(data, 2, ref);
%subplot(1,2,2)
%scatter(U2(:, 2), U2(:, 3), 5, c1, 'filled');     axis tight;
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(U2(data_forplot, 2), U2(data_forplot, 3), 30, 'filled'); axis tight;
axis off;
%title('embedding via reference set design', 'fontsize', 25)

%% check uniformality with linspace S1

theta1 = rand(size(data,1),1)*2*pi; theta1 = sort(theta1);

U = U1(:,2:3);
a1 = mean([max(U(:,1)) min(U(:,1))]);
a2 = mean([max(U(:,2)) min(U(:,2))]);
U = U - [a1 a2]; %centerize
theta_rand = [];
for i = 1:size(data,1)
    theta_rand = [theta_rand; atan_convert(U(i,2), U(i,1))];
end   
theta_rand = sort(theta_rand);

U = U2(:,2:3);
a1 = mean([max(U(:,1)) min(U(:,1))]);
a2 = mean([max(U(:,2)) min(U(:,2))]);
U = U - [a1 a2];
theta_pdf = [];
for i = 1:size(data,1)
    theta_pdf = [theta_pdf; atan_convert(U(i,2), U(i,1))];
end    
theta_pdf = sort(theta_pdf);


%subplot(1,2,1)
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(theta_rand, theta1, 30, 'filled'); axis tight;
hold on; plot(theta1, theta1, '--', 'linewidth', 8); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
%xlabel('theta from embedded S^1', 'fontsize', 25)
%ylabel('theta from uniform S^1', 'fontsize', 25)
%title('embedding via random reference set', 'fontsize', 25)


%subplot(1,2,2)
figure('Renderer', 'painters', 'Position', [10 10 1000 900]);
scatter(theta_pdf, theta1, 30, 'filled'); axis tight;
hold on; plot(theta1, theta1, '--', 'linewidth', 8); axis tight;
xticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks([0.01 pi/2 pi 3*pi/2 2*pi-0.01])
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
ax = gca; ax.FontSize = 35; 
%xlabel('theta from embedded S^1', 'fontsize', 25)
%ylabel('theta from uniform S^1', 'fontsize', 25)
%title('embedding via reference set design', 'fontsize', 25)
