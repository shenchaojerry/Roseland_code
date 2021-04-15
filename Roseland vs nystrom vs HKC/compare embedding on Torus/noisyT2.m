%% 5000 noisy data points with 200 clean landmarks

ref.size = 200;
N = 5000 + ref.size;
c = linspace(1, N-ref.size, N-ref.size);
%% generate T^2 data
a = 0.8; b = .4;
theta = nan(N, 1);
i = 1;
while i <= N
    xvec = 2 * pi * rand(1);
    yvec = (1 / pi) * rand(1);
    fx = (1 + (3 / 5) * cos(xvec)) / (2 * pi);
    if (yvec < fx)
        theta(i) = xvec;
        i = i + 1;
    else
        continue
    end
end
% by surface of revolution formula
phi = 2 * pi * rand(N, 1);
x = (a + b * cos(theta)) .* cos(phi);
y = (a + b * cos(theta)) .* sin(phi);
z = 1.5 * b * sin(theta);
data = sortrows([x, y, z]);
p = size(data, 2);
q = p;

refind = randperm(N);
refind = refind(1:ref.size);
ref.set = data(refind, :);
data(refind, :) = [];

% add noise
data = [data zeros(size(data,1), 100-size(data,2))];
Noise       = randn(size(data))* .1;
data      = data + Noise ;
ref.set = [ref.set zeros(size(ref.set,1), 100-size(ref.set,2))];

figure('Renderer', 'painters', 'Position', [10 10 1200 700]); axis tight;
scatter3(data(:, 1), data(:, 2), data(:, 3), 18, c, 'filled'); 
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
hold on; scatter3(ref.set(:, 1), ref.set(:, 2), ref.set(:, 3), 90, 'r', 'filled'); axis tight;
%ax = gca; ax.FontSize = 40; 
axis off
%title('Noisy Torus and clean subset','fontsize', 20)

% cut into half for better visulization

U = data;
x = U(:,1);
y = U(:,2);
z = U(:,3);
Ux = U(x>0, :);
Uy = U(y>0, :);
Uz = U(z>0, :);

%subplot(3,3,2)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]); axis tight;
scatter3(Ux(:, 1), Ux(:, 2), Ux(:, 3), 40, c(x>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off
%title('Torus vertical cut','fontsize', 20)

%subplot(3,3,3)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]); axis tight;
scatter3(Uz(:, 1), Uz(:, 2), Uz(:, 3), 40, c(z>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
ax = gca; ax.FontSize = 40; 
%title('Torus horizontal cut','fontsize', 20)


%% DM
[U1, ~] = DiffusionMap(data, 1, 3, 100, 1);
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(U1(:, 1), U1(:, 2), U1(:, 3), 25, c, 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off
%title('DM embedding','fontsize', 25)

U = U1;
x = U(:,1);
y = U(:,2);
z = U(:,3);
Ux = U(x>0, :);
Uy = U(y>-0.02, :);
Uz = U(z>0, :);

figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uy(:, 1), Uy(:, 2), Uy(:, 3), 35, c(y>-0.02), 'filled'); axis tight;
%scatter3(Ux(:, 1), Ux(:, 2), Ux(:, 3), 25, c(x>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off

figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uz(:, 1), Uz(:, 2), Uz(:, 3), 35, c(z>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
ax = gca; ax.FontSize = 40; 

%% Nystrom
sample = ref.set;
sample_size = ref.size; 
[V, ~] = Nystrom(data, sample, 3, 0);
V = V(sample_size+1:end, :);
V = V(:, 2:end);
V = normc(V);
 
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(V(:, 1),V(:, 2), V(:, 3), 25, c, 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off

% cut into half for better visulization
U = V;
x = U(:,1);
y = U(:,2);
z = U(:,3);
Ux = U(x>-.005, :);
Uy = U(y>-.005, :);
Uz = U(z>0, :);

%subplot(3,3,5)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uy(:, 1), Uy(:, 2), Uy(:, 3), 35, c(y>-.005), 'filled'); axis tight;
%scatter3(Ux(:, 1), Ux(:, 2), Ux(:, 3), 35, c(x>-.005), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off
%title('Nystrom vertical cut','fontsize', 20)

%subplot(3,3,6)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uz(:, 1), Uz(:, 2), Uz(:, 3), 35, c(z>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
ax = gca; ax.FontSize = 40; 
%title('Nystrom horizontal cut','fontsize', 20)


%% Roseland
[U2, ~] = refDM(data, 3, ref, 0);
U2 = U2(:, 2:end);
U2 = normc(U2);

%subplot(3,3,7)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(U2(:, 1),U2(:, 2), U2(:, 3), 25, c, 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off

% cut into half for better visulization
U = U2;
x = U(:,1);
y = U(:,2);
z = U(:,3);
Ux = U(x>-.005, :);
Uy = U(y>-.005, :);
Uz = U(z>0, :);

%subplot(3,3,8)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uy(:, 1), Uy(:, 2), Uy(:, 3), 35, c(y>-.005), 'filled'); axis tight;
%scatter3(Ux(:, 1), Ux(:, 2), Ux(:, 3), 35, c(x>-.005), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
%ax = gca; ax.FontSize = 40; 
axis off
%title('ref DM vertical cut','fontsize', 20)

%subplot(3,3,9)
figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(Uz(:, 1), Uz(:, 2), Uz(:, 3), 35, c(z>0), 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
ax = gca; ax.FontSize = 40; 
%title('ref DM horizontal cut','fontsize', 20)    


%% Haddad
u_Haddad = HKC(data, 3, ref, 0);
u_Haddad = u_Haddad(:, 2:end);
u_Haddad = normc(u_Haddad);

figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
scatter3(u_Haddad(:, 1),u_Haddad(:, 2), u_Haddad(:, 3), 25, c, 'filled'); axis tight;
xlabel('X', 'fontsize', 40); ylabel('Y', 'fontsize', 40); zlabel('Z', 'fontsize', 40)
axis off

% cut into half for better visulization
U = u_Haddad;
x = U(:,1);
y = U(:,2);
z = U(:,3);
Ux = U(x>-.006, :);
Uy = U(y>-.006, :);
Uz = U(z>0, :);

figure('Renderer', 'painters', 'Position', [10 10 1200 700]);
%scatter3(Uy(:, 1), Uy(:, 2), Uy(:, 3), 35, c(y>-.006), 'filled'); axis tight;
scatter3(Ux(:, 1), Ux(:, 2), Ux(:, 3), 35, c(x>-.006), 'filled'); axis tight;
axis off
%title('ref DM vertical cut','fontsize', 20)


