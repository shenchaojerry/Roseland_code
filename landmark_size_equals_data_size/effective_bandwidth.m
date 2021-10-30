%% effective bandwith

h_ref = 0.01;

%% S1 data
beta = 1;
ref_size = round(2500^beta);
N = 2500 + ref_size;
theta = rand(N,1)*2*pi; 
theta = sort(theta);
data = [cos(theta) sin(theta)];
refind = randperm(N);
refind = refind(1:ref_size);
ref_set = data(refind, :);
data(refind, :) = [];
theta(refind) = [];

% affinity via ref
dist = pdist2(data, ref_set);
W_ref = exp( - dist.^2 / h_ref );
aff_ref = W_ref(1000, :) * W_ref';
aff_ref = aff_ref / ref_size;
C = max(aff_ref);

%% solve for effective bandwidth for dm
dist2 = pdist2(data, data);
d_v = dist2(1000,:).^2;
aff = aff_ref / C;
h_dm = - d_v ./ log(aff);

%%
h_dm = median(h_dm);
dist3 = pdist2(data(1000,:), data);
aff_dm = exp( - dist3.^2 / h_dm);
aff_dm = C * aff_dm;

%% 
hold on
plot(aff_ref, 'linewidth', 2)
plot(aff_dm, 'linewidth', 2)