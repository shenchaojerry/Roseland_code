%% speed & embedding of DM

KNN = 100;
n_t	= 2^7;   %data dim
t	= linspace(-1.5,1.5,n_t);
T1 = [];
U1 = cell(0,0);

count = 0;
for N = [5000 10000 20000 40000 80000]  %data size
    count = count + 1;
    theta	= rand(N,1)*2*pi;
    theta	= sort(theta);
    data = shepp_logan_proj_2(theta,t)';

    %DM
    tic ;
    [u1, ~] = DiffusionMap1(data, 3, KNN, 0);
    t1 = toc;
    T1 = [T1; t1];
    U1{1,count} = u1;
end
