%% speed & embedding of Nystrom, HKC, Roseland, beta=0.5
n_t	= 2^7;   %data dim
t	= linspace(-1.5,1.5,n_t);
T_nys05 = [];
U_nys05 = cell(0,0);
T_ref05 = [];
U_ref05 = cell(0,0);
T_hkc05 = [];
U_hkc05 = cell(0,0);
ref.idx = 0;
beta = 0.5;

count = 0;
for N = [5000 10000 20000 40000 80000 160000 320000 640000]  %data size
    count = count + 1;
    
    T_nys_temp = [];
    T_ref_temp = [];
    T_hkc_temp = [];
    
    for K = 1:20 %run many times
        % get data and subset
        theta	= rand(N,1)*2*pi;
        data = shepp_logan_proj_2(theta,t)';
        ref.size = round(N.^beta);
        theta	= rand(ref.size,1)*2*pi;
        ref.set = shepp_logan_proj_2(theta,t)';
        %
        %nystrom
        tic;
        [u_nys, ~] = Nystrom(data, ref.set, 3);
        u_nys = u_nys(ref.size+1:end, 2:end);
        t1 = toc;
        T_nys_temp = [T_nys_temp t1];
        %}
        %
        %Roseland
        tic;
        [u_ref, ~]  = roseland(data, 3, ref, 0);
        u_ref = u_ref(:,2:end);
        t1 = toc;
        T_ref_temp = [T_ref_temp t1];
        %}
        %
        %HKC
        tic;
        u_hkc05 = HKC(data, 3, ref, 0);
        u_hkc05 = u_hkc05(:,2:end);
        t1 = toc;
        T_hkc_temp = [T_hkc_temp t1];
        %}
    end
    
    T_nys05 = [T_nys05; mean(T_nys_temp)];
    U_nys05{1,count} = u_nys;
    T_ref05 = [T_ref05; mean(T_ref_temp)];
    U_ref05{1,count} = u_ref;
    T_hkc05 = [T_hkc05; mean(T_hkc_temp)];
    U_hkc05{1,count} = u_hkc05;
end
