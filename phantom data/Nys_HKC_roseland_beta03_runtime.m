%% speed & embedding of Nystrom, HKC, Roseland, beta=0.3
n_t	= 2^7;   %data dim
t	= linspace(-1.5,1.5,n_t);
T_nys03 = [];
U_nys03 = cell(0,0);
T_ref03 = [];
U_ref03 = cell(0,0);
T_hkc03 = [];
U_hkc03 = cell(0,0);
ref.idx = 0;
beta = 0.3;

count = 0;
for N = [5000 10000 20000 40000 80000 160000 320000 640000 1280000]  %data size
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

        %nystrom
        tic;
        [u_nys, ~] = Nystrom(data, ref.set, 3);
        u_nys = u_nys(ref.size+1:end, 2:end);
        t1 = toc;
        T_nys_temp = [T_nys_temp t1];
        
        %roseland
        tic;
        [u_ref, ~]  = roseland(data, 3, ref, 0);
        u_ref = u_ref(:,2:end);
        t1 = toc;
        T_ref_temp = [T_ref_temp t1];
        
        %HKC
        tic;
        u_hkc03 = HKC(data, 3, ref, 0);
        u_hkc03 = u_hkc03(:,2:end);
        t1 = toc;
        T_hkc_temp = [T_hkc_temp t1];

    end
    
    T_nys03 = [T_nys03; mean(T_nys_temp)];
    T_ref03 = [T_ref03; mean(T_ref_temp)];
    T_hkc03 = [T_hkc03; mean(T_hkc_temp)];
    U_nys03{1,count} = u_nys;
    U_ref03{1,count} = u_ref;
    U_hkc03{1,count} = u_hkc03;
end
