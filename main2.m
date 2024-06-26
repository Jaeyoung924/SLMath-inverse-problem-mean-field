% close all;
clear all;
rng(0);

%% Set dimension
Total_N = 3000;
N = 200; % number of theta
Total_M = 500;
M = 300; % number of given data
stop_criterion = 10^(-8);
max_iter = 100000;
ht = 0.2; Nt = int64(1/ht); % for U
hs = 2;

disp(['     Total N','     Total M','           N','           M','         h_s'])
disp([Total_N, Total_M, N, M, hs])

%% Set function f(u(t,x),\theta)
f_id = 3; 
switch f_id
    case 1
        disp('f(u,\theta) = ((u-\theta)^+).^2')
    case 2
        disp('f(u,\theta) = log(1+exp(u-\theta))')
    case 3
        disp('f(u,\theta) = exp(-(u-\theta).^2)')
    case 4
        disp('f(u,\theta) = log(1+exp(u-\theta)).^2')
end

%% Given data
Theta_sol = randn(Total_N,1);
% Total_X = rand(1,Total_M);
Total_X = linspace(-5,5,Total_M);

%% Get reference U using Theta_sol
tic
Total_U_ref = Euler2_U(Total_X,Theta_sol,Total_N,Total_M,ht,Nt,f_id); % Euler_U follows o(NM).
toc

Total_U_ref = Total_U_ref(end,:);

Theta = 10*rand(N,1); % Initial data for Theta
idx_per = randperm(Total_M); idx_used = idx_per(1:M);
X = Total_X(idx_used); U_ref = Total_U_ref(idx_used);

%% %% Forward Euler method
s = 0;
L = []; mean_theta = []; var_theta = [];
s_traj = s;
Theta_traj = Theta;
tic
while (s<max_iter)
    U = Euler2_U(X,Theta,N,M,ht,Nt,f_id); % size: (Nt+1)*M
    eta = Euler2_eta(U(end,:),Theta,U_ref,N,M,ht,Nt,f_id); % size: (Nt+1)*M
    grd = (hs*ht)/(N*M)*int_sum(U,Theta,eta,N,M,Nt,f_id);
%     if sum(abs(grd),'all')<stop_criterion
%         break;
%     end
    Theta = Theta-grd; s = s+1;
    mean_theta = [mean_theta, sum(Theta)/N];
    var_theta = [var_theta, var(Theta)];
    L = [L, sum((U(end,:)-U_ref).^2)/M];
    if size(L,2)~=1 && L(end-1)-L(end)<0
%         hs = hs*0.9;
        break;
    end
    if mod(s,1000)==0
        disp(s)
        Theta_traj = [Theta_traj,Theta];
        s_traj = [s_traj,s];
    end
end
%Compute U and L from final Theta
U = Euler2_U(X,Theta,N,M,ht,Nt,f_id); % size: (Nt+1)*M
L = [L, sum((U(end,:)-U_ref).^2)/M];
toc

%% Plotting
% figure; subplot(1,2,1); hold on
% scatter(X,U(end,:),'r'); scatter(X,U_ref,'b');

%% Test
test_X = Total_X(idx_per(end-49:end));
test_sol_U = Total_U_ref(idx_per(end-49:end));
test_U = Euler2_U(test_X,Theta,N,50,ht,Nt,f_id);
test_U = test_U(end,:);

% subplot(1,2,2); hold on
figure; hold on
scatter(test_X,abs(test_U-test_sol_U),'r');
title("test, N="+N)

figure; subplot(1,2,1); plot(1:s,L); title('L');
subplot(1,2,2); hold on
plot(1:s,mean_theta,'r'); plot(1:s,var_theta,'b');
legend('mean of theta','variance of theta');

%% Plot2 
figure;
subplot(1,2,1); hold on
scatter(X,U(end,:),'r'); scatter(X,U_ref,'b');
subplot(1,2,2);
scatter(X,abs(U(end,:)-U_ref));

figure; subplot(1,2,1); histogram(Theta_sol,Total_N);
subplot(1,2,2); histogram(Theta,N);

L_infty = max(abs(U(end,:)-U_ref));
L_1 = norm(U(end,:)-U_ref,1)/M;
L_2 = norm(U(end,:)-U_ref,2)/sqrt(M);
disp(['      L_1 err', '      L_2 err', '    L_inf err'])
disp([L_1,L_2,L_infty])
disp(['  Target mean', '    Est. mean', '   Target std', '     Est. std'])
disp([mean(Theta_sol),mean(Theta),std(Theta_sol),std(Theta)])
disp([' Abs err mean', '  Abs err std'])
disp(abs([mean(Theta)-mean(Theta_sol),std(Theta)-std(Theta_sol)])./abs([mean(Theta_sol),std(Theta_sol)]))
