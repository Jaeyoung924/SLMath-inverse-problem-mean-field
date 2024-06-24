% close all;
clear all;
rng(0);

%% Set dimension
Total_N = 500;
N = 150; % number of theta
Total_M = 200;
M = 100; % number of given data
stop_criterion = 10^(-8);
max_iter = 100000;
ht = 0.2; Nt = int64(1/ht); % for U
hs = 1;

%% Given data
Theta_sol = randn(Total_N,1);
% figure; subplot(1,2,1); histogram(Theta_sol,Total_N);
% Total_X = rand(1,Total_M);
Total_X = linspace(-5,5,Total_M);

%% Get reference U using Theta_sol
tic
Total_U_ref = Euler2_U(Total_X,Theta_sol,Total_N,Total_M,ht,Nt); % Euler_U follows o(NM).
toc

Total_U_ref = Total_U_ref(end,:);

Theta = rand(N,1); % Initial data for Theta
idx_per = randperm(Total_M); idx_used = idx_per(1:M);
X = Total_X(idx_used); U_ref = Total_U_ref(idx_used);

%% %% Forward Euler method
s = 0;
L = [];
tic
while (s<max_iter)
    U = Euler2_U(X,Theta,N,M,ht,Nt); % size: (Nt+1)*M
    eta = Euler2_eta(U(end,:),Theta,U_ref,N,M,ht,Nt); % size: (Nt+1)*M
    grd = (hs*ht)/(N*M)*int_sum(U,Theta,eta,N,M,Nt);
    if sum(abs(grd),'all')<stop_criterion
        break;
    end
    Theta = Theta-grd;
    L = [L, sum((U(end,:)-U_ref).^2)];
    s = s + 1;
    if mod(s,1000)==0
        disp(s)
    end
end
toc

%% Plotting
% subplot(1,2,2); histogram(Theta,N);

% figure; subplot(1,2,1); hold on
% scatter(X,U(end,:),'r'); scatter(X,U_ref,'b');

%% Test
% test_X = Total_X(Total_M-49:Total_M);
% test_sol_U = Total_U_ref(Total_M-49:Total_M);
% test_U = Euler2_U(test_X,Theta,N,50,ht,Nt);
% test_U = test_U(end,:);
% 
% % subplot(1,2,2); hold on
% figure; hold on
% scatter(test_X,test_U,'r'); scatter(test_X,test_sol_U,'b');
% title("test, N="+N)
% 
% figure; plot(1:s,L);

%% Plot2 
figure;
subplot(1,2,1); hold on
scatter(X,U(end,:),'r'); scatter(X,U_ref,'b');
subplot(1,2,2);
scatter(X,abs(U(end,:)-U_ref));

L_infty = max(abs(U(end,:)-U_ref));
L_1 = sum(abs(U(end,:)-U_ref));

%% Plot3
figure;
