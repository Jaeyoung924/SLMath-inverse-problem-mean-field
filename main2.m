% close all;
clear all;
rng(0);

%% Set dimension
Total_N = 3000;
N = 25; % number of theta
Total_M = 1000;
M = 500; % number of given data
stop_criterion = 10^(-8);
max_iter = 1e3;
ht = 0.2; Nt = int64(1/ht); % for U
hs = 100;
n_tests = 30;

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
    case 5
        disp('f(u,\theta) = log(1+(u-\theta).^2)')
    case 6
        disp('f(u,\theta) = 1/(1+(u-\theta).^2)')
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

% --- Uncomment the lines below to run multiple tests ---
Expected_Theta = 0;
Expected_L = 0;
Expected_mean = 0;
Expected_std = 0;
Expected_traj = 0;
% --- Uncomment the lines above to run multiple tests ---

% --- Uncomment the lines below to run multiple tests ---
for tests = 1:n_tests
disp("Test "+tests)
rng(tests);
% --- Uncomment the lines aboves to run multiple tests ---

Theta = rand(N,1); % Initial data for Theta
Theta = sort(Theta);
% Theta = -1 + randn(N,1); % Initial data for Theta
idx_per = randperm(Total_M); idx_used = idx_per(1:M);
X = Total_X(idx_used); U_ref = Total_U_ref(idx_used);

%% %% Forward Euler method
s = 0;
L = []; mean_theta = []; std_theta = [];
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
    std_theta = [std_theta, std(Theta)];
    L = [L, sum((U(end,:)-U_ref).^2)/M];
%     if size(L,2)~=1 && L(end-1)-L(end)<0
%         hs = hs*0.9;
%         break;
%     end
    if mod(s,max_iter/10)==0
        disp(s)
        Theta_traj = [Theta_traj,Theta];
        s_traj = [s_traj,s];
    end
end
%Compute U and L from final Theta
U = Euler2_U(X,Theta,N,M,ht,Nt,f_id); % size: (Nt+1)*M
L = [L, sum((U(end,:)-U_ref).^2)/M];
toc

% --- Uncomment the lines below to run multiple tests ---
Expected_Theta = Expected_Theta + Theta;
Expected_L = Expected_L + L;
Expected_mean = Expected_mean + mean_theta;
Expected_std = Expected_std + std_theta;
Expected_traj = Expected_traj + Theta_traj;
end %for loop over tests
% --- Uncomment the lines above to run multiple tests ---

% --- Uncomment the lines below to run multiple tests ---
Theta = Expected_Theta / n_tests;
U = Euler2_U(X,Theta,N,M,ht,Nt,f_id); % size: (Nt+1)*M
L = Expected_L / n_tests;
mean_theta = Expected_mean / n_tests;
std_theta = Expected_std / n_tests;
Theta_traj = Expected_traj / n_tests;
% --- Uncomment the lines above to run multiple tests ---

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

figure; subplot(1,2,1); semilogy(0:s,L); title('L');
subplot(1,2,2); hold on
plot(1:s,mean_theta,'r'); plot(1:s,std_theta,'b');
legend('mean of theta','variance of theta');

%% Plot2 
figure;
subplot(1,2,1); hold on
scatter(X,U(end,:),'r'); scatter(X,U_ref,'b');
subplot(1,2,2);
scatter(X,abs(U(end,:)-U_ref));

xlin = linspace(-3.5,3.5,1e3);
figure; subplot(1,2,1); histogram(Theta_sol,'Normalization','pdf');
hold on
plot(xlin,1/sqrt(2*pi)*exp(-1/2*xlin.^2),'r-','Linewidth',2)
subplot(1,2,2); histogram(Theta,'Normalization','pdf');
hold on
plot(xlin,1/sqrt(2*pi)*exp(-1/2*xlin.^2),'r-','Linewidth',2)

%Tracjetory plot
[i_end,~] = size(Theta_traj); 
figure();
title("$\theta_i(s)$ Trajectories for $N = $"+N+", $M = $"+M+", $h_s = $"+hs,'Interpreter','Latex')
hold on
for i = 1:i_end
    plot(Theta_traj(i,:),s_traj,'k--o','Linewidth',2)
end
xlabel('$\theta_i(s)$','Interpreter','Latex')
ylabel('$s$','Interpreter','Latex')
grid on
set(gca,'fontsize',12)
set(gca,'linewidth',2)

L_infty = max(abs(U(end,:)-U_ref));
L_1 = norm(U(end,:)-U_ref,1)/M;
L_2 = norm(U(end,:)-U_ref,2)/sqrt(M);
disp(['      L_1 err', '      L_2 err', '    L_inf err'])
disp([L_1,L_2,L_infty])
disp(['  Target mean', '    Est. mean', '   Target std', '     Est. std'])
disp([mean(Theta_sol),mean(Theta),std(Theta_sol),std(Theta)])
disp([' Rel err mean', '  Rel err std'])
disp(abs([mean(Theta)-mean(Theta_sol),std(Theta)-std(Theta_sol)])./abs([mean(Theta_sol),std(Theta_sol)]))
