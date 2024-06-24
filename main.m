close all; clear all;
rng(1.17)

%% Set dimension
Total_N = 400;
N = 100; % number of theta
Total_M = 200;
M = 50; % number of given data
stop_criterion = 0.0001;
h1 = 0.001; Nh = int64(1/h1); % for U
h2 = 01;

%% Given data
Theta_sol = randn(Total_N,1);
figure; subplot(1,2,1); histogram(Theta_sol,Total_N);
Theta = rand(N,1);
Total_X = rand(1,Total_M);

%% Get reference U using Theta_sol
tic
[Total_U_ref,part_U] = Euler_U(Total_X,Theta_sol,Total_N,Total_M,h1,Nh); % Euler_U follows O(NM).
toc

idx_per = randperm(Total_M); idx_used = idx_per(1:M);
X = Total_X(idx_used); U_ref = Total_U_ref(idx_used);

%% Forward Euler method
s = 0;
E = [];
tic
while 1
    [U,part_U] = Euler_U(X,Theta,N,M,h1,Nh);
    E = [E, (1/(2*M))*sum((U-U_ref).^2)];
    grd = (h2/M)*sum((ones(N,1)*(U-U_ref)).*part_U,2);
    if sum(abs(grd))<stop_criterion
        break;
    end
    Theta = Theta-grd;
    s = s + 1;
end
toc

%% Plotting
subplot(1,2,2); histogram(Theta,N);
theta_mean = sum(Theta_sol)/N;

%% Test
test_X = Total_X(Total_M-49:Total_M);
test_sol_U = Total_U_ref(Total_M-49:Total_M);
[test_U,part_U] = Euler_U(test_X,Theta,N,50,h1,Nh);

figure; hold on
scatter(test_X,test_U,'r'); scatter(test_X,test_sol_U,'b');