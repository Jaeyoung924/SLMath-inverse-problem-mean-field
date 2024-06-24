function [U,part_U] = Euler_U(X,Theta,N,M,h1,Nh)
    U = X; % initial data for U, size: 1*M
    part_U = zeros(N,M); % initial data for partial derivative of U, size: N*M
    for i = 1 : Nh
        U = U + (h1/N)*ft_f(U,Theta,N,M);
        part_U = part_U + (h1/N)*(ones(N,1)*deri_uf(U,Theta,N,M)).*part_U + (h1/N)*deri_thetaf(U,Theta,N,M);
    end
end

%% function f and its partial derivative
function f = ft_f(U,Theta,N,M)
%     f = sin(ones(N,1)*U + Theta*ones(1,M)); % size: N*M
%     f = (ones(N,1)*U).*(Theta*ones(1,M)); % to test
    f = ((ones(N,1)*U - Theta*ones(1,M))>0).^2;
    f = sum(f,1); % size: 1*M
end

function deri_uf = deri_uf(U,Theta,N,M)
%     deri_uf = cos(ones(N,1)*U + Theta*ones(1,M)); % size: N*M
%     deri_uf = Theta*ones(1,M); % to test
    deri_uf = 2*((ones(N,1)*U - Theta*ones(1,M))>0);
    deri_uf = sum(deri_uf,1); % size: 1*M
end

function deri_thetaf = deri_thetaf(U,Theta,N,M)
%     deri_thetaf = cos(ones(N,1)*U + Theta*ones(1,M));
%     deri_thetaf = ones(N,1)*U; % to test
    deri_thetaf = -2*((ones(N,1)*U - Theta*ones(1,M))>0);
end