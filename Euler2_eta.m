function eta = Euler2_eta(U_1,Theta,U_ref,N,M,ht,Nt)
    til_eta = zeros(Nt+1,M); til_eta(1,:) = U_1-U_ref;
    for l = 1 : Nt
        til_eta(l+1,:) = til_eta(l,:) - (ht/N)*sum(part_u_f(til_eta(Nt-l+2,:),Theta,N,M),1);
    end
    eta = flipud(til_eta);
end

function part_u_f_value = part_u_f(U,Theta,N,M)
%     part_u_f_value = Theta*ones(1,M);
%     part_u_f_value = 2*((ones(N,1)*U - Theta*ones(1,M))>0);
    part_u_f_value = 1./( exp(Theta*ones(1,M)-ones(N,1)*U)+1 );
%     part_u_f_value = -2*(ones(N,1)*U-Theta*ones(1,M)).*exp(-(ones(N,1)*U - Theta*ones(1,M)).^2);
end