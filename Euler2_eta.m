function eta = Euler2_eta(U_1,Theta,U_ref,N,M,ht,Nt,f_id)
    til_eta = zeros(Nt+1,M); til_eta(1,:) = U_1-U_ref;
    for l = 1 : Nt
        til_eta(l+1,:) = til_eta(l,:) - (ht/N)*sum(part_u_f(til_eta(Nt-l+2,:),Theta,N,M,f_id),1);
    end
    eta = flipud(til_eta);
end

function part_u_f_value = part_u_f(U,Theta,N,M,f_id)
    switch f_id
        case 1
            part_u_f_value = 2*(((ones(N,1)*U - Theta*ones(1,M))>0).*(ones(N,1)*U - Theta*ones(1,M)));
        case 2
            part_u_f_value = 1./( exp(Theta*ones(1,M)-ones(N,1)*U)+1 );
        case 3
            part_u_f_value = -2*(ones(N,1)*U-Theta*ones(1,M)).*exp(-(ones(N,1)*U - Theta*ones(1,M)).^2);
        case 4
            part_u_f_value = 2*(log( 1+exp(ones(N,1)*U - Theta*ones(1,M)) ))./( exp(Theta*ones(1,M)-ones(N,1)*U)+1 );
        case 5
            part_u_f_value = 2*(1./(1+((ones(N,1)*U - Theta*ones(1,M)).^2))).*((ones(N,1)*U - Theta*ones(1,M))) ;
        case 6
            part_u_f_value = -2.*(ones(N,1)*U - Theta*ones(1,M))./(1+(ones(N,1)*U - Theta*ones(1,M)).^2).^2;
    end
end