function grd = int_sum(U,Theta,eta,N,M,Nt,f_id)
    red_U = U(2:Nt+1,:); red_eta = eta(2:Nt+1,:);
    grd = zeros(N,1);
    for l = 1 : Nt
        grd = grd + sum(part_theta_f(red_U(l,:),red_eta(l,:),Theta,N,M,f_id),2);        
    end
end

function part_theta_f_value = part_theta_f(red_U_fixed_time,red_eta_fixed_time,Theta,N,M,f_id)
    switch f_id
        case 1
            part_theta_f_value = -2*(((ones(N,1)*red_U_fixed_time - Theta*ones(1,M))>0).*(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)));
        case 2
            part_theta_f_value = -1./( exp(Theta*ones(1,M) - ones(N,1)*red_U_fixed_time)+1 );
        case 3
            part_theta_f_value = 2*(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).*exp(-(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).^2);
        case 4
            part_theta_f_value = -2*(log( 1+exp(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)) ))./( exp(Theta*ones(1,M) - ones(N,1)*red_U_fixed_time)+1 );
        case 5
            part_theta_f_value = -2*(1./(1+((ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).^2))).*((ones(N,1)*red_U_fixed_time - Theta*ones(1,M)));
        case 6
            part_theta_f_value = 2*(ones(N,1)*red_U_fixed_time - Theta*ones(1,M))./(1+(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).^2).^2;
    end
    part_theta_f_value = part_theta_f_value.*(ones(N,1)*red_eta_fixed_time);

end