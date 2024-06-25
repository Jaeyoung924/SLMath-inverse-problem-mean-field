function grd = int_sum(U,Theta,eta,N,M,Nt)
    red_U = U(2:Nt+1,:); red_eta = eta(2:Nt+1,:);
    grd = zeros(N,1);
    for l = 1 : Nt
        grd = grd + sum(part_theta_f(red_U(l,:),red_eta(l,:),Theta,N,M),2);        
    end
end

function part_theta_f_value = part_theta_f(red_U_fixed_time,red_eta_fixed_time,Theta,N,M)
%     part_theta_f_value = ones(N,1)*red_U_fixed_time;
%     part_theta_f_value = -2*((ones(N,1)*red_U_fixed_time - Theta*ones(1,M))>0);
    part_theta_f_value = -1./( exp(Theta*ones(1,M) - ones(N,1)*red_U_fixed_time)+1 );
%     part_theta_f_value = 2*(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).*exp(-(ones(N,1)*red_U_fixed_time - Theta*ones(1,M)).^2);
    part_theta_f_value = part_theta_f_value.*(ones(N,1)*red_eta_fixed_time);

end