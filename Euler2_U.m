function U = Euler2_U(X,Theta,N,M,ht,Nt,f_id)
    U = zeros(Nt+1,M); U(1,:) = X;
    for i = 1 : Nt
        U(i+1,:) = U(i,:) + (ht/N)*sum(f(U(i,:),Theta,N,M,f_id),1);
    end
end

function f_value = f(U,Theta,N,M,f_id)
    switch f_id
        case 1
            f_value = (((ones(N,1)*U - Theta*ones(1,M))>0).*(ones(N,1)*U - Theta*ones(1,M))).^2;
        case 2
            f_value = log( 1+exp(ones(N,1)*U - Theta*ones(1,M)) );
        case 3
            f_value = exp(-(ones(N,1)*U - Theta*ones(1,M)).^2);
        case 4
            f_value = (log( 1+exp(ones(N,1)*U - Theta*ones(1,M)) )).^2;
    end
end