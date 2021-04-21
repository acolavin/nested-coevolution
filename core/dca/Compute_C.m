function C = Compute_C(Pij,Pi,N,q)
% computes correlation matrix

    C=zeros(N*(q-1),N*(q-1));
    for i=1:N
        for j=1:N
            for alpha=1:q-1
                for beta=1:q-1
                     C(mapkey(i,alpha,q),mapkey(j,beta,q)) = Pij(i,j,alpha,beta) - Pi(i,alpha)*Pi(j,beta);
                end
            end
        end
    end
end