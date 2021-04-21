
function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end