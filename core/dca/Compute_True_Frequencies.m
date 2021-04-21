

function [Pij_true,Pi_true,Meff] = Compute_True_Frequencies(align,M,N,q,theta)
% computes reweighted frequency counts

    W = ones(1,M);
    if( theta > 0.0 )
        W = (1./(1+sum(squareform(pdist(align,'hamm')<theta))));
    end
    Meff=sum(W);

    Pij_true = zeros(N,N,q,q);
    Pi_true = zeros(N,q);

    for j=1:M
        for i=1:N
            Pi_true(i,align(j,i)) = Pi_true(i,align(j,i)) + W(j);
        end
    end
    Pi_true = Pi_true/Meff;

    for l=1:M
        for i=1:N-1
            for j=i+1:N
                Pij_true(i,j,align(l,i),align(l,j)) = Pij_true(i,j,align(l,i),align(l,j)) + W(l);
                Pij_true(j,i,align(l,j),align(l,i)) = Pij_true(i,j,align(l,i),align(l,j));
            end
        end
    end
    Pij_true = Pij_true/Meff;

    scra = eye(q,q);
    for i=1:N
        for alpha=1:q
            for beta=1:q
                Pij_true(i,i,alpha,beta) = Pi_true(i,alpha) * scra(alpha,beta);
            end
        end
    end
end
