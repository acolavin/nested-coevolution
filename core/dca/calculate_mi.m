function [M,s1,s2] = calculate_mi(i,j,P2,P1,q)
% computes mutual information between columns i and j

    M = 0.;
    for alpha=1:q
        for beta = 1:q
             if( P2(i,j,alpha,beta)>0 )
                M = M + P2(i,j,alpha, beta)*log(P2(i,j, alpha, beta) / P1(i,alpha)/P1(j,beta));
            end
        end
    end

    s1=0.;
    s2=0.;
    for alpha=1:q
        if( P1(i,alpha)>0 )
            s1 = s1 - P1(i,alpha) * log(P1(i,alpha));
        end
        if( P1(j,alpha)>0 )
            s2 = s2 - P1(j,alpha) * log(P1(j,alpha));
        end
    end

end