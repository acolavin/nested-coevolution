function W=ReturnW(C,i,j,q)
% extracts coupling matrix for columns i and j

    W = ones(q,q);
    W(1:q-1,1:q-1) = exp( -C(mapkey(i,1:q-1,q),mapkey(j,1:q-1,q)) );

end