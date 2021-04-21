function W=ReturnW_frob_gap(C, i, j, q)

%get matrix of correlations just for sites i and j
W = zeros(q, q); 

%couplings in 1st line & col, corresponding to gap (residue type 1), are set to ZERO (this
%type of residue is the one for which correlations have been suppressed)
W(2:q, 2:q) = C(mapkey(i, 1:q-1, q), mapkey(j, 1:q-1, q));

end
