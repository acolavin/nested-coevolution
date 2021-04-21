function DI = bp_link(i,j,W,P1,q)
% computes direct information

    [mu1, mu2] = compute_mu(i,j,W,P1,q);
    DI = compute_di(i,j,W, mu1,mu2,P1);
    return;
end