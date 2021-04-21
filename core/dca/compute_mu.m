function [mu1,mu2] = compute_mu(i,j,W,P1,q)

    epsilon=1e-4;
    diff =1.0;
    mu1 = ones(1,q)/q;
    mu2 = ones(1,q)/q;
    pi = P1(i,:);
    pj = P1(j,:);

    while ( diff > epsilon )

        scra1 = mu2 * W';
        scra2 = mu1 * W;

        new1 = pi./scra1;
        new1 = new1/sum(new1);

        new2 = pj./scra2;
        new2 = new2/sum(new2);

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) );

        mu1 = new1;
        mu2 = new2;

    end
end