function DI = compute_di(i,j,W, mu1,mu2, Pia)
% computes direct information

    tiny = 1.0e-100;

    Pdir = W.*(mu1'*mu2);
    Pdir = Pdir / sum(sum(Pdir));

    Pfac = Pia(i,:)' * Pia(j,:);

    DI = trace( Pdir' * log( (Pdir+tiny)./(Pfac+tiny) ) );

end
